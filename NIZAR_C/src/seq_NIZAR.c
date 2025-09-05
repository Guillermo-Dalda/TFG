#include "seq_NIZAR.h"

void diversity_phase(Mersenne_Twister* rng, double* X_p, double* x_best, double* x_i, double* x_j, double* x_k, double* x_m, double* P1, double* P2, double* P3, double* T2, double* T3, double* Vj, double* Vk, double*  aux1, double*  aux2, double* alpha, double beta1, double beta2, double deltaj, double deltak, int* lambda, int dim){
    //P1
    construct_P1(P1, x_best, x_m, aux1, aux2, beta1, lambda, dim);
    if (lambda[0] == 1){
        //S1
        construct_S1(X_p, P1, x_i, x_j, x_k, aux1, aux2, beta1, beta2, lambda, dim);
    }
    else {
        //T3
        construct_T3(T3, x_best, x_i, lambda, dim);
        //Vj y Vk
        transfer_map(x_j, T3, alpha[0], Vj, dim);
        transfer_map(x_k, T3, alpha[0], Vk, dim);
        //T2
        construct_T2(rng, T2, T3, P1, alpha, beta2, lambda, dim);
        //P2
        construct_P2(P2, T2, x_m, aux1, beta1, lambda, dim);
        //P3
        construct_P3(rng, P3, P1, T2, T3, alpha, beta2, lambda, dim);
        //S2
        construct_S2(X_p, P2, P3, Vj, Vk, aux1, aux2, deltaj, deltak, lambda, dim);
    }
}

void overlap_phase(Mersenne_Twister* rng, double* X_p, double* x_best, double* x_i, double* x_j, double* x_k, double* x_m, double* P1, double* aux1, double* aux2, double* B1,  double* B2, int* lambda, int dim){
    if(vector_comparison(X_p, x_best, dim) || vector_comparison(x_i, x_best, dim) || (getRandomNumber(rng) / (double)UINT_MAX) <= 0.25){
        for (int i = 0; i < dim; i++){
            B1[i] = (getRandomNumber(rng) / (double)UINT_MAX) * 2 - 1;
            B2[i] = (getRandomNumber(rng) / (double)UINT_MAX) * 2 - 1;
        }
        reconstruct_P1(P1, x_best, x_m, aux1, aux2, B1, lambda, dim);
        reconstruct_S1(X_p, P1, x_i, x_j, x_k, aux1, aux2, B1, B2, lambda, dim);
    }
}

void check_bounds(double* X_p, double* x_i, double* lower_bound, double* upper_bound, int dim){
    for (int i = 0; i < dim; i++){
        if(X_p[i] > upper_bound[i] || X_p[i] < lower_bound[i])
            X_p[i] = x_i[i];
    }
}

void generate_alphas(Mersenne_Twister* rng, double* alpha){
    alpha[0] = getRandomNumber(rng) / (double)UINT_MAX;
    alpha[1] = getRandomNumber(rng) / (double)UINT_MAX - 3 / 8;
    alpha[2] = getRandomNumber(rng) / (double)UINT_MAX - 1 / 4;
    alpha[3] = getRandomNumber(rng) / (double)UINT_MAX - 3 / 8;
}

void generate_lambdas(Mersenne_Twister* rng, int* lambda){
    for(int i = 0; i < 8; i++){
        lambda[i] = getRandomNumber(rng) % 2;
    }
}

void get_JKM_betas_deltas(Mersenne_Twister* rng, int n, int i, int* j, int* k, int* m, double* beta1, double* beta2, double* deltaj, double* deltak){
    *j = getRandomNumber(rng) % n;
    *k = getRandomNumber(rng) % n;
    *m = getRandomNumber(rng) % n;
    while (i == *j)
        *j = getRandomNumber(rng) % n;
    while (i == *k || *j == *k)
        *k = getRandomNumber(rng) % n;
    while (i == *m | *j == *m || *k == *m)
        *m = getRandomNumber(rng) % n;
    
    *beta1 = getRandomNumber(rng) / (double)UINT_MAX;
    *beta2 = getRandomNumber(rng) / (double)UINT_MAX;
    *deltaj = *beta1 * ((*j % 2 == 0) ? 1 : -1);
    *deltak = *beta2 * ((*k % 2 == 0) ? 1 : -1);
}

int current_best_solution(int n, int dim, double* solutions){
    int best = 0;
    double min = solutions[0];
    
    for (int i = 1; i < n; i++){
        if (solutions[i] < min){
            min = solutions[i];
            best = i;
        }
    }
    return best;
}

void evaluate_fitness(int dim, double* x, double* solution, double (*objective_function)(const double*, const int)){
    *solution = (*objective_function)(x, dim);
}

double* NOA(Mersenne_Twister* rng, int n, int stopping_criteria, int dim, double* lower_bound, double* upper_bound, double (*objective_function)(const double*, const int)){
    //initialize pointers
    double** x = (double**)malloc(n * sizeof(double*));
    x[0] = (double*)malloc(n * dim * sizeof(double));
    double* solutions = malloc(n * sizeof(double*));
    double X_p[dim];
    double P1[dim];
    double P2[dim];
    double P3[dim];
    double T2[dim];
    double T3[dim];
    double Vj[dim];
    double Vk[dim];
    double B1[dim];
    double B2[dim];
    double aux1[dim];
    double aux2[dim];
    double alpha[4];
    int lambda[8];
    
    //initial population
    for (int i = 0; i < n; i++){
        x[i] = x[0] + (i * dim);
        
        for (int j = 0; j < dim; j++){
            x[i][j] = lower_bound[j] + (getRandomNumber(rng) / (double)UINT_MAX) * (upper_bound[j] - lower_bound[j]);
        }
        //fitness value
        evaluate_fitness(dim, x[i], &solutions[i], objective_function);
    }
    
    int best = 0;
    int current = 0;
    //main loop
    while (current < stopping_criteria){
        generate_alphas(rng, alpha);
        generate_lambdas(rng, lambda);
        best = current_best_solution(n, dim, solutions);

        for (int i = 0; i < n; i++){
            //inicializacion de randoms
            int j,k,m;
            double beta1, beta2, deltaj, deltak;
            get_JKM_betas_deltas(rng, n, i, &j, &k, &m, &beta1, &beta2, &deltaj, &deltak);
            //Diversity phase
            diversity_phase(rng, X_p, x[best], x[i], x[j], x[k], x[m], P1, P2, P3, T2, T3, Vj, Vk, aux1, aux2, alpha, beta1, beta2, deltaj, deltak, lambda, dim);
            //Overlap phase
            overlap_phase(rng, X_p, x[best], x[i], x[j], x[k], x[m], P1, aux1, aux2, B1, B2, lambda, dim);//Eliminar argumentos innecesarios
            //Make sure there are no elements that exceed the function bounds
            check_bounds(X_p, x[i], lower_bound, upper_bound, dim);
            //Evaluate new solution
            double result_aux;
            evaluate_fitness(dim, X_p, &result_aux, objective_function);
            //Save the better solution
            if (result_aux < solutions[i]){
                memcpy(x[i], X_p, dim * sizeof(double));
                solutions[i] = result_aux;
            }
        }
        current+=n;
    }
    best = current_best_solution(n, dim, solutions);
    double* result = malloc((dim + 1) * sizeof(double*));
    memcpy(result, x[best], dim * sizeof(double));
    result[dim] = solutions[best];
    free(x[0]);
    free(x);
    free(solutions);
    
    return result;
}