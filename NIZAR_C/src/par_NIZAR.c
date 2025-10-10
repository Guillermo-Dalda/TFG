#include "par_NIZAR.h"

#ifdef MPI
    #define MASTER 0
    #define NORMAL_TAG 1
#endif

typedef struct params{
    Mersenne_Twister rng;
	int id, n_threads, n, dim, stopping_criteria, rank, num_procs;
	double* lower_bound;
    double* upper_bound;
    double* ini;
    double** x;
    double* solutions;
    bool* update;
    double (*objective_function)(const double*, const int);
    double* alpha;
    int* lambda;
} task;

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

void get_JKM_betas_deltas(Mersenne_Twister* rng, int n, int i, int* j, int* k, int* m, double* beta1, double* beta2, double* deltaj, double* deltak){
    *j = getRandomNumber(rng) % n;
    *k = getRandomNumber(rng) % n;
    *m = getRandomNumber(rng) % n;
    while (i == *j)
        *j = getRandomNumber(rng) % n;
    while (i == *k || *j == *k)
        *k = getRandomNumber(rng) % n;
    while (i == *m || *j == *m || *k == *m)
        *m = getRandomNumber(rng) % n;
    
    *beta1 = getRandomNumber(rng) / (double)UINT_MAX;
    *beta2 = getRandomNumber(rng) / (double)UINT_MAX;
    *deltaj = *beta1 * ((*j % 2 == 0) ? 1 : -1);
    *deltak = *beta2 * ((*k % 2 == 0) ? 1 : -1);
}

void evaluate_fitness(int dim, double* x, double* solution, double (*objective_function)(const double*, const int)){
    *solution = (*objective_function)(x, dim);
}


#ifdef MPI
    void* parallel_task(void* param){
        task tasks = *((task*) param);
        Mersenne_Twister rng = tasks.rng;
        int id = tasks.id;
        int n_threads = tasks.n_threads;
        int n = tasks.n;
        int dim = tasks.dim;
        int rank = tasks.rank;
        int num_procs = tasks.num_procs;
        int stopping_criteria = tasks.stopping_criteria;
        double** shared_x = tasks.x;
        double* lower_bound = tasks.lower_bound;
        double* upper_bound = tasks.upper_bound;
        double* solutions = tasks.solutions;
        bool* update = tasks.update;
        double* alpha = tasks.alpha;
        int* lambda = tasks.lambda;
        pthread_mutex_lock(&mutex);
        unsigned int seed = getRandomNumber(&rng);
        pthread_mutex_unlock(&mutex);
        Create_MersenneTwister(&rng, seed);
        
        int aux = n / n_threads;
        int offset = n % n_threads;
        int start = id * aux + (id < offset ? id : offset);
        int end = start + aux + (id < offset);
        //initialize population
        for(int i = start; i < end; i++){
            shared_x[i] = tasks.ini + (i * dim);
            
            if (rank == MASTER){
                for (int j = 0; j < dim; j++){
                    shared_x[i][j] = lower_bound[j] + (getRandomNumber(&rng) / (double)UINT_MAX) * (upper_bound[j] - lower_bound[j]);
                }
                //fitness value
                evaluate_fitness(dim, shared_x[i], &(solutions[i]), tasks.objective_function);
            }
        }
        pthread_barrier_wait(&barrier);
        if (rank == MASTER && id == 0){
            for (int i = 1; i < num_procs; i++){
                MPI_Send(shared_x[0], n * dim, MPI_DOUBLE, i, NORMAL_TAG, MPI_COMM_WORLD);
                MPI_Send(&(solutions[0]), n, MPI_DOUBLE, i, NORMAL_TAG, MPI_COMM_WORLD);
            }
        }
        else if (id == 0){
            MPI_Recv(shared_x[0], n * dim, MPI_DOUBLE, MASTER, NORMAL_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(solutions, n, MPI_DOUBLE, MASTER, NORMAL_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        pthread_barrier_wait(&barrier);
        double** x = (double**)malloc(n * sizeof(double*));
        x[0] = (double*)malloc(n * dim * sizeof(double));
        for (int i = 1; i < n; ++i) {
            x[i] = x[0] + i * dim;
        }
        memcpy(x[0], shared_x[0], n * dim * sizeof(double));

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
        int counts[num_procs];
        int disps[num_procs];
        int counts_sol[num_procs];
        int disps_sol[num_procs];
        int current = 0;
        int best = 0 ;

        for (int i = 0; i < num_procs; i++){   
            int div = n / num_procs;
            offset = n % num_procs;
            disps_sol[i] = i * div + (i < offset ? i : offset);
            int proc_end = disps[i] + div + (i < offset);
            counts_sol[i] = proc_end - disps[i];
            counts[i] = counts_sol[i] * dim;
            disps[i] = disps_sol[i] * dim;
        }

        aux = counts_sol[rank] / n_threads;
        offset = counts_sol[rank] % n_threads;
        start = disps_sol[rank] + id * aux + (id < offset ? id : offset);
        end = start + aux + (id < offset);

        while (current < stopping_criteria){
            if (id == 0){
                generate_alphas(&rng, alpha);
                generate_lambdas(&rng, lambda);
            }
            pthread_barrier_wait(&barrier);
            update[id] = false;
            best = current_best_solution(n, dim, solutions);
            for(int i = start; i < end; i++){
                //inicializacion de randoms
                int j,k,m;
                double beta1, beta2, deltaj, deltak;
                get_JKM_betas_deltas(&rng, n, i, &j, &k, &m, &beta1, &beta2, &deltaj, &deltak);
                //Diversity phase
                diversity_phase(&rng, x[i], shared_x[best], shared_x[i], shared_x[j], shared_x[k], shared_x[m], P1, P2, P3, T2, T3, Vj, Vk, aux1, aux2, alpha, beta1, beta2, deltaj, deltak, lambda, dim);
                //Overlap phase
                overlap_phase(&rng, x[i], shared_x[best], shared_x[i], shared_x[j], shared_x[k], shared_x[m], P1, aux1, aux2, B1, B2, lambda, dim);
                //Make sure there are no elements that exceed the function bounds
                check_bounds(x[i], shared_x[i], lower_bound, upper_bound, dim);
            }
            pthread_barrier_wait(&barrier);
            for (int i = start; i < end; i++){
                double result_aux;
                evaluate_fitness(dim, x[i], &result_aux, tasks.objective_function);
                if (result_aux < solutions[i]){
                    update[id] = true;
                    solutions[i] = result_aux;
                    memcpy(shared_x[i], x[i], dim * sizeof(double));
                }
            }
            pthread_barrier_wait(&barrier);
            if (id == 0){
                MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, shared_x[0], counts, disps, MPI_DOUBLE, MPI_COMM_WORLD);
                MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &(solutions[0]), counts_sol, disps_sol, MPI_DOUBLE, MPI_COMM_WORLD);
            }
            current += n;
        }
        free(x[0]);
        free(x);
    }

    double* NOA(Mersenne_Twister* rng, int rank, int num_procs, int n_threads, int n, int stopping_criteria, int dim, double* lower_bound, double* upper_bound, double (*objective_function)(const double*, const int)){
        //initialize pointers
        double** x = (double**)malloc(n * sizeof(double*));
        x[0] = (double*)malloc(n * dim * sizeof(double));
        double* solutions = malloc(n * sizeof(double));
        bool* update = malloc(n_threads * sizeof(bool*));
        double alpha[4];
        int lambda[8];
        int best = 0;
        int current = 0;
    
        //initial population
        pthread_t threads[n_threads];
        struct params* tasks = malloc(n * sizeof(struct params));
        if (tasks == NULL)
            printf("Array of structs not reserved");
        pthread_barrier_init(&barrier, 0, n_threads);

        for(int i = 0; i < n_threads; i++){
            tasks[i].rng = *rng;
            tasks[i].id = i;
            tasks[i].n_threads = n_threads;
            tasks[i].n = n;
            tasks[i].dim = dim;
            tasks[i].rank = rank;
            tasks[i].num_procs = num_procs;
            tasks[i].stopping_criteria = stopping_criteria;
            tasks[i].lower_bound = lower_bound;
            tasks[i].upper_bound = upper_bound;
            tasks[i].ini = x[0];
            tasks[i].x = x;
            tasks[i].solutions = solutions;
            tasks[i].update = update;
            tasks[i].objective_function = objective_function;
            tasks[i].alpha = alpha;
            tasks[i].lambda = lambda;
            
            pthread_create(&threads[i], 0, parallel_task, &(tasks[i]));
        }
        for(int i = 0; i < n_threads; i++){
            pthread_join(threads[i], 0);
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
#else
    void* parallel_task(void* param){
        task tasks = *((task*) param);
        Mersenne_Twister rng = tasks.rng;
        int id = tasks.id;
        int n_threads = tasks.n_threads;
        int n = tasks.n;
        int dim = tasks.dim;
        int stopping_criteria = tasks.stopping_criteria;
        double** shared_x = tasks.x;
        double* lower_bound = tasks.lower_bound;
        double* upper_bound = tasks.upper_bound;
        double* solutions = tasks.solutions;
        double* alpha = tasks.alpha;
        int* lambda = tasks.lambda;
        pthread_mutex_lock(&mutex);
        unsigned int seed = getRandomNumber(&rng);
        pthread_mutex_unlock(&mutex);
        Create_MersenneTwister(&rng, seed);
        
        int aux = n / n_threads;
        int offset = n % n_threads;
        int start = id * aux + (id < offset ? id : offset);
        int end = start + aux + (id < offset);
        //initialize population
        for(int i = start; i < end; i++){
            shared_x[i] = tasks.ini + (i * dim);
            
            for (int j = 0; j < dim; j++){
                shared_x[i][j] = lower_bound[j] + (getRandomNumber(&rng) / (double)UINT_MAX) * (upper_bound[j] - lower_bound[j]);
            }
            //fitness value
            evaluate_fitness(dim, shared_x[i], &(solutions[i]), tasks.objective_function);
        }
        double** x = (double**)malloc(n * sizeof(double*));
        x[0] = (double*)malloc(n * dim * sizeof(double));
        for (int i = 1; i < n; ++i) {
            x[i] = x[0] + i * dim;
        }
        memcpy(x[0], shared_x[0], n * dim * sizeof(double));

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
        int current = 0;
        int best = 0 ;

        while (current < stopping_criteria){
            if (id == 0){
                generate_alphas(&rng, alpha);
                generate_lambdas(&rng, lambda);
            }
            pthread_barrier_wait(&barrier);
            best = current_best_solution(n, dim, solutions);
            for(int i = start; i < end; i++){
                //inicializacion de randoms
                int j,k,m;
                double beta1, beta2, deltaj, deltak;
                get_JKM_betas_deltas(&rng, n, i, &j, &k, &m, &beta1, &beta2, &deltaj, &deltak);
                //Diversity phase
                diversity_phase(&rng, x[i], shared_x[best], shared_x[i], shared_x[j], shared_x[k], shared_x[m], P1, P2, P3, T2, T3, Vj, Vk, aux1, aux2, alpha, beta1, beta2, deltaj, deltak, lambda, dim);
                //Overlap phase
                overlap_phase(&rng, x[i], shared_x[best], shared_x[i], shared_x[j], shared_x[k], shared_x[m], P1, aux1, aux2, B1, B2, lambda, dim);
                //Make sure there are no elements that exceed the function bounds
                check_bounds(x[i], shared_x[i], lower_bound, upper_bound, dim);
            }
            pthread_barrier_wait(&barrier);
            for (int i = start; i < end; i++){
                double result_aux;
                evaluate_fitness(dim, x[i], &result_aux, tasks.objective_function);
                if (result_aux < solutions[i]){
                    solutions[i] = result_aux;
                    memcpy(shared_x[i], x[i], dim * sizeof(double));
                }
            }
            current+=n;
        }
        free(x[0]);
        free(x);
    }

    double* NOA(Mersenne_Twister* rng, int n_threads, int n, int stopping_criteria, int dim, double* lower_bound, double* upper_bound, double (*objective_function)(const double*, const int)){
        //initialize pointers
        double** x = (double**)malloc(n * sizeof(double*));
        x[0] = (double*)malloc(n * dim * sizeof(double));
        double* solutions = malloc(n * sizeof(double));
        double alpha[4];
        int lambda[8];
        int best = 0;
        int current = 0;

        //initial population
        pthread_t threads[n_threads];
        struct params* tasks = malloc(n * sizeof(struct params));
        if (tasks == NULL)
            printf("Array of structs not reserved");
        pthread_barrier_init(&barrier, 0, n_threads);

        for(int i = 0; i < n_threads; i++){
            tasks[i].rng = *rng;
            tasks[i].id = i;
            tasks[i].n_threads = n_threads;
            tasks[i].n = n;
            tasks[i].dim = dim;
            tasks[i].stopping_criteria = stopping_criteria;
            tasks[i].lower_bound = lower_bound;
            tasks[i].upper_bound = upper_bound;
            tasks[i].ini = x[0];
            tasks[i].x = x;
            tasks[i].solutions = solutions;
            tasks[i].objective_function = objective_function;
            tasks[i].alpha = alpha;
            tasks[i].lambda = lambda;
            
            pthread_create(&threads[i], 0, parallel_task, &(tasks[i]));
        }
        for(int i = 0; i < n_threads; i++){
            pthread_join(threads[i], 0);
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

#endif
