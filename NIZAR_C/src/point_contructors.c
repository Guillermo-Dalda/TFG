#include "point_contructors.h"

void reconstruct_S1(double* X_p, double* P1, double* x_i, double* x_j, double* x_k, double* aux1, double* aux2, double* B1, double* B2, int* lambda, int dim){
    if (lambda[1] == 1){
        dif(x_i,  x_j, aux1, dim);
        dif(x_i, x_k, aux2, dim);
        vector_elements_product(B1, aux1, aux1, dim);
        vector_elements_product(B2, aux2, aux2, dim);
        sum(aux1, aux2, X_p, dim);
        sum(P1, X_p, X_p, dim);
    }
    else{
        dif(P1,  x_j, aux1, dim);
        dif(P1, x_k, aux2, dim);
        vector_elements_product(B1, aux1, aux1, dim);
        vector_elements_product(B2, aux2, aux2, dim);
        dif(aux1, aux2, X_p, dim);
        sum(x_i, X_p, X_p, dim);
    }
}

void reconstruct_P1(double* P1, double* x_best, double* x_m, double* aux1, double* aux2, double* B1, int* lambda, int dim){
    if (lambda[2] == 1)
        memcpy(P1, x_m, dim * sizeof(double));
    //Else construct T1 as follows...
    else if (lambda[5] == 1){
        sum(x_best, x_m, P1, dim);
        mul(0.5, P1, P1, dim);
    }
    else{
        vector_elements_product(B1, x_m, aux1, dim);
        for (int i = 0; i < dim; i++){
            aux2[i] = 1;
        }
        dif(aux2, B1, aux2, dim);
        vector_elements_product(aux2, x_best, aux2, dim);
        sum(aux1, aux2, P1, dim);
    }
}

void construct_S1(double* X_p, double* P1, double* x_i, double* x_j, double* x_k, double* aux1, double* aux2, double beta1, double beta2, int* lambda, int dim){
    if (lambda[1] == 1){
        dif(x_i,  x_j, aux1, dim);
        dif(x_i, x_k, aux2, dim);
        mul(beta1, aux1, aux1, dim);
        mul(beta2, aux2, aux2, dim);
        sum(aux1, aux2, X_p, dim);
        sum(P1, X_p, X_p, dim);
    }
    else{
        dif(P1,  x_j, aux1, dim);
        dif(P1, x_k, aux2, dim);
        mul(beta1, aux1, aux1, dim);
        mul(beta2, aux2, aux2, dim);
        dif(aux1, aux2, X_p, dim);
        sum(x_i, X_p, X_p, dim);
    }
}

void construct_P1(double* P1, double* x_best, double* x_m, double* aux1, double* aux2, double beta1, int* lambda, int dim){
    if (lambda[2] == 1)
        memcpy(P1, x_m, dim * sizeof(double));
    //Else construct T1 as follows...
    else if (lambda[5] == 1){
        sum(x_best, x_m, P1, dim);
        mul(0.5, P1, P1, dim);
    }
    else{
        mul(beta1, x_m, aux1, dim);
        mul((1 - beta1), x_best, aux2, dim);
        sum(aux1, aux2, P1, dim);
    }
}

void construct_S2(double* X_p, double* P2, double* P3, double* Vj, double* Vk, double* aux1, double* aux2, double deltaj, double deltak, int* lambda, int dim){
    if (lambda[1] == 1){
        dif(P3, Vj, aux1, dim);
        dif(P3, Vk, aux2, dim);
        mul(deltaj, aux1, aux1, dim);
        mul(deltak, aux2, aux2, dim);
        sum(aux1, aux2, X_p, dim);
        sum(P2, X_p, X_p, dim);
    }
    else{
        dif(P2, Vj, aux1, dim);
        dif(P2, Vk, aux2, dim);
        mul(deltaj, aux1, aux1, dim);
        mul(deltak, aux2, aux2, dim);
        dif(aux1, aux2, X_p, dim);
        sum(P3, X_p, X_p, dim);
    }
}

void construct_P2(double* P2, double* T2, double* x_m, double* aux1, double beta1, int* lambda, int dim){
    if (lambda[3] == 1){
        for (int i = 0; i < dim; i++){
            aux1[i] = beta1;
        }
        sum(x_m, aux1, P2, dim);
    }
    else {
        memcpy(P2, T2, dim * sizeof(double));
    }
}

void construct_P3(Mersenne_Twister* rng, double* P3, double* P1, double* T2, double* T3, double* alpha, double beta2, int* lambda, int dim){
    if (lambda[4] == 1){
        distribute(rng, T3, P1, P3, dim);
        transfer_map(P3, T2, alpha[2], P3, dim);
        dilation_map(P3, beta2, alpha[1], P3, dim);
    }
    else{
        distribute(rng, P1, T3, P3, dim);
        transfer_map(P3, T2, alpha[2], P3, dim);
        dilation_map(P3, beta2, alpha[1], P3, dim);
    }
}

void construct_T2(Mersenne_Twister* rng, double* T2, double* T3, double* P1, double* alpha, double beta2, int* lambda, int dim){
    if (lambda[6] == 1){
        replace(rng, T3, P1, T2, dim);
        translation_map(T2, beta2, alpha[3], T2, dim);
    }
    else{
        scramble(rng, T3, P1, T2, dim);
        translation_map(T2, beta2, alpha[3], T2, dim);
    }
}

void construct_T3(double* T3, double* x_best, double* x_i, int* lambda, int dim){
    if (lambda[7] == 1)
        memcpy(T3, x_best, dim * sizeof(double));
    else
        memcpy(T3, x_i, dim * sizeof(double));
}