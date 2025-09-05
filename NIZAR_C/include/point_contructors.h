#include "effective_transformation_mappings.h"
#include "effective_mixing_mapping.h"


void reconstruct_S1(double* X_p, double* P1, double* x_i, double* x_j, double* x_k, double* aux1, double* aux2, double* B1, double* B2, int* lambda, int dim);
void reconstruct_P1(double* P1, double* x_best, double* x_m, double* aux1, double* aux2, double* B1, int* lambda, int dim);
void construct_S1(double* X_p, double* P1, double* x_i, double* x_j, double* x_k, double* aux1, double* aux2, double beta1, double beta2, int* lambda, int dim);
void construct_P1(double* P1, double* x_best, double* x_m, double* aux1, double* aux2, double beta1, int* lambda, int dim);
void construct_S2(double* X_p, double* P2, double* P3, double* Vj, double* Vk, double* aux1, double* aux2, double deltaj, double deltak, int* lambda, int dim);
void construct_P2(double* P2, double* T2, double* x_m, double* aux1, double beta1, int* lambda, int dim);
void construct_P3(Mersenne_Twister* rng, double* P3, double* P1, double* T2, double* T3, double* alpha, double beta2, int* lambda, int dim);
void construct_T2(Mersenne_Twister* rng, double* T2, double* T3, double* P1, double* alpha, double beta2, int* lambda, int dim);
void construct_T3(double* T3, double* x_best, double* x_i, int* lambda, int dim);