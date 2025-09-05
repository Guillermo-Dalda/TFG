#include "point_contructors.h"
#include <stdbool.h>
#include <limits.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

static void diversity_phase(Mersenne_Twister* rng,double* X_p, double* x_best, double* x_i, double* x_j, double* x_k, double* x_m, double* P1, double* P2, double* P3, double* T2, double* T3, double* Vj, double* Vk, double* aux1, double* aux2, double* alpha, double beta1, double beta2, double deltaj, double deltak, int* lambda, int dim);
static void overlap_phase(Mersenne_Twister* rng,double* X_p, double* x_best, double* x_i, double* x_j, double* x_k, double* x_m, double* P1, double* aux1, double* aux2, double* B1,  double* B2, int* lambda, int dim);
static void check_bounds(double* X_p, double* x_i, double* lower_bound, double* upper_bound, int dim);

static void generate_alphas(Mersenne_Twister* rng,double* alpha);
static void generate_lambdas(Mersenne_Twister* rng,int* lambda);
static void get_JKM_betas_deltas(Mersenne_Twister* rng, int n, int i, int* j, int* k, int* m, double* beta1, double* beta2, double* deltaj, double* deltak);
static int current_best_solution(int n, int dimension, double* solutions);
static void evaluate_fitness(int dimension, double* x, double* solution, double (*objective_function)(const double*, const int));

double* NOA(Mersenne_Twister* rng,int n, int stopping_criteria, int  dim, double* lower_bound, double* upper_bound, double (*objective_function)(const double*, const int));