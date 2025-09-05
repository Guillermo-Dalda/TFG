#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

#ifndef EPSILON
#define EPSILON 1e-12
#endif

bool vector_comparison(double* u, double* v, int len);
void vector_elements_product(double* u, double* v, double* result, int len);
void mul(double a, double* v, double* result, int len);
void sum(double* u, double* v, double* result, int len);
void dif(double* u, double* v, double* result, int len);