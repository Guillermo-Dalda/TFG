#include "vector_operations.h"

bool vector_comparison(double* u, double* v, int dim){
    for (int i = 0; i < dim; i++){
        if (fabs(u[i] - v[i]) > EPSILON)
            return false;
    }
    return true;
}

void vector_elements_product(double* u, double* v, double* result, int dim){
    for(int i = 0; i < dim; i++){
        result[i] = u[i] * v[i];
    }
}

void mul(double a, double* v, double* result, int dim){
    for(int i = 0; i < dim; i++){
        result[i] = a * v[i];
    }
}

void sum(double* u, double* v, double* result, int dim){
    for(int i = 0; i < dim; i++){
        result[i] = u[i] + v[i];
    }
}

void dif(double* u, double* v, double* result, int dim){
    for(int i = 0; i < dim; i++){
        result[i] = u[i] - v[i];
    }
}