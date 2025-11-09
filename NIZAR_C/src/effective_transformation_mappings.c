#include "effective_transformation_mappings.h"

void translation_map(double* X, double r, double alpha, double* result, int dim){
    if (alpha <= 0.5)
        if (X != result)
            memcpy(result, X, dim * sizeof(double));
    else{
        double r2 = r*r;
        for (int i = 0; i < dim; i++){
            result[i] = round(X[i] + r2);
        }
    }
}

void dilation_map(double* X, double r, double alpha, double* result, int dim){
    if (alpha <= 0.5)
        if (X != result)
            memcpy(result, X, dim * sizeof(double));
    else{
        mul(r*r, X, result, dim);
        result = round(result);
    }
}

void transfer_map(double* Xi, double* Xj, double alpha, double* result, int dim){
    if (alpha <= 0.5)
        if (Xi != result)
            memcpy(result, Xi, dim * sizeof(double));
    else
        if (Xj != result)
            memcpy(result, Xj, dim * sizeof(double));

}
