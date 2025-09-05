#include "benchmark_functions.h"

double sphere(const double* x, const int dimension){
    double result = 0.0;
    for (int i = 0; i < dimension; i++)
        result = result + x[i] * x[i];
    
    return result;
}

double quartic(const double* x, const int dimension){
    double result = 0.0;
    for (int i = 0; i < dimension; i++)
        result = result + i * x[i] * x[i] * x[i] * x[i];
    
    return result;
}

double powell_sum(const double* x, const int dimension){
    double result = 0.0;
    for (int i = 0; i < dimension; i++){
        result = result + pow(fabs(x[i]),i+1);
    }
    
    return result;
}

double sum_squares(const double* x, const int dimension){
    double result = 0.0;
    for (int i = 0; i < dimension; i++)
        result = result + i * x[i] * x[i];
    
    return result;
}

double schwefel_2_20(const double* x, const int dimension){
    double result = 0.0;
    for (int i = 0; i < dimension; i++){
        result = result + fabs(x[i]);
    }
    
    return result;
}

double stepint(const double* x, const int dimension){
    double result = 0.0;
    for (int i = 0; i < dimension; i++)
        result = result + floor(x[i]);
    result = result + 25;
    
    return result;
}

double ridge(const double* x, const int dimension){
    double result = 0.0;
    for (int i = 1; i < dimension; i++)
        result = result + x[i] * x[i];
    result = x[0] + sqrt(result);
    
    return result;
}

double neumaier_N3(const double* x, const int dimension){
    double result = 0.0, sum1 = 0.0, sum2 = 0.0;
    sum1 = sum1 + (x[0] - 1) * (x[0] - 1);
    for (int i = 1; i < dimension; i++){
        sum1 = sum1 + (x[i] - 1) * (x[i] - 1);
        sum2 = sum2 + x[i] * x[i-1];
    }
    result = sum1 - sum2;
    
    return result;
}

double ackley_N2(const double* x, const int dimension){
    double result = -200 * exp(-0.02 * sqrt(x[0] * x[0] + x[1] * x[1]));
    
    return result;
}

double shekel_10(const double* x, const int dimension){
    double result = 0.0;
    for (int i = 0; i < 10; i++){
        double sum = 0.0;

        for (int j = 0; j < 4; j++){
            sum = sum + (x[j] - a[i][j]) * (x[j] - a[i][j]);
        }
        result = result + 1 / (sum + c[i]);
    }
    
    return -result;
}