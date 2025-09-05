#include <stdio.h>
#include <math.h>

//shekel_10 matrix a
static const double a[10][4] = {
    {4.0, 4.0, 4.0, 4.0},
    {1.0, 1.0, 1.0, 1.0},
    {8.0, 8.0, 8.0, 8.0},
    {6.0, 6.0, 6.0, 6.0},
    {3.0, 7.0, 3.0, 7.0},
    {2.0, 9.0, 2.0, 9.0},
    {5.0, 5.0, 3.0, 3.0},
    {8.0, 1.0, 8.0, 1.0},
    {6.0, 2.0, 6.0, 2.0},
    {7.0, 3.6, 7.0, 3.6},
};
//shekel_10 matrix c
static const double c[10] = {0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3, 0.7, 0.5, 0.5};

double sphere(const double* x, const int dimension);
double quartic(const double* x, const int dimension);
double powell_sum(const double* x, const int dimension);
double sum_squares(const double* x, const int dimension);
double schwefel_2_20(const double* x, const int dimension);
double stepint(const double* x, const int dimension);
double ridge(const double* x, const int dimension);
double neumaier_N3(const double* x, const int dimension);
double ackley_N2(const double* x, const int dimension);
double shekel_10(const double* x, const int dimension);