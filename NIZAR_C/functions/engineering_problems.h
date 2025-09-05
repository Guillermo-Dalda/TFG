#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>

static const int N = 50;

static void invert_matrix(double matrix[N][N], double inverse[N][N]);
static double matrix_operations();

double pressure_vessel_design(const double* x, const int dimension);
double tension_compression_spring_design(const double* x, const int dimension);
double artificial_function(const double* x, const int dimension);