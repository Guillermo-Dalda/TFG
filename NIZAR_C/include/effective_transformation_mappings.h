#include "vector_operations.h"
#include <stdlib.h>
#include <string.h>

void translation_map(double* X, double r, double alpha, double* result, int len);
void dilation_map(double* X, double r, double alpha, double* result, int len);
void transfer_map(double* Xi, double* Xj, double alpha, double* result, int len);