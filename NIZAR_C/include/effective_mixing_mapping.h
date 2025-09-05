#include "vector_operations.h"
#include "mersenne.h"
#include <stdlib.h>

void replace(Mersenne_Twister* rng, double* Xin, double* Xtg, double* result, int len);
void scramble(Mersenne_Twister* rng, double* Xin, double* Xtg, double* result, int len);
void distribute(Mersenne_Twister* rng, double* Xin, double* Xtg, double* result, int len);