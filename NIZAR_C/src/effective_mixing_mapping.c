#include "effective_mixing_mapping.h"

void replace(Mersenne_Twister* rng, double* Xin, double* Xtg, double* result, int dim){
    for (int i = 0; i < dim; i++){
        if (getRandomNumber(rng) % 2 == 1)
            result[i] = Xtg[i];
        else
            result[i] = Xin[i];
    }
}

void scramble(Mersenne_Twister* rng, double* Xin, double* Xtg, double* result, int dim){
    int index[dim];
    for (int i = 0; i < dim; i++)
        index[i] = i;
    for (int i = dim-1; i >= 0; i--){
        int j = getRandomNumber(rng) % (i+1);
        int aux = index[j];
        index[j] = index[i];
        index[i] = aux;
    }
    for (int i = 0; i < dim; i++){
        if (getRandomNumber(rng) % 2 == 1)
            result[i] = Xtg[index[i]];
        else
            result[i] = Xin[i];
    }
}

void distribute(Mersenne_Twister* rng, double* Xin, double* Xtg, double* result, int dim){
    int pos_tg = getRandomNumber(rng) % dim;
    for (int i = 0; i < dim; i++){
        if (getRandomNumber(rng) % 2 == 1)
            result[i] = Xtg[pos_tg];
        else
            result[i] = Xin[i];
    }
}