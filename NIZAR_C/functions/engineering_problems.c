#include "engineering_problems.h"

void invert_matrix(double matrix[N][N], double inverse[N][N]) {
    double aux[N][N];

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            aux[i][j] = matrix[i][j];

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (i == j)    
                inverse[i][j] = 1.0;
            else
                inverse[i][j] = 0.0;
        }
    }

    for (int i = 0; i < N; i++){
        double pivot = aux[i][i];
        for (int j = 0; j < N; j++){
            aux[i][j] /= pivot;
            inverse[i][j] /= pivot;
        }

        for (int j = 0; j < N; j++){
            if (i != j){
                double factor = aux[j][i];
                for (int k = 0; k < N; k++){
                    aux[j][k] -= factor * aux[i][k];
                    inverse[j][k] -= factor * inverse[i][k];
                }
            }
        }
    }
}

double matrix_operations(){
    double matrix[N][N], inverse[N][N], aux[N][N];
    double determinant = 1;

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++){
            if (i == j){
                matrix[i][i] = 1.0;
                aux[i][i] = 1.0;
            }
            else{
                matrix[i][j] = 0.0;
                aux[i][j] = 0.0;
            }
        }

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            for (int k = 0; k < N; k++){
                matrix[i][k] += aux[j][k] * j;
            }
        }
    }

    invert_matrix(matrix, inverse);
    
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            aux[i][j] = 0.0;
    
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                aux[i][j] += matrix[i][k] * inverse[k][j];
    
    for (int i = 0; i < N; i++){
        determinant *= aux[i][i];
    }
    
    return determinant;
}

double pressure_vessel_design(const double* x, const int dimension){
    double result = 0.6224 * x[0] * x[2] * x[3] + 1.7781 * x[1] * x[2]*x[2] + 3.1661 * x[0]*x[0] * x[3] + 19.84 * x[0]*x[0] * x[2];
    double g1 = -x[0] + 0.0193 * x[2];
    double g2 = -x[1] + 0.00954 * x[2];
    double g3 = -M_PI * x[2]*x[2] * x[3] - (4.0 / 3.0) * M_PI * x[2]*x[2]*x[2] + 1296000;
    double g4 = x[3] - 240;
    
    //if (g1 > 0 || g2 > 0 || g3 > 0 || g4 > 0)
        return result + (10000 + g1) * (g1>0) + (10000 + g2) * (g2>0) + (10000 + g3) * (g3>0) + (10000 + g4) * (g4>0);
    //return result;
}

double tension_compression_spring_design(const double* x, const int dimension){
    double result = (x[2] + 2) * x[0]*x[0] * x[1];
    double g1 = 1 - (x[1]*x[1]*x[1] * x[2]) / (71785 * x[0]*x[0]*x[0]*x[0]);
    double g2 = (4 * x[1]*x[1] - x[0] * x[1]) / (12566 * (x[1] * x[0]*x[0]*x[0] - x[0]*x[0]*x[0]*x[0])) + 1 / (5108 * x[0]*x[0]) - 1;
    double g3 = 1 - (140.45 * x[0]) / (x[1]*x[1] * x[2]);
    double g4 = (x[1] + x[0]) / 1.5 - 1;
    
    return result + g1 * (g1>0) + g2 * (g2>0) + g3 * (g3>0) + g4 * (g4>0);
}

double artificial_function(const double* x, const int dimension){
    double result = (x[2] + 2) * x[0]*x[0] * x[1];
    double g1 = 1 - (x[1]*x[1]*x[1] * x[2]) / (71785 * x[0]*x[0]*x[0]*x[0]);
    double g2 = (4 * x[1]*x[1] - x[0] * x[1]) / (12566 * (x[1] * x[0]*x[0]*x[0] - x[0]*x[0]*x[0]*x[0])) + 1 / (5108 * x[0]*x[0]) - 1;
    double g3 = 1 - (140.45 * x[0]) / (x[1]*x[1] * x[2]);
    double g4 = (x[1] + x[0]) / 1.5 - 1;
    double determinant = matrix_operations();

    return determinant * result + g1 * (g1>0) + g2 * (g2>0) + g3 * (g3>0) + g4 * (g4>0);
}