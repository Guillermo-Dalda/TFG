/* *********************************** */
//  Author: Guillermo Dalda del Olmo
//  
//
/* *********************************** */

#include "engineering_problems.h"
#include "benchmark_functions.h"
#include "seq_NIZAR.h"
#include <time.h>
#include <sys/time.h>

int main(int argc, char* argv[]){
    Mersenne_Twister rng;
    int seed = (int)time(0);//5489
    double (*objective_function[13])(const double*, const int);
    objective_function[0] = sphere;
    objective_function[1] = quartic;
    objective_function[2] = powell_sum;
    objective_function[3] = sum_squares;
    objective_function[4] = schwefel_2_20;
    objective_function[5] = stepint;
    objective_function[6] = ridge;
    objective_function[7] = neumaier_N3;
    objective_function[8] = ackley_N2;
    objective_function[9] = shekel_10;
    objective_function[10] = pressure_vessel_design;
    objective_function[11] = tension_compression_spring_design;
    objective_function[12] = artificial_function;

    if (argc < 7){
        printf("arguments: function, population size, stopping criteria, dimension, lower bound 1, upper bound 1, lower bound 2, upper bound 2..., rand seed(optional)\n");
    }
    else{
        int function = atoi(argv[1]);
        int population_size = atoi(argv[2]);
        if (population_size < 4){
            printf("El tamaño de población ha de ser minimo 4");
        }
        else{
            int stopping_criteria = atoi(argv[3]);
            int dimension = atoi(argv[4]);
            
            switch (function){
                case 0:
                    printf("\nSphere:\n");
                    break;
                case 1:
                    printf("\nQuartic:\n");
                    break;
                case 2:
                    printf("\nPowell_sum:\n");
                    break;
                case 3:
                    printf("\nSum_squares:\n");
                    break;
                case 4:
                    printf("\nSchwefel's 2.20:\n");
                    break;
                case 5:
                    printf("\nStepint:\n");
                    break;
                case 6:
                    printf("\nRidge:\n");
                    break;
                case 7:
                    printf("\nNeumaier's N. 3:\n");
                    break;
                case 8:
                    printf("\nAckley N. 2:\n");
                    break;
                case 9:
                    printf("\nShekel 10:\n");
                    break;
                case 10:
                    printf("\nPVD:\n");
                    break;
                case 11:
                    printf("\nTCSD:\n");
                    break;
                case 12:
                    printf("\nArtificial function:\n");
                    break;
                default:
                    printf("Functions: 0-12\n");
                    return 1;
            }

            double lower_bound[dimension];
            double upper_bound[dimension];
            //if not given lower bound and upper bound for every dimension the first bounds will apply for all dimensions
            if (argc < 2 * dimension + 5){
                for (int i = 0; i < dimension; i++){
                    lower_bound[i] = atof(argv[5]);
                    upper_bound[i] = atof(argv[6]);
                }
            }
            //if bounds are given for all dimensions
            else if (argc >= 2 * dimension + 5)
                for (int i = 0; i < dimension; i++){
                    lower_bound[i] = atof(argv[5+i*2]);
                    upper_bound[i] = atof(argv[6+i*2]);
                }
            if (argc > 2 * dimension + 6)
                seed = atoi(argv[argc-1]);
            Create_MersenneTwister(&rng, seed);
                
            struct timeval start, end;
            gettimeofday(&start, NULL);
            
            double* solution = NOA(&rng, population_size, stopping_criteria, dimension, lower_bound, upper_bound, objective_function[function]);

            gettimeofday(&end, NULL);
			double time_elapsed = (end.tv_sec-start.tv_sec) + (end.tv_usec-start.tv_usec) * 0.000001;
			            
            printf("Solution:\n");
            for(int i = 0; i < dimension; i++)
                printf("%lf ",solution[i]);
            printf("\nOptimum: %lf\n", solution[dimension]);
            printf("Time: %lf\n", time_elapsed);
            free(solution);
        }
    }
    return 0;
}