#include "engineering_problems.h"
#include "benchmark_functions.h"
#include "par_NIZAR.h"
#include <time.h>
#include <sys/time.h>

#ifdef MPI
    #define MASTER 0
    #define NORMAL_TAG 1
    #define DESTROY_TAG 9

    void end_process(int numProcs){
        for(int i = 1; i<numProcs; i++){
            MPI_Send(0, 0, MPI_INT, i, DESTROY_TAG, MPI_COMM_WORLD);
        }
    }
#endif

int main(int argc, char* argv[]){
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

    #ifdef MPI
        MPI_Init(&argc, &argv);
        int rank, num_procs;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
        
        if (rank == MASTER){
            if (argc < 8){
                printf("arguments: function, number of threads, population size, stopping criteria, dimension, lower bound 1, upper bound 1, lower bound 2, upper bound 2..., rand seed(optional)\n");
                end_process(num_procs);
            }
            else{
                int function = atoi(argv[1]);
                int n_threads = atoi(argv[2]);
                int population_size = atoi(argv[3]);
                if (population_size < 4){
                    printf("El tamaño de población ha de ser minimo 4");
                    end_process(num_procs);
                }
                else{
                    int stopping_criteria = atoi(argv[4]);
                    int dimension = atoi(argv[5]);
                    
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
                            end_process(num_procs);
                            MPI_Finalize();
                            return 1;
                    }

                    double lower_bound[dimension];
                    double upper_bound[dimension];
                    //if not given lower bound and upper bound for every dimension the first bounds will apply for all dimensions
                    if (argc < 2 * dimension + 6){
                        for (int i = 0; i < dimension; i++){
                            lower_bound[i] = atof(argv[6]);
                            upper_bound[i] = atof(argv[7]);
                        }
                    }
                    //if bounds are given for all dimensions
                    else if (argc >= 2 * dimension + 6)
                        for (int i = 0; i < dimension; i++){
                            lower_bound[i] = atof(argv[6+i*2]);
                            upper_bound[i] = atof(argv[7+i*2]);
                        }
                    Mersenne_Twister rng;
                    if (argc > 2 * dimension + 6)
                        seed = atoi(argv[argc-1]);
                    Create_MersenneTwister(&rng, seed);
                    
                    int data[6];
                    data[0] = seed;
                    data[1] = function;
                    data[2] = n_threads;
                    data[3] = population_size;
                    data[4] = stopping_criteria;
                    data[5] = dimension;
                    struct timeval start, end;
                    gettimeofday(&start, NULL);
                    
                    for (int i = 1; i < num_procs; i++){
                        MPI_Send(&(data[0]), 6, MPI_INT, i, NORMAL_TAG, MPI_COMM_WORLD);
                        MPI_Send(&(lower_bound[0]), dimension, MPI_DOUBLE, i, NORMAL_TAG, MPI_COMM_WORLD);
                        MPI_Send(&(upper_bound[0]), dimension, MPI_DOUBLE, i, NORMAL_TAG, MPI_COMM_WORLD);
                    }
                    double* solution = NOA(&rng, rank, num_procs, n_threads, population_size, stopping_criteria, dimension, lower_bound, upper_bound, objective_function[function]);
                    
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
        }
        else{
            struct timeval start, end;
            gettimeofday(&start, NULL);
            MPI_Status status;
            MPI_Probe(MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            
            if (status.MPI_TAG == NORMAL_TAG){
                int buffer[6];
                MPI_Recv(buffer, 6, MPI_INT, MASTER, NORMAL_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                int seed = buffer[0];
                int function = buffer[1];
                int n_threads = buffer[2];
                int population_size = buffer[3];
                int stopping_criteria = buffer[4];
                int dimension = buffer[5];
                double lower_bound[dimension];
                double upper_bound[dimension];

                MPI_Recv(lower_bound, dimension, MPI_DOUBLE, MASTER, NORMAL_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(upper_bound, dimension, MPI_DOUBLE, MASTER, NORMAL_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                Mersenne_Twister rng;
                Create_MersenneTwister(&rng, seed+rank);
                
                double* solution = NOA(&rng, rank, num_procs, n_threads, population_size, stopping_criteria, dimension, lower_bound, upper_bound, objective_function[function]);
                free(solution);
            }
        }
        MPI_Finalize();
    #else
        if (argc < 8){
            printf("arguments: function, number of threads, population size, stopping criteria, dimension, lower bound 1, upper bound 1, lower bound 2, upper bound 2..., rand seed(optional)\n");
        }
        else{
            int function = atoi(argv[1]);
            int n_threads = atoi(argv[2]);
            int population_size = atoi(argv[3]);
            if (population_size < 4){
                printf("El tamaño de población ha de ser minimo 4");
            }
            else{
                int stopping_criteria = atoi(argv[4]);
                int dimension = atoi(argv[5]);
                
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
                if (argc < 2 * dimension + 6){
                    for (int i = 0; i < dimension; i++){
                        lower_bound[i] = atof(argv[6]);
                        upper_bound[i] = atof(argv[7]);
                    }
                }
                //if bounds are given for all dimensions
                else if (argc >= 2 * dimension + 6)
                    for (int i = 0; i < dimension; i++){
                        lower_bound[i] = atof(argv[6+i*2]);
                        upper_bound[i] = atof(argv[7+i*2]);
                    }
                Mersenne_Twister rng;
                if (argc > 2 * dimension + 6)
                    seed = atoi(argv[argc-1]);
                Create_MersenneTwister(&rng, seed);
                    
                struct timeval start, end;
                gettimeofday(&start, NULL);
                
                double* solution = NOA(&rng, n_threads, population_size, stopping_criteria, dimension, lower_bound, upper_bound, objective_function[function]);

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
    #endif
    return 0;
}