#!/bin/bash
opt="compile"
processes="2"
threads="2"
n="25"
eval="35000"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --opt) opt="$2"; shift 2 ;;
    --p) processes="$2"; shift 2 ;;
    --t) threads="$2"; shift 2 ;;
    --n) n="$2"; shift 2 ;;
    --eval) eval="$2"; shift 2 ;;
  esac
done

if [[ "$opt" == "compile" ]]; then
  #sequential compilation
  source=$(find src -name '*.c' ! -name 'par_NIZAR.c' ! -name 'par_main.c')
  gcc -O2 -o bin/seq_optimizer functions/*.c $source -I include -I functions -lm
  #pthread veersion compilation
  source=$(find src -name '*.c' ! -name 'seq_NIZAR.c' ! -name 'seq_main.c')
  gcc -O2 -o bin/pthr_optimizer functions/*.c $source -I include -I functions -lm
  #MPI+pthread veersion compilation
  source=$(find src -name '*.c' ! -name 'seq_NIZAR.c' ! -name 'seq_main.c')
  mpicc -O2 -o bin/mpi_optimizer functions/*.c $source -I include -I functions -lm  -DMPI

elif [[ "$opt" == "seq" ]]; then
  ./bin/seq_optimizer 0 $n $eval 50 -100 100
  ./bin/seq_optimizer 1 $n $eval 50 -1.28 1.28
  ./bin/seq_optimizer 2 $n $eval 50 -1 1
  ./bin/seq_optimizer 3 $n $eval 50 -10 10
  ./bin/seq_optimizer 4 $n $eval 50 -100 100
  ./bin/seq_optimizer 5 $n $eval 50 -5.12 5.12
  ./bin/seq_optimizer 6 $n $eval 50 -5 5
  ./bin/seq_optimizer 7 $n $eval 15 -100 100
  ./bin/seq_optimizer 8 $n $eval 2 -32 32
  ./bin/seq_optimizer 9 $n $eval 4 0 10
  ./bin/seq_optimizer 10 $n $eval 4 0 100 0 100 10 200 10 200
  ./bin/seq_optimizer 11 $n $eval 3 0.05 2 0.25 1.3 2 15
  ./bin/seq_optimizer 12 $n $eval 3 0.05 2 0.25 1.3 2 15

elif [[ "$opt" == "pthr" ]]; then
  ./bin/pthr_optimizer 0 $threads $n $eval 50 -100 100
  ./bin/pthr_optimizer 1 $threads $n $eval 50 -1.28 1.28
  ./bin/pthr_optimizer 2 $threads $n $eval 50 -1 1
  ./bin/pthr_optimizer 3 $threads $n $eval 50 -10 10
  ./bin/pthr_optimizer 4 $threads $n $eval 50 -100 100
  ./bin/pthr_optimizer 5 $threads $n $eval 50 -5.12 5.12
  ./bin/pthr_optimizer 6 $threads $n $eval 50 -5 5
  ./bin/pthr_optimizer 7 $threads $n $eval 15 -100 100
  ./bin/pthr_optimizer 8 $threads $n $eval 2 -32 32
  ./bin/pthr_optimizer 9 $threads $n $eval 4 0 10
  ./bin/pthr_optimizer 10 $threads $n $eval 4 0 100 0 100 10 200 10 200
  ./bin/pthr_optimizer 11 $threads $n $eval 3 0.05 2 0.25 1.3 2 15
  ./bin/pthr_optimizer 12 $threads $n $eval 3 0.05 2 0.25 1.3 2 15

elif [[ "$opt" == "mpi" ]]; then
  mpirun --bind-to none -n $processes ./bin/mpi_optimizer 0 $threads $n $eval 50 -100 100
  mpirun --bind-to none -n $processes ./bin/mpi_optimizer 1 $threads $n $eval 50 -1.28 1.28
  mpirun --bind-to none -n $processes ./bin/mpi_optimizer 2 $threads $n $eval 50 -1 1
  mpirun --bind-to none -n $processes ./bin/mpi_optimizer 3 $threads $n $eval 50 -10 10
  mpirun --bind-to none -n $processes ./bin/mpi_optimizer 4 $threads $n $eval 50 -100 100
  mpirun --bind-to none -n $processes ./bin/mpi_optimizer 5 $threads $n $eval 50 -5.12 5.12
  mpirun --bind-to none -n $processes ./bin/mpi_optimizer 6 $threads $n $eval 50 -5 5
  mpirun --bind-to none -n $processes ./bin/mpi_optimizer 7 $threads $n $eval 15 -100 100
  mpirun --bind-to none -n $processes ./bin/mpi_optimizer 8 $threads $n $eval 2 -32 32
  mpirun --bind-to none -n $processes ./bin/mpi_optimizer 9 $threads $n $eval 4 0 10
  mpirun --bind-to none -n $processes ./bin/mpi_optimizer 10 $threads $n $eval 4 0 100 0 100 10 200 10 200
  mpirun --bind-to none -n $processes ./bin/mpi_optimizer 11 $threads $n $eval 3 0.05 2 0.25 1.3 2 15
  mpirun --bind-to none -n $processes ./bin/mpi_optimizer 12 $threads $n $eval 3 0.05 2 0.25 1.3 2 15
fi
