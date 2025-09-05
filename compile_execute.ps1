param($opt = 'seq', $p = 2, $t = 2, $n = 25, $eval = 35000)

if ($opt -eq 'compile'){
    #sequential compilation
    $source = Get-ChildItem src\* | Where-Object {$_.Name -notmatch "par_NIZAR.c" -and $_.Name -notmatch "par_main.c"}
    gcc -O2 -o bin\seq_optimizer functions\* $source -I include -I functions
    #pthread veersion compilation
    $source = Get-ChildItem src\* | Where-Object {$_.Name -notmatch "seq_NIZAR.c" -and $_.Name -notmatch "seq_main.c"}
    gcc -O2 -o bin\pthr_optimizer functions\* $source -I include -I functions
    #MPI+pthread veersion compilation
    $include="C:\Program Files (x86)\Microsoft SDKs\MPI\Include"
    $lib="C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64"
    $source = Get-ChildItem src\* | Where-Object {$_.Name -notmatch "seq_NIZAR.c" -and $_.Name -notmatch "seq_main.c"}
    mpicc -O2 -o bin\mpi_optimizer functions\* $source -I include -I functions -DMPI
}
if ($opt -eq 'seq'){
    .\bin\seq_optimizer.exe 0 $n $eval 50 -100 100
    .\bin\seq_optimizer.exe 1 $n $eval 50 -1.28 1.28
    .\bin\seq_optimizer.exe 2 $n $eval 50 -1 1
    .\bin\seq_optimizer.exe 3 $n $eval 50 -10 10
    .\bin\seq_optimizer.exe 4 $n $eval 50 -100 100
    .\bin\seq_optimizer.exe 5 $n $eval 50 -5.12 5.12
    .\bin\seq_optimizer.exe 6 $n $eval 50 -5 5
    .\bin\seq_optimizer.exe 7 $n $eval 15 -100 100
    .\bin\seq_optimizer.exe 8 $n $eval 2 -32 32
    .\bin\seq_optimizer.exe 9 $n $eval 4 0 10
    .\bin\seq_optimizer.exe 10 $n $eval 4 0 100 0 100 10 200 10 200
    .\bin\seq_optimizer.exe 11 $n $eval 3 0.05 2 0.25 1.3 2 15
    .\bin\seq_optimizer.exe 12 $n $eval 3 0.05 2 0.25 1.3 2 15
}
if ($opt -eq 'pthr'){
    .\bin\pthr_optimizer.exe 0 $t $n $eval 50 -100 100
    .\bin\pthr_optimizer.exe 1 $t $n $eval 50 -1.28 1.28
    .\bin\pthr_optimizer.exe 2 $t $n $eval 50 -1 1
    .\bin\pthr_optimizer.exe 3 $t $n $eval 50 -10 10
    .\bin\pthr_optimizer.exe 4 $t $n $eval 50 -100 100
    .\bin\pthr_optimizer.exe 5 $t $n $eval 50 -5.12 5.12
    .\bin\pthr_optimizer.exe 6 $t $n $eval 50 -5 5
    .\bin\pthr_optimizer.exe 7 $t $n $eval 15 -100 100
    .\bin\pthr_optimizer.exe 8 $t $n $eval 2 -32 32
    .\bin\pthr_optimizer.exe 9 $t $n $eval 4 0 10
    .\bin\pthr_optimizer.exe 10 $t $n $eval 4 0 100 0 100 10 200 10 200
    .\bin\pthr_optimizer.exe 11 $t $n $eval 3 0.05 2 0.25 1.3 2 15
    .\bin\pthr_optimizer.exe 12 $t $n $eval 3 0.05 2 0.25 1.3 2 15
}
if ($opt -eq 'mpi'){    
    mpiexec -n $p .\bin\mpi_optimizer.exe 0 $t $n $eval 50 -100 100
    mpiexec -n $p .\bin\mpi_optimizer.exe 1 $t $n $eval 50 -1.28 1.28
    mpiexec -n $p .\bin\mpi_optimizer.exe 2 $t $n $eval 50 -1 1
    mpiexec -n $p .\bin\mpi_optimizer.exe 3 $t $n $eval 50 -10 10
    mpiexec -n $p .\bin\mpi_optimizer.exe 4 $t $n $eval 50 -100 100
    mpiexec -n $p .\bin\mpi_optimizer.exe 5 $t $n $eval 50 -5.12 5.12
    mpiexec -n $p .\bin\mpi_optimizer.exe 6 $t $n $eval 50 -5 5
    mpiexec -n $p .\bin\mpi_optimizer.exe 7 $t $n $eval 15 -100 100
    mpiexec -n $p .\bin\mpi_optimizer.exe 8 $t $n $eval 2 -32 32
    mpiexec -n $p .\bin\mpi_optimizer.exe 9 $t $n $eval 4 0 10
    mpiexec -n $p .\bin\mpi_optimizer.exe 10 $t $n $eval 4 0 100 0 100 10 200 10 200
    mpiexec -n $p .\bin\mpi_optimizer.exe 11 $t $n $eval 3 0.05 2 0.25 1.3 2 15
    mpiexec -n $p .\bin\mpi_optimizer.exe 12 $t $n $eval 3 0.05 2 0.25 1.3 2 15
}