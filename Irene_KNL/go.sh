#!/bin/bash
#MSUB -q knl
#MSUB -n 64
#MSUB -N 1
#MSUB -c 1
#MSUB -r test
#MSUB -T 300                  # Elapsed time limit in seconds
#MSUB -o out_%I.o              # Standard output. %I is the job id
#MSUB -e err_%I.e              # Error output. %I is the job id

module purge
module load mpi
module load fftw3

ccc_mprun ./sc_compiled/flow36

######### submit script ccc_msub
######### check queue with ccc_mpp -u user
##### N : number of nodes (64 cores x nodes)
##### n : number of MPI processes per node
##### c : Hyperthreaing
