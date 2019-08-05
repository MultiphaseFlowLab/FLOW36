#!/bin/bash

#SBATCH -N1 -n1                             # n tasks on N node (total number of tasks)
#SBATCH --time=5:00                         # time limits: hh:mm:ss
#SBATCH --error error.err                   # std-error file
#SBATCH --output output.out                 # std-output file
#SBATCH --job-name=name                     # job name (for squeue)
#SBATCH	--mem=30000

# load modules
module purge
module load intel
module load mpich
module load FFTW

# compiler and mpirun must be from same MPI implementation
/apps/software/impi/2017.1.132-iccifort-2017.1.132-GCC-5.4.0-2.26/bin64/mpirun -n NUMTASKS ./sc_compiled/flow36
