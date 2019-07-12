#!/bin/bash

#SBATCH -N1 -n4                             # n tasks on N node (total number of tasks)
#SBATCH --time=5:00                         # time limits: hh:mm:ss
#SBATCH --error error.err                   # std-error file
#SBATCH --output output.out                 # std-output file
#SBATCH --account=ac560tp                   # account number
#SBATCH --reservation=challenge
#SBATCH --job-name=name                     # job name (for squeue)

# load modules
module purge
module load pgi/19.4
module load mpi/pgi_openmpi/19.4
module load fftw3/3.3.4
#module load intel
#module load mpi/intel_mpi
#module load fftw3/3.3.4


mpirun -n NUMTASKS ./sc_compiled/flow36


# submit script with qsub

# ncpus : number of CPUs per node
# mpiprocs : number of MPI processes per node
# mem : RAM memory per node
