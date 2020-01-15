#!/bin/bash

#SBATCH -N1 -n1                             # n tasks on N node (total number of tasks)
#SBATCH --time=5:00                         # time limits: hh:mm:ss
#SBATCH --error error.err                   # std-error file
#SBATCH --output output.out                 # std-output file
#SBATCH --account=try19_soligo            # account number
#SBATCH --partition=dvd_usr_prod             # partition to be used
#SBATCH --job-name=name                     # job name (for squeue)
#SBATCH --mem=3000

# load modules
module purge
module load gnu
module load openmpi
module load fftw

mpirun -n NUMTASKS ./sc_compiled/flow36


# submit script with qsub

# ncpus : number of CPUs per node
# mpiprocs : number of MPI processes per node
# mem : RAM memory per node
