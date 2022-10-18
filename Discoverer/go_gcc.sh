#!/bin/bash
#####SBATCH --account="    "
#SBATCH --job-name="Pr4_gcc"
#SBATCH --time=24:00:00
#SBATCH --nodes=1      ##adjust
#SBATCH --ntasks-per-node=128
#SBATCH --ntasks-per-core=1 # That guarantees every MPI taks will be bind to one CPU core which is very effective
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --partition=cn  # You should use this partition name of Discoverer

# load modules
module purge
module load gcc/11/latest
module load mpich/3/gcc/latest
module load fftw/3/latest-gcc-mpich

mpirun  ./sc_compiled/flow36 #not necessary to specify task on Discoverer
