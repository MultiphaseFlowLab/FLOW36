#!/bin/bash
#####SBATCH --account="    "
#SBATCH --job-name="flo36_test"
#SBATCH --time=00:05:00
#SBATCH --nodes=1      ##adjust
#SBATCH --ntasks-per-node=128
#SBATCH --ntasks-per-core=1 # That guarantees every MPI taks will be bind to one CPU core which is very effective
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --partition=cn  # You should use this partition name of Discoverer

# load modules
module purge
module load nvidia
module load nvhpc-nompi/latest
module load gcc/11/latest
module load openmpi/4/nvidia/latest
module load fftw/3/latest-nvidia-openmpi

mpirun  ./sc_compiled/flow36 #not necessary to specify task on Discoverer
