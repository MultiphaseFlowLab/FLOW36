#!/bin/bash
#SBATCH --account="IscrB_TORNADO"
#SBATCH --job-name="flo36gpu_test"
#SBATCH --time=00:25:00
#SBATCH --nodes=1      ##adjust
#SBATCH --ntasks-per-node=32
#SBATCH --output=test.out
#SBATCH -p boost_usr_prod
#SBATCH --error=test.err

# to avoid perl warning
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
# load modules
module purge
module load gcc/11.3.0
module load openmpi/4.1.4--gcc--11.3.0-cuda-11.8
module load fftw/3.3.10--openmpi--4.1.4--gcc--11.3.0

#if using HPC-SDK (OPENPMPI) use (CUDA-aware already enabled):
mpirun -n NUMTASKS ./sc_compiled/flow36

# submit script with sbatch
