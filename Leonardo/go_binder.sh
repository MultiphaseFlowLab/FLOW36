#!/bin/bash
#SBATCH --account="IscrB_BMO"
#SBATCH --job-name="flo36gpu_test"
#SBATCH --time=00:25:00
#SBATCH --nodes=1      ##adjust
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
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
chmod 777 binder.sh
mpirun -n NUMTASKS --map-by node:PE=8 --rank-by core  ./binder.sh ./sc_compiled/flow36

# submit script with sbatch
