#!/bin/bash
#SBATCH --account="Ppp4x_D069"
#SBATCH --job-name="flo36gpu_test"
#SBATCH --time=00:05:00
#SBATCH --nodes=1      ##adjust
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4   ###4 GPUs per node on 4 MPI tasks
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --partition=m100_usr_prod

# to avoid perl warning
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
# load modules
module purge
module load hpc-sdk

#if using HPC-SDK (OPENPMPI) use (CUDA-aware already enabled):
mpirun -n NUMTASKS --map-by socket ./sc_compiled/flow36

# submit script with sbatch
