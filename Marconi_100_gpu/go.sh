#!/bin/bash

#SBATCH --account="IscrB_BLUES"
#SBATCH --job-name="flo36gpu_test"
#SBATCH --time=00:05:00
#SBATCH --nodes=1
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
module load profile/advanced
module load hpc-sdk/2021--binary
module load cuda/11.3

mpirun -n NUMTASKS ./sc_compiled/flow36


# submit script with sbatch