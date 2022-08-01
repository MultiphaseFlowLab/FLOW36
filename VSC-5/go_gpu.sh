#!/bin/sh
#SBATCH -J "F36GPU"
#SBATCH -N 1
#SBATCH --partition=gpu_a100_dual
#SBATCH --qos goodluck
#SBATCH --gres=gpu:2
spack load nvhpc
mpirun -n NUMTASKS ./sc_compiled/flow36
