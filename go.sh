#!/bin/sh
#SBATCH -J "F36GPU"
#SBATCH -N 1
#SBATCH --partition=gpu_a100_dual
#SBATCH --qos goodluck
#SBATCH --gres=gpu:2
#SBATCH --time=00-00:03:00
spack load nvhpc@22.5
mpirun -n NUMTASKS ./sc_compiled/flow36
