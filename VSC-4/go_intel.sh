#!/bin/bash
####### 48 cores x node (2 x Intel Xeon)#####
#SBATCH -J TEST
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH --time=00-00:05:00
#SBATCH --partition=skylake_0096
#SBATCH --qos=skylake_0096

# load modules
module load intel
module load intel-mpi
module load mkl

mpirun -n NUMTASKS ./sc_compiled/flow36
