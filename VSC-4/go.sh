#!/bin/bash
####### 96 cores x node (2 x Intel)
#SBATCH -J TEST
#SBATCH -N 1
#SBATCH --ntasks-per-node=96
######SBATCH --ntasks-per-core=2 HT ON/OFF
#SBATCH --time=00-00:05:00
#SBATCH --partition=mem_0096
#SBATCH --qos=mem_0096

# load modules
module load intel
module load intel-mpi
module load fftw
# or (problem with .mod files)
#module load gcc/5.3 intel-mpi/5.1.3 fftw/3.3.4-DP

mpirun -n NUMTASKS ./sc_compiled/flow36
