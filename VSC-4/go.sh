#!/bin/bash
####### 48 cores x node (2 x Intel Xeon)#####
#SBATCH -J TEST
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
####### 13/12/2019 problems with HT (no scalability)
######SBATCH --ntasks-per-core=2 HT ON/OFF
#SBATCH --time=00-00:05:00
#SBATCH --partition=mem_0096
#SBATCH --qos=mem_0096

# load modules
module load intel/19.0.5
module load intel-mpi/2019.7.pre
module load fftw/3.3.8-intel-19.0.5.281-un2sutg
# or (problem with .mod files)
#module load gcc/5.3 intel-mpi/5.1.3 fftw/3.3.4-DP

mpirun -n NUMTASKS ./sc_compiled/flow36
