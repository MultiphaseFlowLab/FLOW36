#!/bin/bash
####### 48 cores x node (2 x Intel Xeon)#####
#SBATCH -J TEST
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH --time=00-00:05:00
#SBATCH --partition=skylake_0096
#SBATCH --qos=skylake_0096

# load modules
module purge
module load intel/19.1.3
module load intel-mpi/2019.10.317-intel-19.1.3.304-x276qb5
module load fftw


mpirun -n NUMTASKS ./sc_compiled/flow36
