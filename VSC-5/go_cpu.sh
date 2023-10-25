#!/bin/sh
#SBATCH -J  ^=  
#SBATCH -N 1
#SBATCH --partition=zen3_0512
#SBATCH --qos zen3_0512
#SBATCH --time=00-00:03:00

spack load gcc@12.2.0
#spack load nvhpc@22.9
spack load mpich@4.1.1
spack load /xgvooar

mpirun -n NUMTASKS ./sc_compiled/flow36