#!/bin/sh
#SBATCH -J "F36CPU"
#SBATCH -N 1
#SBATCH --partition=zen3_0512
#SBATCH --qos goodluck
#SBATCH --time=00-00:03:00
pack unload --all
spack load gcc@11.2
spack load /ltuvxjf  #FFTW3
mpirun -n NUMTASKS ./sc_compiled/flow36
