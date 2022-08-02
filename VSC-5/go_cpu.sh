#!/bin/sh
#SBATCH -J "F36CPU"
#SBATCH -N 1
#SBATCH --partition=zen3_0512
#SBATCH --qos goodluck
#SBATCH --time=00-00:03:00
spack unload --all
spack load gcc@11.2
spack load fftw@3.3.10%gcc@11.2.0
spack load /btrp7nc
mpirun -n NUMTASKS ./sc_compiled/flow36
