#!/bin/sh
#SBATCH -J "F36CPU"
#SBATCH -N 1
#SBATCH --partition=zen3_0512
#SBATCH --qos goodluck

mpirun -n NUMTASKS ./sc_compiled/flow36
