#!/bin/sh
#SBATCH -J "F36CPU"
#SBATCH -N 1
#SBATCH --partition=zen3_0512
#SBATCH --qos zen3_0512
#SBATCH --time=00-00:03:00
pack unload --all
module load openmpi/4.1.4-gcc-12.2.0-sugs3ze
module load fftw/3.3.10-gcc-12.2.0-42q2cmu
mpirun -n NUMTASKS ./sc_compiled/flow36
