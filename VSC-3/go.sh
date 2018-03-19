#!/bin/bash

#SBATCH -J rhor-
#SBATCH -N 32
#SBATCH --ntasks-per-node=32
#SBATCH --ntasks-per-core=2
#SBATCH --time=03-00:00:00
#SBATCH --no-requeue
####SBATCH --partition=mem_0128
####SBATCH --qos=devel_0128

####SBATCH --mail-type=END
####SBATCH --mail-user=aaa.mail.com

# load modules
module load intel/16 intel-mpi/5.1.3 fftw/3.3.4-DP
# or (problem with .mod files)
#module load gcc/5.3 intel-mpi/5.1.3 fftw/3.3.4-DP

mpirun -n NUMTASKS ./sc_compiled/flow36
