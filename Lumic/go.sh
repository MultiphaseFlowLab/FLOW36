#!/bin/bash
#####SBATCH --account="    "
#SBATCH --job-name="flo36_test"
#SBATCH --time=00:05:00
#SBATCH --nodes=1      ##adjust
#SBATCH --ntasks-per-node=128
#SBATCH --ntasks-per-core=1 
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --partition=

# load modules
module purge
module load PrgEnv-cray
module load craype-x86-milan
module load cray-fftw

mpirun  -np NUMTASKS ./sc_compiled/flow36 #not necessary to specify task on Discoverer
