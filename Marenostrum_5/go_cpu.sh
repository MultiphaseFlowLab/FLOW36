#!/bin/bash
#SBATCH --account="ehpc87"
#SBATCH --job-name="flow36"
#SBATCH --time=00:05:00
#SBATCH --ntasks=112
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --qos=gp_debug

##
module purge
module load gcc/13.2.0
module load openmpi/4.1.5-gcc
module load fftw/3.3.10-gcc-ompi

#automatically set from tasks
srun ./sc_compiled/flow36

# submit script with sbatch
