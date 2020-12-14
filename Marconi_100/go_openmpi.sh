#!/bin/bash

#SBATCH --account="IscrB_LUPIN"
#SBATCH --job-name="surf36test"
#SBATCH --time=00:05:00
#SBATCH --nodes=2                       ###modify according to necessity
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=32            ###modify according to necessity
#SBATCH --cpus-per-task=1
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --partition=m100_usr_prod

# load modules
module purge
module load profile/advanced
module load gnu
module load cuda/10.1
module load openmpi/4.0.3--gnu--8.4.0
module load fftw/3.3.8--gnu--8.4.0

mpirun -n 64 ./sc_compiled/flow36


# submit script with sbatch
