#!/bin/bash

#SBATCH --account="IscrB_LUPIN"
#SBATCH --job-name="surf36test"
#SBATCH --time=00:05:00
#SBATCH --nodes=1                       ###modify according to necessity
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=32            ###modify according to necessity
#SBATCH --cpus-per-task=1
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --partition=m100_usr_prod

# to avoid perl warning
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
# load modules
module purge
module load profile/advanced
module load profile/base
module load gnu/8.4.0
module load cuda/10.1
module load spectrum_mpi/10.3.1--binary
module load fftw/3.3.8--spectrum_mpi--10.3.1--binary

mpirun -n 32 ./sc_compiled/flow36


# submit script with sbatch
