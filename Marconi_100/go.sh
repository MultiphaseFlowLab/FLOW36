#!/bin/bash

#SBATCH --account="......."
#SBATCH --job-name="test"
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --partition=m100_usr_prod

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
