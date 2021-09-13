#!/bin/bash

#SBATCH --account="IscrC_INFLO"
#SBATCH --job-name="march"
#SBATCH --time=00:01:00
#SBATCH --nodes=1                      ###modify according to necessity
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=24            ###modify according to necessity
#SBATCH --cpus-per-task=1
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --partition="g100_usr_prod"

# to avoid perl warning
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
# load modules
module purge
module load profile/advanced
module load profile/base
module load intel/oneapi-2021--binary
module load intelmpi/oneapi-2021--binary
module load fftw/3.3.9--intelmpi--oneapi-2021--binary

srun -n NUMTASKS ./sc_compiled/flow36


# submit script with sbatch
