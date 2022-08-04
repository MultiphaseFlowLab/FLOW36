#!/bin/bash
#SBATCH --account="project_465000168"
#SBATCH --job-name="flow36_test"
#SBATCH --time=00:05:00
#SBATCH --nodes=1      ##adjust
#SBATCH --ntasks-per-node=128
#SBATCH --ntasks-per-core=1
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --partition=small
#SBATCH --mem=40G

# load modules
module purge
module load PrgEnv-gnu
module load craype-x86-milan
module load cray-fftw

srun   ./sc_compiled/flow36
