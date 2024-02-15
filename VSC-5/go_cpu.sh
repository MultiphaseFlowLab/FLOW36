#!/bin/sh
#
#SBATCH	--mail-type=ALL
#SBATCH --mail-user=davide.procacci@tuwien.ac.at
#SBATCH -J "F36CPU"
#SBATCH -N 1
#SBATCH --ntasks-per-node=NUMTASKS
#SBATCH --mem=60GB
#SBATCH --partition=zen3_0512
#SBATCH --qos zen3_0512
##SBATCH --qos zen3_0512_devel
#SBATCH --time=03-00:00:00

spack unload --all
spack load gcc@12.2.0
spack load mpich@4.1.1
spack load /xgvooar
module list
export LD_LIBRARY_PATH=$LIBRARY_PATH:$LD_LIBRARY_PATH
#module purge
#module load --auto openmpi/4.1.4-gcc-12.2.0-sugs3ze
#module load --auto fftw/3.3.10-gcc-12.2.0-42q2cmu
#module list

mpirun -np NUMTASKS ./sc_compiled/flow36
