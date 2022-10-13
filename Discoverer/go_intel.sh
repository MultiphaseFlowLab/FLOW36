#!/bin/bash
#####SBATCH --account="    "
#SBATCH --job-name="intel_test"
#SBATCH --time=00:05:00
#SBATCH --nodes=1      ##adjust
#SBATCH --ntasks-per-node=128
#SBATCH --ntasks-per-core=1 # That guarantees every MPI taks will be bind to one CPU core which is very effective
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --partition=cn  # You should use this partition name of Discoverer

#Use Infiniband for node-2-node
export UCX_NET_DEVICES=mlx5_0:1


# load modules
module purge
module load intel compiler-rt/latest
module load compiler/latest
module load mpi/latest
module load mkl/latest

mpirun  ./sc_compiled/flow36 #not necessary to specify task on Discoverer
