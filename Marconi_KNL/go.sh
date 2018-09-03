#!/bin/bash

#SBATCH -N1 -n1                             # n tasks on N node (total number of tasks)
#SBATCH --time=5:00                         # time limits: hh:mm:ss
#SBATCH --error error.err                   # std-error file
#SBATCH --output output.out                 # std-output file
#SBATCH --account=IscrB_TURBINLA            # account number
#SBATCH --partition=knl_usr_dbg             # partition to be used, queue dbg: debug, prod, bprod: bigprod
#SBATCH --job-name=name                     # job name (for squeue)
#SBATCH	--mem=30000

# load modules
module purge
module load env-knl
module load profile/global
module load intel/pe-xe-2017--binary
module load intelmpi/2017--binary
module load fftw/3.3.5--intelmpi--2017--binary
# or
#module load gnu/6.1.0
#module load openmpi/1-10.3--gnu--6.1.0
#module load fftw/3.3.4--openmpi--1-10.3--gnu--6.1.0

mpirun -n NUMTASKS ./sc_compiled/flow36


# submit script with qsub

# ncpus : number of CPUs per node
# mpiprocs : number of MPI processes per node
# mem : RAM memory per node
