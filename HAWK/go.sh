# submit script with "qsub -koed go.sh" to have STDOUT and STDERR in $PBS_O_WORKDIR at runtime. Otherwise they are copied there only at the end o the job
# guide at https://www.altair.com/pdfs/pbsworks/PBSUserGuide19.2.3.pdf

#!/bin/bash
# AMD Rome Epyc: 2X AMD Epyc Rome 7742, 64 cores/socket, 128 cores/node, support hyperthreading 2x

# select=<number of nodes>:<node_resource_variable=type>
#PBS -l select=1:node_type=rome:node_type_mem=256gb:mpiprocs=128:ncpus=128:mem=250gb
#PBS -l walltime=00:20:00
#PBS -o output.out
#PBS -e error.err
#PBS -N job_name
# queues: single (single node), normal (2 to 4096 nodes), test (1 to 384 nodes, max walltime 25 mins), interactive
###PBS -q test

# Change to the direcotry that the job was submitted from
cd $PBS_O_WORKDIR

module purge
# load modules
module load gcc
module load openmpi
# check aocl for AMD-optimized FFTW, see https://kb.hlrs.de/platforms/index.php/Libraries(Hawk)
module load fftw

# Launch the parallel mpi application (compiled with intel mpi) to the allocated compute nodes
mpirun -np NUMTASKS ./sc_compiled/flow36


# ncpus : number of CPUs per node
# mpiprocs : number of MPI processes per node
# mem : RAM memory per node
