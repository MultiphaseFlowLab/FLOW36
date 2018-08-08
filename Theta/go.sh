#!/bin/bash
#COBALT -t 30
#COBALT -n 2
#COBALT --attrs=mcdram=cache
#COBALT -A ATPESC2018
#COBALT -q training
#COBALT --jobname name

module load fftw
module load craype-hugepages16M


export n_nodes=$COBALT_JOBSIZE
export n_mpi_ranks_per_node=64
export n_mpi_ranks=$(($n_nodes * $n_mpi_ranks_per_node))
##export n_openmp_threads_per_rank=4
##export n_hyperthreads_per_core=2
##export n_hyperthreads_skipped_between_ranks=4
aprun -n $n_mpi_ranks -N $n_mpi_ranks_per_node  ./sc_compiled/flow36  # for cache mcdram
##aprun -n $n_mpi_ranks -N $n_mpi_ranks_per_node /usr/bin/numactl -m 1 ./sc_compiled/flow36  # for flat mcdram
