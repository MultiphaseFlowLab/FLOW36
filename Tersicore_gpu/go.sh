NVARCH=Linux_x86_64; export NVARCH
NVCOMPILERS=/opt/nvidia/hpc_sdk; export NVCOMPILERS
MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/23.1/compilers/man; export MANPATH
PATH=$NVCOMPILERS/$NVARCH/23.1/compilers/bin:$PATH; export PATH
export PATH=$NVCOMPILERS/$NVARCH/23.1/comm_libs/mpi/bin:$PATH
export MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/23.1/comm_libs/mpi/man
mpirun -n NUMTASKS ./sc_compiled/flow36
