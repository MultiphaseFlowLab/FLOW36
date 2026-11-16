#!/bin/bash
# Marge (local cluster): one node, 2x RTX PRO 6000 Blackwell -> NYCPU*NZCPU = 2 (one rank per GPU)
# run from set_run: ./go.sh  (or nohup ./go.sh > test.out 2> test.err &)
# the #SBATCH lines are used only if Marge has Slurm (sbatch go.sh), bash ignores them
#SBATCH --job-name="flo36gpu"
#SBATCH --ntasks=NUMTASKS
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:2
#SBATCH --output=test_%j.out
#SBATCH --error=test_%j.err

# to avoid perl warning
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

# same environment as compile.sh
# NVHPC 26.5 from the HPC SDK install, with its bundled MPI (comm_libs/mpi = HPC-X Open MPI)
NVHPC_VERSION=26.5
NVARCH=Linux_x86_64; export NVARCH
NVCOMPILERS=/opt/nvidia/hpc_sdk; export NVCOMPILERS
export MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/$NVHPC_VERSION/compilers/man
export PATH=$NVCOMPILERS/$NVARCH/$NVHPC_VERSION/compilers/bin:$PATH
export PATH=$NVCOMPILERS/$NVARCH/$NVHPC_VERSION/comm_libs/mpi/bin:$PATH
export MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/$NVHPC_VERSION/comm_libs/mpi/man
if ! mpif90 --version 2>&1 | grep -qi nvfortran; then
  echo "mpif90 ($(which mpif90)) is not the NVHPC one: check NVCOMPILERS and NVHPC_VERSION"
  exit 1
fi

# the code uses paths relative to set_run (./sc_compiled, ./results, ./initial_fields)
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
if [ -d ./set_run/sc_compiled ]; then cd ./set_run; fi

# the GPU of each rank is set in main.f90 (acc_set_device_num, local rank -> device)
chmod +x ./binder_marge.sh

# PROFILE=1: profile every rank with nsys, one report per rank (report_<rank>.nsys-rep)
PROFILE=0
PROFDELAY=30
PROFDUR=1

echo "Marge run, NUMTASKS tasks, dir $(pwd), host $(hostname)"
echo "mpirun: $(which mpirun)"
nvidia-smi -L

if [ "$PROFILE" == "1" ]; then
  PROFCMD="nsys profile -t cuda,nvtx --sample=none --cpuctxsw=none --delay=$PROFDELAY --duration=$PROFDUR -o report_%q{OMPI_COMM_WORLD_RANK}"
else
  PROFCMD=""
fi

mpirun -np NUMTASKS --bind-to socket ./binder_marge.sh $PROFCMD ./sc_compiled/flow36
