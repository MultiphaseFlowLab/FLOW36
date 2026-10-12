#!/bin/bash
#SBATCH --account="IscrB_EXCEED"
#SBATCH --job-name="flo36gpu_test"
#SBATCH --time=00:15:00   ## profiling needs time to write the reports after the run
#SBATCH --ntasks=NUMTASKS   ## total MPI tasks = NYCPU*NZCPU (set by compile.sh), nodes = ntasks/4
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4   ###4 GPUs per node on 4 MPI tasks
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --partition=boost_usr_prod
##SBATCH --qos=boost_qos_dbg   ## debug QOS: max 8 nodes, 30 min

# to avoid perl warning
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
# load modules
module purge
module load nvhpc/25.11
# MPI: HPC-X (CUDA-aware Open MPI) bundled with NVHPC 25.11, CUDA 12.9 build
NVHPC_DIR=$(dirname "$(readlink -f "$(which nvfortran)")")/../..
source "$NVHPC_DIR/comm_libs/12.9/hpcx/hpcx-2.25.1/hpcx-init.sh"
hpcx_load

# the code uses paths relative to set_run (./sc_compiled, ./results, ./initial_fields):
# run from set_run, also when the job is submitted from the main folder (sbatch set_run/go.sh)
cd "${SLURM_SUBMIT_DIR:-.}"
if [ -d ./set_run/sc_compiled ]; then cd ./set_run; fi

#if using HPC-SDK, CUDA-aware already enabled):
# 4 consecutive ranks per node (one y-group with NYCPU=4: xz2yz/yz2xz stay on NVLink);
# binder_leo.sh sets the InfiniBand card of each GPU (UCX_NET_DEVICES=mlx5_<local rank>:1)
chmod +x ./binder_leo.sh

# PROFILE=1: profile every rank with nsys from nvhpc/25.11, one report per rank (report_<rank>.nsys-rep).
# Low-overhead setup: only CUDA + NVTX traced (no MPI/OpenACC interception, no CPU sampling);
# MPI waits appear as gaps between kernels. Recording window: from PROFDELAY seconds after launch
# (startup + a few steps, see test.out) for PROFDUR seconds; analyse a step in the middle of the
# window (the first recorded step pays the start-up cost of the tracing).
PROFILE=0
PROFDELAY=30
PROFDUR=1

# MPI over UCX only (InfiniBand + GPUDirect RDMA): if UCX cannot start, the run stops with an
# error instead of silently falling back to TCP (ob1/tcp, ~0.3 GB/s between nodes)
MPIOPT="--mca pml ucx"
# UCX_INFO=1: print the UCX protocols selected for each message size and memory type (look for
# rc/dc_mlx5 zcopy/get/put on cuda memory = GPUDirect RDMA; cuda_copy/host staging = no GDR).
# The tables go to test.err/test.out; set UCX_INFO=0 for production runs.
UCX_INFO=0
if [ "$UCX_INFO" == "1" ]; then
  export UCX_PROTO_INFO=y
  MPIOPT="$MPIOPT -x UCX_PROTO_INFO=y"
fi

if [ "$PROFILE" == "1" ]; then
  # the binder sets UCX_NET_DEVICES and then starts nsys, which profiles flow36 directly
  mpirun -n NUMTASKS --map-by ppr:4:node $MPIOPT ./binder_leo.sh nsys profile -t cuda,nvtx --sample=none --cpuctxsw=none \
    --delay=$PROFDELAY --duration=$PROFDUR -o report_%q{OMPI_COMM_WORLD_RANK} ./sc_compiled/flow36
else
  mpirun -n NUMTASKS --map-by ppr:4:node $MPIOPT ./binder_leo.sh ./sc_compiled/flow36
fi

# submit script with sbatch
