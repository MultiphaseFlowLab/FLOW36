#!/bin/bash
#SBATCH --account="r2025r02-287-users"
#SBATCH --job-name="flo36gpu_test"
#SBATCH --time=00:15:00
#SBATCH --ntasks=NUMTASKS   ## total MPI tasks = NYCPU*NZCPU (set by compile.sh), nodes = ntasks/4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=32   ## 128 cores per node / 4 ranks
#SBATCH --gres=gpu:4   ###4 GPUs per node on 4 MPI tasks
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --partition=gpu
# Vega charges whole nodes (1 node hour = 128 core hours): always use 4 ranks (4 GPUs) per node

# to avoid perl warning
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
# load modules: NVHPC 24.9 with HPC-X 2.20 (CUDA-aware Open MPI), CUDA 12.6
module purge
module load NVHPC-HPCX-CUDA12/24.9

# the code uses paths relative to set_run (./sc_compiled, ./results, ./initial_fields):
# run from set_run, also when the job is submitted from the main folder (sbatch set_run/go.sh)
cd "${SLURM_SUBMIT_DIR:-.}"
if [ -d ./set_run/sc_compiled ]; then cd ./set_run; fi

# 4 consecutive ranks per node (one y-group with NYCPU=4: xz2yz/yz2xz stay on the node);
# the GPU of each local rank is set in main.f90 (acc_set_device_num);
# binder_vega.sh only sets the InfiniBand card of each rank (when configured)
chmod +x ./binder_vega.sh

# PROFILE=1: profile every rank with nsys, one report per rank (report_<rank>.nsys-rep)
PROFILE=0
PROFDELAY=30
PROFDUR=1

# MPI over UCX only: if UCX cannot start, stop with an error instead of falling back to TCP
MPIOPT="--mca pml ucx"
# UCX_INFO=1: print the UCX protocols selected for each message size and memory type
UCX_INFO=0
if [ "$UCX_INFO" == "1" ]; then
  export UCX_PROTO_INFO=y
  MPIOPT="$MPIOPT -x UCX_PROTO_INFO=y"
fi

# LAUNCHER=srun: Slurm starts the ranks (PMIx), no mpirun daemons (default on Vega)
# LAUNCHER=mpirun: HPC-X mpirun
LAUNCHER=srun

# diagnostics: printed before the launch, so an empty log means the script never started
echo "Job $SLURM_JOB_ID on $SLURM_JOB_NODELIST, $SLURM_NTASKS tasks, launcher $LAUNCHER, dir $(pwd)"
echo "mpirun: $(which mpirun)"
srun --mpi=list 2>&1 | tr '\n' ' '; echo
nvidia-smi -L
ls -l ./sc_compiled/flow36

if [ "$PROFILE" == "1" ]; then
  PROFCMD="nsys profile -t cuda,nvtx --sample=none --cpuctxsw=none --delay=$PROFDELAY --duration=$PROFDUR -o report_%q{SLURM_PROCID}"
else
  PROFCMD=""
fi

if [ "$LAUNCHER" == "srun" ]; then
  srun --mpi=pmix --cpus-per-task=32 --cpu-bind=cores ./binder_vega.sh $PROFCMD ./sc_compiled/flow36
else
  mpirun -n NUMTASKS --map-by ppr:4:node:PE=32 $MPIOPT ./binder_vega.sh $PROFCMD ./sc_compiled/flow36
fi

# submit script with sbatch
