#!/bin/bash
#SBATCH --account="r2025r02-287-users"
#SBATCH --job-name="flo36gpu_test"
#SBATCH --time=00:15:00
#SBATCH --ntasks=NUMTASKS   ## total MPI tasks = NYCPU*NZCPU (set by compile.sh), nodes = ntasks/4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:4   ###4 GPUs per node on 4 MPI tasks
#SBATCH --output=test_%j.out
#SBATCH --error=test_%j.err
#SBATCH --partition=gpu
# Vega charges whole nodes (1 node hour = 128 core hours): always use 4 ranks (4 GPUs) per node
# setup as in MHIT36 (channel/go_vega.sh), which runs on Vega

# to avoid perl warning
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
# modules: NVHPC 24.9 with HPC-X 2.20 (CUDA-aware Open MPI), CUDA 12.6, loaded only if the
# job did not get them from the login environment (the job environment is not always exported)
echo "$(date) loading modules"
if ! mpif90 --version 2>&1 | grep -qi nvfortran; then
  module purge
  module load NVHPC-HPCX-CUDA12/24.9 NVHPC/24.9
fi
echo "$(date) modules done"
if ! mpif90 --version 2>&1 | grep -qi nvfortran; then
  echo "NVHPC not loaded (module load NVHPC-HPCX-CUDA12/24.9 NVHPC/24.9 failed)"
  exit 1
fi

# the code uses paths relative to set_run (./sc_compiled, ./results, ./initial_fields):
# run from set_run, also when the job is submitted from the main folder (sbatch set_run/go.sh)
cd "${SLURM_SUBMIT_DIR:-.}"
if [ -d ./set_run/sc_compiled ]; then cd ./set_run; fi

# the GPU of each local rank is set in main.f90 (acc_set_device_num);
# binder_vega.sh sets the InfiniBand card of each rank
chmod +x ./binder_vega.sh

# PROFILE=1: profile every rank with nsys, one report per rank (report_<rank>.nsys-rep)
PROFILE=0
PROFDELAY=30
PROFDUR=1

# LAUNCHER=mpirun: HPC-X mpirun as in MHIT36 (default)
# LAUNCHER=srun: Slurm starts the ranks (PMIx)
LAUNCHER=mpirun

# extra mpirun options, e.g. "--mca pml ucx" (MPI over UCX only, no TCP fallback)
MPIOPT=""
# UCX_INFO=1: print the UCX protocols selected for each message size and memory type
UCX_INFO=0
if [ "$UCX_INFO" == "1" ]; then
  export UCX_PROTO_INFO=y
  MPIOPT="$MPIOPT -x UCX_PROTO_INFO=y"
fi

# diagnostics: printed before the launch, so an empty log means the script never started
echo "Job $SLURM_JOB_ID on $SLURM_JOB_NODELIST, $SLURM_NTASKS tasks, launcher $LAUNCHER, dir $(pwd)"
echo "mpirun: $(which mpirun)"
nvidia-smi -L
ls -l ./sc_compiled/flow36

if [ "$PROFILE" == "1" ]; then
  PROFCMD="nsys profile -t cuda,nvtx --sample=none --cpuctxsw=none --delay=$PROFDELAY --duration=$PROFDUR -o report_%q{OMPI_COMM_WORLD_RANK}"
else
  PROFCMD=""
fi

if [ "$LAUNCHER" == "srun" ]; then
  srun --mpi=pmix --cpus-per-task=8 --cpu-bind=cores ./binder_vega.sh $PROFCMD ./sc_compiled/flow36
else
  # --bind-to none: Open MPI's own core binding picked cores outside the Slurm allocation
  # (hwloc_set_cpubind error on gn27); the Slurm cgroup already confines the ranks
  mpirun -np NUMTASKS --rank-by core --bind-to none $MPIOPT ./binder_vega.sh $PROFCMD ./sc_compiled/flow36
fi

# submit script with sbatch
