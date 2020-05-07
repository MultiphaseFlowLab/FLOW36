#!/bin/bash
#SBATCH --account="IscrB_LUPIN"
#SBATCH --partition=m100_usr_prod
#SBATCH --job-name="surf36test"
#SBATCH --time=00:05:00
#SBATCH --nodes=1			###modify according to necessity
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=4		###modify according to necessity
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:4			###same as ntasks per node
#SBATCH --error="err.err"
#SBATCH --output="out.out"
#SBATCH --qos=m100_qos_dbg    


module load profile/base   
module load gnu/8.4.0   
module load cuda/10.1   
module load spectrum_mpi/10.3.1--binary   
module load fftw/3.3.8--spectrum_mpi--10.3.1--binary

##module load profile/advanced   
##module load gnu/8.4.0   
##module load cuda/10.1   
##module load openmpi/4.0.3--gnu--8.4.0   
##module load fftw/3.3.8--gnu--8.4.0  


##export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
##export CRAY_CUDA_MPS=1
##export MPICH_RDMA_ENABLED_CUDA=1
#srun ./sc_compiled/surf36_gpu.exe > a.out
##mpirun -n 4 --map-by socket ./surf_gpu.exe > a.out

mpirun -np 4 --map-by socket ./sc_compiled/surf36_gpu.exe > a.out		###modify according to necessity

##SPECTRUM MPI NOTE
##not srun but:
##mpirun -n 4 --map-by socket:PE=4 ./myprogram
##CUDA AWARE MPI flag -gpu
