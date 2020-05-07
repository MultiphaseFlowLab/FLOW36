#Flow 36 GPU version

This version of Flow 36 uses the CUFFT library to perform the FFTs and the DCT, both in forward and backward directions, on the GPU.

The initialisation is based on the CPU so the result shall be fully compatible, except for the different approximation.

To build the GPU version, set the variable GPU_RUN to 1 in the compile.sh.

Marconi 100:
the file go.sh is not automatically compiled, not straightforward!
Each node has 4 Tesla V100 GPUs. Set number of tasks per node and number of GPUs according to your necessities (always the same in pure flow solver)
Modify also the mpi wrapper creation (mpirun -np ..)

Launch with SBATCH go.sh

#Notes:

For what concerns the DCT, the data is passed, prepared and transformed on the GPU; it has to be transposed once to achieve memory coalescent readings and a second time when passing it back to the host;

instead, the FFT in x and y are performed with the original memory layout.

#Warning: multi-GPU version

output is valid after one timestep up to 4 MPI tasks (2x2)
