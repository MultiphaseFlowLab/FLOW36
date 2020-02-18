#Flow 36 GPU version

This version of Flow 36 uses the CUFFT library to perform the FFTs and the DCT, both in forward and backward directions, on the GPU.

The initialisation is based on the CPU so the result shall be fully compatible, except for the different approximation.

To build the GPU-version, just set the flag GPU_RUN to 1 in the compile.sh and select the machine GPU_local (#12)

#Notes:

For what concerns the DCT, the data is passed, prepared and transformed on the GPU; it has to be transposed once to achieve memory coalescent readings and a second time when passing it back to the host;

instead, the FFT in x and y are performed with the original memory layout.

#Warning: multi-GPU version

The initialization of the GPUs has not been tested for a multiple MPI task execution!!!
