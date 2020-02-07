#Flow 36 GPU version

This version of Flow 36 implements the CUFFT library to perform the FFTs and the DCT, both in forward and backward directions, on the GPU.
The initialisation is based on the CPU so the result shall be fully compatible, with regard to the device tollerance.
To build the GPU-version, just set to one the flag GPU_RUN in the compile.sh
#Notes:
For what concerns the DCT, the data is passed and fully manipuled on the GPU; it has to be transposed once to achieve memory coalescent readings and complete the operations in a performing way;
instead, the FFT in x and y are performed with the original memory layout.
#Warning: multi-GPU version
The initialization of the GPUs has not been tested for a multiple MPI task execution!!!
