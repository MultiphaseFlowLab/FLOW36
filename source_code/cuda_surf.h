#ifndef CUDA_SURF
#define CUDA_SURF

//__global__ void k_mirror_c(int N, cufftDoubleReal *in, cufftDoubleComplex *out);

__global__ void k_sec_copy(cufftDoubleComplex *out, cufftDoubleComplex *im, int size);

__global__ void k_mirror(int N, cufftDoubleReal *in, cufftDoubleReal *out);

__global__ void manipulate(cufftDoubleComplex *a, int nz);

__global__ void saxpy(cufftDoubleComplex *a);

#endif
