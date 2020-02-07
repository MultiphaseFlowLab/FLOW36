#ifndef CUDA_PHYS_SPEC
#define CUDA_PHYS_SPEC

void cuda_phys_to_spectral(double *d_uur, cufftDoubleComplex *d_batch);

extern "C" void h_ffty_many_fwd(double *in_r, double *in_c, double *out_r, double *out_c, int aliasing);

extern "C" void h_fftx_fwd(double *in_r, double *out_r, double *out_c, int aliasing);

extern "C" void h_ffty_fwd(double *in_r, double *in_c, double *out_r, double *out_c, int aliasing);

extern "C" void h_chebyshev_fwd(double *in_r, double *in_c, double *out_r, double *out_c, int aliasing);

#endif

