#ifndef CUDA_SPEC_PHYS
#define CUDA_SPEC_PHYS

void cuda_spectral_to_phys(double *real_d,double *imm_d,double *out_d);

extern "C" void h_ffty_many_bwd(double *in_r, double *in_c, double *out_r, double *out_c, int aliasing);

extern "C" void h_chebyshev_back(double *in_r, double *in_c, double *out_r, double *out_c, int aliasing);

extern "C" void h_ffty_back(double *in_r, double *in_c, double *out_r, double *out_c, int aliasing);

extern "C" void h_fftx_back(double *in_r, double *in_c, double *out_r, int aliasing);

#endif
