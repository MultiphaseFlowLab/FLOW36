#ifndef CUDA_PHYS_SPEC_FG
#define CUDA_PHYS_SPEC_FG

extern "C" void h_fftxfwd_fg(double *in_r, double *out_r, double *out_c, int aliasing);

extern "C" void h_chebfwd_fg(double *in_r, double *in_c, double *out_r, double *out_c, int aliasing);

extern "C" void h_fftymanyfwd_fg(double *in_r, double *in_c, double *out_r, double *out_c, int aliasing);

#endif

