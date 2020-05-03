#ifndef CUDA_SPEC_PHYS_FG
#define CUDA_SPEC_PHYS_FG

extern "C" void h_chebback_fg(double *in_r, double *in_c, double *out_r, double *out_c, int aliasing);

extern "C" void h_fftxback_fg(double *in_r, double *in_c, double *out_r, int aliasing);

extern "C" void h_fftymanybwd_fg(double *in_r, double *in_c, double *out_r, double *out_c, int aliasing);


#endif
