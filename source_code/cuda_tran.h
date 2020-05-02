#ifndef CUDA_TRAN
#define CUDA_TRAN

void __global__ k_alias_1st_cmp(cufftDoubleComplex *a, int al_low, int nx, int dim);

void __global__ k_alias_1st(double *a, double *b, int al_low, int nx, int dim);

__global__ void k_merge_cmp(cufftDoubleComplex *out, double *re, double *im, int size);

__global__ void k_sec_separate(cufftDoubleComplex *in,
		                       double *re, double *im,
		                       int size);

__global__ void k_mixed_prod(double *a, double *b, double *out, int size);

__device__ double atomic_max(double* address, double val);

__global__ void k_courant(double *u, double *v, double *w, double *dz_vec, double *result,
		                  double dx, double dy, double dt,
		                  int size, int fstart, int nx, int spy);

__global__ void k_norm_real(double *a,
		                    double n, int size);

__global__ void k_norm_cmp(cufftDoubleComplex *a,
		                    double n, int size);

__global__ void k_sep_cmp(cufftDoubleComplex *in, double *out, double norm, int size);

__global__ void k_sec_copy(cufftDoubleComplex *out,
		                   cufftDoubleComplex *im,
		                   int size);

__global__ void k_manip(double *a, int nz, int dim);

__global__ void k_manip_cmp (cufftDoubleComplex *a, int nz, int size);

__global__ void k_mirr_bigtime(double *out1, double *out2,
		                       double *in1,  double *in2,
                               int nx, int nx_big, int size);

__global__ void k_mirror(double *out,
		                 double *in,
		                 int nx, int nx_big);

__global__ void k_mirror_2n(double *out,
		                    double *in,
                            int nx, int nx_big);

void __global__ k_alias_small(double *a, int ali, int nx);

void __global__ k_cmp_alias_small(cufftDoubleComplex *a, int ali, int nx);

void __global__ k_cmp_alias(cufftDoubleComplex *a, int al_low, int al_val, int nx, int tot_size);

void __global__ cmp_alias(cufftDoubleComplex *a, int ali, int nx);

__global__ void k_t102(double *out,
                       double *in,
                       int nx, int ny, int nz);

__global__ void k_cmp_t102(cufftDoubleComplex *out,
                           cufftDoubleComplex *in,
                           int nx, int ny, int nz);

__global__ void k_cmp_t210(cufftDoubleComplex *out,
		                   cufftDoubleComplex *in,
		                   int nx, int ny, int nz);

__global__ void k_cmp_t021(cufftDoubleComplex *out,
                           cufftDoubleComplex *in,
                           int nx, int ny, int nz);

#endif
