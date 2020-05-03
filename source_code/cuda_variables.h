#include <complex.h>
#include <cufft.h>
#include <stdio.h>

#ifndef CUDA_VARIABLES
#define CUDA_VARIABLES 1

#ifndef VAR_DECLS
# define _DECL extern
# define _INIT(x)
#else
# define _DECL
# define _INIT(x) = x
#endif


//GPU version variables-------------------------------------------------------------
//integers
_DECL int spx      _INIT(0);
_DECL int nsx      _INIT(0);
_DECL int nx       _INIT(0);
_DECL int npx      _INIT(0);
_DECL int spy      _INIT(0);
_DECL int fpy      _INIT(0);
_DECL int npy      _INIT(0);
_DECL int ny       _INIT(0);
_DECL int nz       _INIT(0);
_DECL int fpz      _INIT(0);
_DECL int npz      _INIT(0);

_DECL int aliasing _INIT(0);
_DECL int dblk     _INIT(0);
_DECL int ok       _INIT(0);
_DECL int nbx      _INIT(0);
_DECL int nby      _INIT(0);

//fg dims
_DECL int npsix    _INIT(0);
_DECL int fpypsi   _INIT(0);
_DECL int fpzpsi   _INIT(0);

_DECL int spxpsi   _INIT(0);
_DECL int spypsi   _INIT(0);
_DECL int npsiz    _INIT(0);

_DECL int npsiy    _INIT(0);


//doubles
_DECL double xl       _INIT(0);
_DECL double yl       _INIT(0);
_DECL double dx       _INIT(0);
_DECL double dy       _INIT(0);
_DECL double cou_val  _INIT(0);
_DECL double dt       _INIT(0);

//memory check flag
_DECL int check_mem _INIT(0);

//mpi comm variables
_DECL int nGPUs _INIT(0);
_DECL int idGPU _INIT(0);
_DECL int ilGPU _INIT(0);




//HOST ARRAYS-----------------------------------------------------------------------
//integer HOST pointers
//_DECL int *fstart    _INIT(NULL);
//
////double HOST pointers
//_DECL double *ur     _INIT(NULL);
//_DECL double *uc     _INIT(NULL);
//_DECL double *h_uout _INIT(NULL);
//
//_DECL double *vr     _INIT(NULL);
//_DECL double *vc     _INIT(NULL);
//_DECL double *h_vout _INIT(NULL);
//
//_DECL double *wr     _INIT(NULL);
//_DECL double *wc     _INIT(NULL);
//_DECL double *h_wout _INIT(NULL);
//
//_DECL double *uu     _INIT(NULL);
//
//_DECL double *z      _INIT(NULL);
//_DECL double *dz_h   _INIT(NULL);

//mixed products OUTPUT
//_DECL double *h_uur _INIT(NULL);
//_DECL double *h_uuc _INIT(NULL);
//
//_DECL double *h_uvr _INIT(NULL);
//_DECL double *h_uvc _INIT(NULL);
//
//_DECL double *h_uwr _INIT(NULL);
//_DECL double *h_uwc _INIT(NULL);
//
//_DECL double *h_vvr _INIT(NULL);
//_DECL double *h_vvc _INIT(NULL);
//
//_DECL double *h_vwr _INIT(NULL);
//_DECL double *h_vwc _INIT(NULL);
//
//_DECL double *h_wwr _INIT(NULL);
//_DECL double *h_wwc _INIT(NULL);




//GPU ARRAYS------------------------------------------------------------------------
//GPU integer arrays
//_DECL int *fstart_d  _INIT(NULL);

//GPU realpr arrays
_DECL double *res    _INIT(NULL);

//velocities
_DECL double *ur_d   _INIT(NULL);
_DECL double *uc_d   _INIT(NULL);
_DECL double *d_uout _INIT(NULL);

_DECL double *psir_d   _INIT(NULL);
_DECL double *psic_d   _INIT(NULL);

//_DECL double *vr_d   _INIT(NULL);
//_DECL double *vc_d   _INIT(NULL);
//_DECL double *d_vout _INIT(NULL);

//_DECL double *wr_d   _INIT(NULL);
//_DECL double *wc_d   _INIT(NULL);
//_DECL double *d_wout _INIT(NULL);




//_DECL double *test_r _INIT(NULL);
//_DECL cufftDoubleComplex *d_uuc2 _INIT(NULL);

//_DECL double *dz_d    _INIT(NULL);

//mixed product terms
//_DECL double *d_uu  _INIT(NULL);
//_DECL double *d_uv  _INIT(NULL);
//_DECL double *d_uw  _INIT(NULL);
//_DECL double *d_vv  _INIT(NULL);
//_DECL double *d_vw  _INIT(NULL);
//_DECL double *d_ww  _INIT(NULL);

//mixed products in spectral space
//_DECL cufftDoubleComplex *d_uuc _INIT(NULL);
//_DECL cufftDoubleComplex *d_uvc _INIT(NULL);
//_DECL cufftDoubleComplex *d_uwc _INIT(NULL);
//_DECL cufftDoubleComplex *d_vvc _INIT(NULL);
//_DECL cufftDoubleComplex *d_vwc _INIT(NULL);
//_DECL cufftDoubleComplex *d_wwc _INIT(NULL);



//subroutine operational arrays-----------------------------------------------------------------
//double arrays
//_DECL double *d_urs   _INIT(NULL);
//_DECL double *d_ucs   _INIT(NULL);
_DECL double *d_uopr  _INIT(NULL);
_DECL double *d_uopc  _INIT(NULL);

//fg
_DECL double *d_psir  _INIT(NULL);
_DECL double *d_psic  _INIT(NULL);

//_DECL double *d_uext  _INIT(NULL);

//cufftDoubleComplex arrays
_DECL cufftDoubleComplex *d_batch   _INIT(NULL);
_DECL cufftDoubleComplex *d_batch_c _INIT(NULL);

//fg
_DECL cufftDoubleComplex *d_phic    _INIT(NULL);
_DECL cufftDoubleComplex *d_phic_c  _INIT(NULL);

//_DECL cufftDoubleComplex *d_uf  _INIT(NULL);
//_DECL cufftDoubleComplex *d_uf2 _INIT(NULL);

//_DECL cufftDoubleComplex *h_check   _INIT(NULL);
//_DECL double *h_che2 _INIT(NULL);
//_DECL double *h_che22 _INIT(NULL);







//CUFFT PLANS
_DECL cufftHandle plan_z     _INIT(0);
_DECL cufftHandle plan_y     _INIT(0);
_DECL cufftHandle plan_x     _INIT(0);

_DECL cufftHandle plan_x_fwd _INIT(0);
_DECL cufftHandle plan_y_fwd _INIT(0);
_DECL cufftHandle plan_z_fwd _INIT(0);

//CUFFT PLAN MANY
_DECL cufftHandle plan_y_many _INIT(0);

//CUFFT plans fg
_DECL cufftHandle plan_x_psi      _INIT(0);
_DECL cufftHandle plan_x_fwd_psi  _INIT(0);
_DECL cufftHandle plan_y_many_psi _INIT(0);
_DECL cufftHandle plan_z_psi      _INIT(0);


#endif

