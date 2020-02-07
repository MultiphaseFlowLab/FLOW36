#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cufft.h>
#include <complex.h>

//#define VAR_DECLS
#include "cuda_variables.h"
#include "cuda_surf.h"
#include "cuda_spec_phys.h"
#include "cuda_tran.h"
#include "cuda_phys_spec.h"

#define MODULO(X,Y) (((X)%(Y))>=0) ? ((X)%(Y)) : ((X)%(Y) + (Y))


//__global__ void k_mirror_c(int N, cufftDoubleReal *in, cufftDoubleComplex *out)
//{
//	  int index = blockIdx.x*blockDim.x + threadIdx.x;
//	  //per avere performance va fatto con passaggio in shared memory per la scrittura!!!
//	  if(index < N && index > 0)
//	  {
//	    out[index].x=in[index];
//	    out[index].y=0.0e0;
//	    out[ 2*(N -1) - index ].x =  in[index];
//	    out[ 2*(N -1) - index ].y = 0.0e0;
//	  }
//	  out[0].x=in[0];
//	  out[0].y=0.0e0;
//}

//__global__ void k_mirror(int N, cufftDoubleReal *in, cufftDoubleReal *out)
//{
//	  int index = blockIdx.x*blockDim.x + threadIdx.x;
//	  //per avere performance va fatto con passaggio in shared memory per la scrittura!!!
//	  if(index < N && index > 0)
//	  {
//	    out[index]=in[index];
//	    out[ 2*(N -1) - index ] =  in[index];
//	  }
//	  out[0]=in[0];
//}


//__global__ void saxpy(cufftDoubleComplex *in, double *a, double *b, int dim)
//{
//  int index = blockIdx.x*blockDim.x + threadIdx.x;
//  if (index < dim)
//  {
//    a[index] = in[index].x;
//    b[index] = in[index].y;
////  printf ("%d %le \n",i,a[i]);
//  }
//}

extern "C" void h_mainloop_()
{

  //integers
  extern int spx_f;
  extern int spy_f;
  extern int nz_f;
  extern int nx_f;

  //doubles
  extern double xl_f;
  extern double yl_f;
  extern double dt_f;


  //double arrays
  extern double *ur_f;
  extern double *uc_f;
  extern double *uout_f;

  extern double *vr_f;
  extern double *vc_f;
  extern double *vout_f;

  extern double *wr_f;
  extern double *wc_f;
  extern double *wout_f;

  extern double *z_f;

  //integer arrays
  extern int *fstart_f;



  //local integers
  int i,j,k,ind;
  int ok = 0;

  //pass sizes
  spx = spx_f;
  spy = spy_f;
  nz  = nz_f;
  nx  = nx_f;
  xl  = xl_f;
  yl  = yl_f;
  dt  = dt_f;

  //pass control parameters
  aliasing = 1;
  dblk = 128;

  //pass integer arrays
  z      = z_f;
  fstart = fstart_f;

  //calculate dz
  dz_h = (double *)malloc(sizeof(double) * nz);
  dz_h[0]=z[0]-z[1];
  dz_h[nz-1]=z[nz-2]-z[nz-1];
  for (k=1;k<nz-1;k++)
  {
	  dz_h[k]=(z[k-1]-z[k+1])/2.0e0;
  }

  //calculate largest size for double arrays
  int d_dim = nx * nz * spy;
  if (spx*spy*2*(nz-1) > nx * nz * spy)
  {
    d_dim = spx * spy * 2 * (nz-1);
  }


  //pass realpr arrays
  ur = ur_f;
  uc = uc_f;

  vr = vr_f;
  vc = vc_f;

  wr = wr_f;
  wc = wc_f;

  //allocate on GPU
  //integer arrays
  ok = ok + cudaMalloc((void **)&dz_d,      sizeof(double)*nz);
  ok = ok + cudaMalloc((void **)&fstart_d, sizeof(int)*nz);


  //double arrays
  ok = ok + cudaMalloc((void **)&ur_d,    sizeof(double)*d_dim);
  ok = ok + cudaMalloc((void **)&uc_d,    sizeof(double)*d_dim);
  ok = ok + cudaMalloc((void **)&d_uout,  sizeof(double)*nx*spy*nz);

  ok = ok + cudaMalloc((void **)&d_uopr,  sizeof(double)*d_dim);
  ok = ok + cudaMalloc((void **)&d_uopc,  sizeof(double)*d_dim);
  ok = ok + cudaMalloc((void **)&d_uext,  sizeof(double)*d_dim);

  ok = ok + cudaMalloc((void **)&vr_d,    sizeof(double)*spx*nz*spy);
  ok = ok + cudaMalloc((void **)&vc_d,    sizeof(double)*spx*nz*spy);
  ok = ok + cudaMalloc((void **)&d_vout,  sizeof(double)*nx*spy*nz);

  ok = ok + cudaMalloc((void **)&wr_d,    sizeof(double)*spx*nz*spy);
  ok = ok + cudaMalloc((void **)&wc_d,    sizeof(double)*spx*nz*spy);
  ok = ok + cudaMalloc((void **)&d_wout,  sizeof(double)*nx*spy*nz);
  if (ok!=0) printf ("Allocation of double arrays failed\n");

  //allocate output arrays for Chebyshev inverse DCT
  ok = ok + cudaMalloc((void **)&d_batch,  sizeof(cufftDoubleComplex)*nz*spx*spy);
  ok = ok + cudaMalloc((void **)&d_batch_c,sizeof(cufftDoubleComplex)*nz*spx*spy);
  if (ok!=0) printf("Output arrays for Chebyshev not allocated correctly!!\n");

  //mixed products
  ok = ok + cudaMalloc((void **)&d_uu,   sizeof(double)*nx*nz*spy);
  ok = ok + cudaMalloc((void **)&d_uuc,  sizeof(cufftDoubleComplex)*nz*spx*spy);
  ok = ok + cudaMalloc((void **)&d_uuc2,  sizeof(cufftDoubleComplex)*nz*spx*spy);

  if (ok!=0) printf("Failed allocation of mixed product arrays!\n");


  //checks
  h_check=(cufftDoubleComplex *)malloc(sizeof(cufftDoubleComplex) * spx*spy*nz);
  h_che2 =(double *)malloc(sizeof(double) * nx*spy*2*(nz-1));
  h_che22 =(double *)malloc(sizeof(double) * nx*spy*2*(nz-1));


  //output
//  h_uout=(double *)malloc(sizeof(double) * nx*spy*nz);
//  h_vout=(double *)malloc(sizeof(double) * nx*spy*nz);
//  h_wout=(double *)malloc(sizeof(double) * nx*spy*nz);
  h_uur = (double *)malloc(sizeof(double) * nx*spy*nz);
  h_uuc = (double *)malloc(sizeof(double) * nx*spy*nz);
  if (ok!=0) printf("Allocation of output array failed!!\n");

  //create plans for cufft transformations
  cufftPlan1d(&plan_z, (nz-1)*2, CUFFT_D2Z, spx*spy);
  cufftPlan1d(&plan_y,      spy, CUFFT_Z2Z, spx*nz);
  cufftPlan1d(&plan_x,       nx, CUFFT_Z2D, spy*nz);
  ok = cudaGetLastError();
  if (ok!=0) printf("Error in creating the backward plans!!\n");

  cufftPlan1d(&plan_x,  nx, CUFFT_D2Z, spy*nz);
  cufftPlan1d(&plan_y, spy, CUFFT_Z2Z, spx*nz);
  ok = cudaGetLastError();
  if (ok!=0) printf("Error in creating the forward plans !!\n");



  //copy to the GPU
  //doubles
  ok = ok + cudaMemcpy(ur_d,ur,sizeof(double)*spx*nz*spy,cudaMemcpyHostToDevice);
  ok = ok + cudaMemcpy(uc_d,uc,sizeof(double)*spx*nz*spy,cudaMemcpyHostToDevice);

  ok = ok + cudaMemcpy(vr_d,vr,sizeof(double)*spx*nz*spy,cudaMemcpyHostToDevice);
  ok = ok + cudaMemcpy(vc_d,vc,sizeof(double)*spx*nz*spy,cudaMemcpyHostToDevice);

  ok = ok + cudaMemcpy(wr_d,wr,sizeof(double)*spx*nz*spy,cudaMemcpyHostToDevice);
  ok = ok + cudaMemcpy(wc_d,wc,sizeof(double)*spx*nz*spy,cudaMemcpyHostToDevice);

  //integers z
  ok = ok + cudaMemcpy(dz_d,dz_h,sizeof(double)*nz,cudaMemcpyHostToDevice);

  if(ok!=0) printf("Failed to copy the data to the GPU!!\n");

  //show memory usage
  size_t free_byte ;
  size_t total_byte ;
  cudaError_t cuda_status = cudaMemGetInfo( &free_byte, &total_byte ) ;

  if (cudaSuccess != cuda_status)
  {
    printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status) );
  }
  else
  {
    double free_db  = (double)free_byte ;
    double total_db = (double)total_byte ;
    double used_db  = total_db - free_db ;
    printf("\n GPU memory usage: used %5.2f GB, total (free+used) %5.2f GB, percentage %5.2f\n",used_db/(1073741824),total_db/(1073741824),used_db/total_db*100);
  }












//INPUT is ordered as xzy
  cuda_spectral_to_phys(ur_d,uc_d,d_uout);
  cuda_spectral_to_phys(vr_d,vc_d,d_vout);
  cuda_spectral_to_phys(wr_d,wc_d,d_wout);
  //WARNING: OUTPUT is ordered as xyz

//INTERMEZZO
  //  ok = ok + cudaMemcpy(h_che2,d_uout,sizeof(double)*nx*nz*spy,cudaMemcpyDeviceToHost);//store d_uout in physical space


  //WARNING: INPUT is ordered as xyz
//  cuda_phys_to_spectral(d_uout,d_uuc);
  //WARNING: OUTPUT is ordered as zxy

  //zxy to xzy
  //012    102


////  ok = ok + cudaMemcpy(h_check,d_uuc,sizeof(cufftDoubleComplex)*spx*nz*spy,cudaMemcpyDeviceToHost);
//
//  //TEST
//  int TILE_DIM =32;
  //  int BLOCK_ROWS = 8;
//  nbx = nz/TILE_DIM;
//  nby = spx/TILE_DIM;
//  if (nz   % TILE_DIM != 0) nbx = nbx + 1;
//  if (spx  % TILE_DIM != 0) nby = nby + 1;
//  dim3 tt(nbx,nby,spy);
//  dim3 tB(TILE_DIM, BLOCK_ROWS,1);
//  k_cmp_t102<<<tt,tB>>>(d_uuc2,d_uuc,nz,spx,spy);
//  ok = cudaGetLastError();
//  if(ok!=0)printf("ERROR HERE\n");
//  //END
//
//  saxpy<<<(spx*spy*nz+128-1)/128,128>>>(d_uuc2,ur_d,uc_d,spx*spy*nz);
//  ok = cudaGetLastError();
////	cudaMemcpy2D (resultReal, 1 * sizeof(resultReal),complex_vec_d, 2 * sizeof(complex_vec_d),sizeof(complex_vec_d), spx*nz*spy, cudaMemcpyDeviceToHost);
//  if (ok!=0) printf("ERROR HERE 2\n");
//
//
////INPUT is ordered as xzy
//  cuda_spectral_to_phys(ur_d,uc_d,d_uout);
////WARNING: OUTPUT is ordered as zxy
//
//
//  ok = ok + cudaMemcpy(h_che22,d_uout,sizeof(double)*nx*nz*spy,cudaMemcpyDeviceToHost);//store d_uout in physical space
//
//
//
//  double errore;
//  for (i=0;i<nz*nx*spy;i++)
//  {
//	  printf("original %20.16lf \n",h_che22[i]);
////	  printf("original %20.16lf \n",h_che2[i]);
//
//
////	  printf("original %20.16lf %20.16lf \n",ur[i*spx],uc[i*spx]);
////	  printf("original %20.16lf %20.16lf \n",h_check[i].x,h_check[i].y);
////	  printf("hola %20.16lf %20.16lf - result %20.16lf %20.16lf\n",ur[i*spx],uc[i*spx],h_check[i].x,h_check[i].y);
//	  errore = (h_che2[i*spx]-h_che22[i])/h_che2[i*spx]*100.0e0;
////      printf("Error on real component %5.4le  !!\n", errore);
//  }
//INTERMEZZO

  k_mixed_prod<<<(spy*nx*nz+dblk-1)/dblk,dblk>>>(d_uout, d_uout, d_uu, nx*spy*nz);
//  k_mixed_prod<<<(spy*nx*nz+dblk-1)/dblk,dblk>>>(d_uout, d_vout, d_uv, nx*spy*nz);
//  k_mixed_prod<<<(spy*nx*nz+dblk-1)/dblk,dblk>>>(d_uout, d_wout, d_uw, nx*spy*nz);
//  k_mixed_prod<<<(spy*nx*nz+dblk-1)/dblk,dblk>>>(d_vout, d_vout, d_vv, nx*spy*nz);
//  k_mixed_prod<<<(spy*nx*nz+dblk-1)/dblk,dblk>>>(d_vout, d_wout, d_vw, nx*spy*nz);
//  k_mixed_prod<<<(spy*nx*nz+dblk-1)/dblk,dblk>>>(d_wout, d_wout, d_ww, nx*spy*nz);

  ok = cudaGetLastError();
  if (ok!=0) printf("===============>error in call kernel k_mixed_prod for uu not called!! \n");

  //calling nz blocks with nx*spy threads; each block receives in shared memory the proper value for dz
  //this kernel overwrites the u velocity therefore attention should be paid

  dx = xl / (double)(nx-1);
  dy = yl / (double)(spy-1);

  ok = ok + cudaMalloc ((void **) &res,	         sizeof(double)*1);
  k_courant <<<nz,nx*spy,nx*spy*sizeof(double)>>>(d_uout, d_vout, d_wout, dz_d, res,
		                                          dx, dy, dt,
		                                          nx*nz*spy, fstart[1], nx, spy); //how to pass the constant dz to shared memory?? pass the whole dz vector
  ok = cudaGetLastError();
  if (ok!=0) printf("===============>error in call kernel k_courant not called!! \n");

  ok = ok + cudaMemcpy(&cou_val,res,sizeof(double)*1,cudaMemcpyDeviceToHost);
  printf("Courant check: %20.16f\n",cou_val);



  //WARNING: INPUT is ordered as xyz
  cuda_phys_to_spectral(d_uu,d_uuc);
//  cuda_phys_to_spectral(d_uv,d_uvc);
//  cuda_phys_to_spectral(d_uw,d_uwc);
//  cuda_phys_to_spectral(d_vv,d_vvc);
//  cuda_phys_to_spectral(d_vw,d_vwc);
//  cuda_phys_to_spectral(d_ww,d_wwc);

  //WARNING: OUTPUT is ordered as zxy


  //copy back the mixed product in spectral space
//  ok = ok + cudaMemcpy2D (&h_uur[0], sizeof(h_uur[0]),&d_uuc[0].x, 1 * sizeof(d_uuc[0]),sizeof(d_uuc[0]), spx*nz*spy, cudaMemcpyDeviceToHost);
//  ok = ok + cudaMemcpy2D (&h_uuc[0], sizeof(h_uuc[0]),&d_uuc[0].y, 1 * sizeof(d_uuc[0]),sizeof(d_uuc[0]), spx*nz*spy, cudaMemcpyDeviceToHost);




//  ok = ok + cudaMemcpy(h_vout,d_vout,sizeof(double)*nx*nz*spy,cudaMemcpyDeviceToHost);
//  ok = ok + cudaMemcpy(h_wout,d_wout,sizeof(double)*nx*nz*spy,cudaMemcpyDeviceToHost);

  //  uout_f = h_uout;
//  uout_f = h_uout;
//  vout_f = h_vout;
//  wout_f = h_wout;



  //---------------------------------------------free memory--------------------------------------------------
//  free(ur);
//  free(uc);
  //free cuda arrays
//  free(h_uout);
//  free(h_res);

  free(h_check);
  free(h_che2);
  free(dz_h);

  ok = 0;

  ok = ok + cudaFree(d_batch);
  ok = ok + cudaFree(d_batch_c);


  ok = ok + cudaFree(ur_d);
  ok = ok + cudaFree(uc_d);
  ok = ok + cudaFree(d_uout);

  ok = ok + cudaFree(d_uopr);
  ok = ok + cudaFree(d_uopc);
  ok = ok + cudaFree(d_uext);

  //  ok = ok + cudaFree(vr_d);
//  ok = ok + cudaFree(vc_d);
//  ok = ok + cudaFree(d_vout);

//  ok = ok + cudaFree(wr_d);
//  ok = ok + cudaFree(wc_d);
//  ok = ok + cudaFree(d_wout);

  ok = ok + cudaFree(d_uu);

  ok = ok + cudaFree(dz_d);
  ok = ok + cudaFree(fstart_d);
  if (ok!=0) printf("Failure in freeing the arrays!!\n");

  ok = ok + cufftDestroy(plan_x);
  ok = ok + cufftDestroy(plan_y);
  ok = ok + cufftDestroy(plan_z);
//  ok = ok + cufftDestroy(plan_z_fwd);

  if (ok!=0) printf("Error in planes destruction!!\n");

  cudaDeviceReset();

}


////extra pieces
//cufftDoubleComplex *h_res=(cufftDoubleComplex*)malloc(sizeof(cufftDoubleComplex) * spx*nz*spy);
//ok = ok + cudaMemcpy (h_res,u1_d, sizeof(cufftDoubleComplex) * spx*nz*spy, cudaMemcpyDeviceToHost);
//for (i=0;i<nz;i++)
//{
//	  !printf("cuda %18.12le abs_ind %d\n",h_res[i].x,i);
//}
//
//
//
//  double *h_uout=(double*)malloc(sizeof(double) * nx*spy*nz);
//
//  ok = ok + cudaMemcpy (h_uout,d_uout,nz*nx*spy*sizeof(double),cudaMemcpyDeviceToHost); //uout order is xyz
//  k=15;
//  j=15;
////  printf("output value absolute index x y z\n");
//  for (k=0;k<nz;k++)//z last index
//    for (j=0;j<spy;j++)//y second index
//      for (i=0;i<nx;i++)//x first index
//      {
//        ind = i + nx * (j + k * spy);
//        printf("output %20.16le %d\n",h_uout[ind],ind+1);
//      }

//ok = ok + cudaMemcpy2D (&ur_d[0].x, sizeof(ur_d[0]),ur, 1 * sizeof(ur[0]),sizeof(ur[0]), spx*nz*spy, cudaMemcpyHostToDevice);
//ok = ok + cudaMemcpy2D (&uc_d[0].x, sizeof(uc_d[0]),uc, 1 * sizeof(uc[0]),sizeof(uc[0]), spx*nz*spy, cudaMemcpyHostToDevice);



