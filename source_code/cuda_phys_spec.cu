#include <cuda.h>
#include <stdio.h>
#include <math.h>
#include <cufft.h>
#include <complex.h>

#include "cuda_variables.h"
#include "cuda_phys_spec.h"
#include "cuda_tran.h"

extern "C" void h_ffty_many_fwd(double *in_r, double *in_c, double *out_r, double *out_c, int aliasing)
//-------------------------------------------------------------------------------------
//
//     Subroutine h_ffty_many_fwd from physical to spectral
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{
  dblk = 128;

  ok = ok + cudaMemcpy(ur_d,in_r,sizeof(double)*nsx*npz*ny,cudaMemcpyHostToDevice);
  ok = ok + cudaMemcpy(uc_d,in_c,sizeof(double)*nsx*npz*ny,cudaMemcpyHostToDevice);
  if (ok!=0) printf("Error in copying to device!! ffty fwd\n");

  k_merge_cmp<<<(nsx*ny*npz+dblk-1)/dblk,dblk>>>(d_batch,ur_d,uc_d,nsx*ny*npz);
  if (ok!=0) printf("===============>error in call kernel k_merge_cmp not called!! ffty_fwd\n");
//  ok = ok + cudaMemcpy2D(d_batch,       2 * sizeof(d_batch), in_r, sizeof(in_r), sizeof(in_r), spx*nz*spy, cudaMemcpyHostToDevice);
//  ok = ok + cudaMemcpy2D(&d_batch[0].y, 2 * sizeof(d_batch), in_c, sizeof(in_c), sizeof(in_c), spx*nz*spy, cudaMemcpyHostToDevice);
//  if (ok!=0) printf("Error in copying to device!!\n");


  //BACKWARDS FFT IN Y
  cufftExecZ2Z(plan_y_many, d_batch, d_batch, CUFFT_FORWARD);//this can be done inplace
  ok = cudaGetLastError();
  if (ok!=0) printf("Error in cufftExecZ2Z forward ffty many!!\n");



  //copy back to host
  k_sec_separate<<<(nsx*ny*npz+dblk-1)/dblk,dblk>>>(d_batch,ur_d,uc_d, nsx*ny*npz);
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call kernel k_sec_separate not called!! ffty_fwd\n");

  ok = ok + cudaMemcpy(out_r,ur_d,sizeof(double)*nsx*npz*ny,cudaMemcpyDeviceToHost);
  ok = ok + cudaMemcpy(out_c,uc_d,sizeof(double)*nsx*npz*ny,cudaMemcpyDeviceToHost);
  if (ok!=0) printf("Error in copying to host!!\n");

}



extern "C" void h_fftx_fwd(double *in_r, double *out_r, double *out_c, int aliasing)
//-------------------------------------------------------------------------------------
//
//     Subroutine h_fftx_fwd from physical to spectral
//     WARNING: memory management needs improvement, could consider working on in-place transposition
//     or batched FFTs that need no transposition, despite the current approach grants coalescence for
//     the other kernels
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{
  dblk = 128;

  //copy to device
  ok = ok + cudaMemcpy(ur_d,in_r,sizeof(double)*nx*fpz*fpy,cudaMemcpyHostToDevice);
  if (ok!=0) printf("Error in copying to device!! fftx_fwd\n");


  //execute
  cufftExecD2Z(plan_x_fwd, ur_d, d_batch);
  ok = cudaGetLastError();
  if (ok!=0) printf("Error in forward fft D2Z!!\n");


  //inplace dealiasing
  if(aliasing)
  {
    int al_low = floor(2.0/3.0*(nx/2+1)); //arrays start from zero
	int al_up  = (nx/2+1)-1;//
	k_cmp_alias<<<((al_up-al_low+1)*npy*npz+dblk-1)/dblk,dblk>>>(d_batch, al_low, al_up-al_low+1, (nx/2+1), (nx/2+1)*npy*npz);
	ok = cudaGetLastError();
	if (ok!=0) printf("===============>error in call kernel k_cmp_alias_y not called!! fftx_fwd\n");
  }


  //copy back to host
  k_sec_separate<<<((nx/2+1)*npy*npz+dblk-1)/dblk,dblk>>>(d_batch,ur_d,uc_d,(nx/2+1)*npy*npz);
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call kernel k_sec_separate not called!! fftx_fwd\n");
  ok = ok + cudaMemcpy(out_r,ur_d,sizeof(double)*(nx/2+1)*npz*npy,cudaMemcpyDeviceToHost);
  ok = ok + cudaMemcpy(out_c,uc_d,sizeof(double)*(nx/2+1)*npz*npy,cudaMemcpyDeviceToHost);

//  ok = ok + cudaMemcpy2D(out_r, sizeof(out_r),d_batch,       2 * sizeof(d_batch),sizeof(d_batch), spx*nz*spy, cudaMemcpyDeviceToHost);
//  ok = ok + cudaMemcpy2D(out_c, sizeof(out_c),&d_batch[0].y, 2 * sizeof(d_batch),sizeof(d_batch), spx*nz*spy, cudaMemcpyDeviceToHost);
  if (ok!=0) printf("Error in copying to host!! fftx_fwd\n");


}//end subroutine h_fftx_fwd
//*************************************************************************************
//
//
//
//
//
//
//
//
//
//
//*************************************************************************************
extern "C" void h_ffty_fwd(double *in_r, double *in_c, double *out_r, double *out_c, int aliasing)
//-------------------------------------------------------------------------------------
//
//     Subroutine h_ffty_fwd from physical to spectral
//     WARNING: memory management needs improvement, could consider working on in-place transposition
//     or batched FFTs that need no transposition, despite the current approach grants coalescence for
//     the other kernels
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{
  dblk = 128;
  //a little bit heavy for local GPU, should be good on Teslas!
  int dthr = 8;
  dim3 threadsPerBlock(dthr, dthr, dthr);

  ok = ok + cudaMemcpy(ur_d,in_r, sizeof(double)*spx*npz*ny, cudaMemcpyHostToDevice);
  ok = ok + cudaMemcpy(uc_d,in_c, sizeof(double)*spx*npz*ny, cudaMemcpyHostToDevice);
  if (ok!=0) printf("Error in copying to device!! ffty fwd\n");

  k_merge_cmp<<<(spx*ny*npz+dblk-1)/dblk,dblk>>>(d_batch, ur_d, uc_d, spx*ny*npz);
  if (ok!=0) printf("===============>error in call kernel k_merge_cmp not called!! ffty_fwd\n");


  //BACKWARDS FFT IN Y
  cufftExecZ2Z(plan_y_many, d_batch, d_batch, CUFFT_FORWARD);//this can be done inplace
  ok = cudaGetLastError();
  if (ok!=0) printf("Error in cufftExecZ2Z forward ffty many!!\n");


  //perform inplace aliasing
  if (aliasing == 1)
  {
	//3rd dimension aliasing
	int al_low  = floor(2.0/3.0*(ny/2+1))-1;//first aliased position
	int al_up   = ny-floor(2.0/3.0*(ny/2));//last aliased position
    //grid for working on the x pencils
    dim3 grid_aly((spx+dthr-1)/dthr,(al_up-al_low+dthr-1)/dthr,(npz+dthr-1)/dthr);
	k_alias_3rd<<<grid_aly,threadsPerBlock>>>(d_batch, al_low, al_up-al_low, spx, ny, npz);
    ok = cudaGetLastError();
    if(ok!=0) printf("===============>error in call kernel k_alias_3rd not called!! \n");
  }


  //copy back to host
  k_sec_separate<<<(spx*ny*npz+dblk-1)/dblk,dblk>>>(d_batch, ur_d, uc_d, spx*ny*npz);
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call kernel k_sec_separate not called!! ffty_fwd\n");


  ok = ok + cudaMemcpy(out_r, ur_d, sizeof(double)*spx*npz*ny, cudaMemcpyDeviceToHost);
  ok = ok + cudaMemcpy(out_c, uc_d, sizeof(double)*spx*npz*ny, cudaMemcpyDeviceToHost);
  if (ok!=0) printf("Error in copying to host!! ffty_fwd\n");


}//end subroutine h_ffty_fwd
//*************************************************************************************
//
//
//
//
//
//
//
//
//
//
//*************************************************************************************
extern "C" void h_chebyshev_fwd(double *in_r, double *in_c, double *out_r, double *out_c, int aliasing)
//-------------------------------------------------------------------------------------
//
//     Subroutine h_ffty_fwd from physical to spectral
//     WARNING: memory management needs improvement, could consider working on in-place transposition
//     or batched FFTs that need no transposition, despite the current approach grants coalescence for
//     the other kernels
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{
  dblk = 128;

  ok = ok + cudaMemcpy(ur_d,in_r,sizeof(double)*nsx*nz*npy,cudaMemcpyHostToDevice);
  ok = ok + cudaMemcpy(uc_d,in_c,sizeof(double)*nsx*nz*npy,cudaMemcpyHostToDevice);
  if (ok!=0) printf("Error in copying to device!! dct_fwd\n");


  //grid,block sizes for xzy to zxy transposition
  int TILE_DIM=32;
  int BLOCK_ROWS=8;
  nbx = nsx/TILE_DIM;
  nby = nz/TILE_DIM;
  if (nsx % TILE_DIM != 0) nbx = nbx + 1;
  if (nz  % TILE_DIM != 0) nby = nby + 1;
  dim3 tgrid1(nbx,nby,npy);
  dim3 tBlock(TILE_DIM, BLOCK_ROWS,1);

  //transpose to gain coalescence in z
  k_t102<<<tgrid1,tBlock>>>(d_uopr, ur_d, nsx, nz, npy);//#OUTPUT,INPUT,SIZES
  k_t102<<<tgrid1,tBlock>>>(d_uopc, uc_d, nsx, nz, npy);
  ok = cudaGetLastError();
  if (ok!=0) printf ("===============>error in call kernel k_t102 not called!! \n");


  //make input for cufftD2Z even symmetrical ##format OUTPUT,INPUT,sizes
  k_mirr_bigtime<<<(nz*nsx*npy+dblk-1)/dblk,dblk>>>(ur_d, uc_d, d_uopr, d_uopc, nz, 2*(nz-1), nsx*npy*nz);

//
//  k_mirror<<<spx*spy,32*(nz/32+1),nz*sizeof(double)>>>(ur_d, d_uopr,
//		                                               nz, 2*(nz-1));
//  k_mirror<<<spx*spy,32*(nz/32+1),nz*sizeof(double)>>>(uc_d, d_uopc,
//		                                               nz, 2*(nz-1));
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call kernel k_mirr_bigtime not called - 2!! \n");


  //EXECUTE THE DCT
  cufftExecD2Z(plan_z, ur_d, d_batch);//batched
  cufftExecD2Z(plan_z, uc_d, d_batch_c);
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in cufftExecD2Z - 2!! \n");
  //This transform can also be performed using the Z2D or the Z2Z; each gives different result: WARNING!


  //warning:GPU - the following two kernels should be merged in one kernel
  //Merge the output
  k_sec_copy<<<(spx*spy*nz+dblk-1)/dblk,dblk>>>(d_batch, d_batch_c, nz*spx*spy);//#OUTPUT,INPUT,size -spx*spy,32*(nz/32+1)
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call kernel k_sec_copy not called!! - 2 \n");


  //normalize the result
  int tot_size = spy * spx * nz;
  double norm_fact = 2.0e0 / (double)(nz-1);
  k_norm_cmp<<<(spy*spx*nz+dblk-1)/dblk,dblk>>>(d_batch, norm_fact, tot_size);



  //manipulate last terms
  k_manip_cmp<<<(nz+dblk-1)/dblk,dblk>>>(d_batch, nz, tot_size);

  //inplace aliasing
  if(aliasing)
  {
    int al_low  = floor(2.0*(double)(nz)/3.0);//first aliased position
    int al_up   = nz-1;//last aliased position
    k_cmp_alias<<<((al_up-al_low+1)*spy*spx+dblk-1)/dblk,dblk>>>(d_batch, al_low, al_up-al_low+1, nz, spx*spy*nz);
    ok = cudaGetLastError();
    if (ok!=0) printf("===============>error in kernel k_cmp_alias_y not called - 2!! \n");
  }

    //transpose back from zxy to xzy
    dim3 tgrid2(nby,nbx,spy);
    k_cmp_t102<<<tgrid2,tBlock>>>(d_batch_c, d_batch, nz, spx, spy);//#OUTPUT,INPUT,SIZES
    ok = cudaGetLastError();
    if (ok!=0) printf("===============>error in kernel k_cmp_t102 not called!! \n");


    //copy back to host
    k_sec_separate<<<(spx*spy*nz+dblk-1)/dblk,dblk>>>(d_batch_c,ur_d,uc_d,spx*spy*nz);
    ok = cudaGetLastError();
    if(ok!=0) printf("===============>error in call kernel k_sec_separate not called!! fftx_fwd\n");

    ok = ok + cudaMemcpy(out_r,ur_d,sizeof(double)*spx*nz*spy,cudaMemcpyDeviceToHost);
    ok = ok + cudaMemcpy(out_c,uc_d,sizeof(double)*spx*nz*spy,cudaMemcpyDeviceToHost);


//    ok = ok + cudaMemcpy2D(out_r, sizeof(out_r),d_batch_c,       2 * sizeof(d_batch_c),sizeof(d_batch_c), spx*nz*spy, cudaMemcpyDeviceToHost);
//    ok = ok + cudaMemcpy2D(out_c, sizeof(out_c),&d_batch_c[0].y, 2 * sizeof(d_batch_c),sizeof(d_batch_c), spx*nz*spy, cudaMemcpyDeviceToHost);
    if (ok!=0) printf("Error in copying to host!! dct_fwd\n");

}//end subroutine h_chebyshev_fwd
//*************************************************************************************
//
//
//
//
//
//
//
//
//
//
//*************************************************************************************
//*************************************************************************************
//
//
//
//
//
//
//
//
//
//
//*************************************************************************************
//*************************************************************************************
//
//
//
//
//
//
//
//
//
//
//*************************************************************************************
////EXTRA PIECES-------------------------------------------------------------------------------------------------------
//  cufftDoubleComplex *h_ress=(cufftDoubleComplex*)malloc(sizeof(cufftDoubleComplex) * spx*nz*spy);
//  ok = ok + cudaMemcpy(h_ress,d_uucmp,sizeof(cufftDoubleComplex)*spx*nz*spy,cudaMemcpyDeviceToHost);
//  for (int i=0;i<nz;i++)//<nz*spx*spy;i++)
//  {
//    printf("inside %20.16lf %20.16lf %d\n",h_ress[i].x,h_ress[i].y,i+1);
//  }
//
//cufftDoubleComplex *h_res=(cufftDoubleComplex*)malloc(sizeof(cufftDoubleComplex) * spx*nz*spy);
//ok = ok + cudaMemcpy(h_res,d_uf2,sizeof(cufftDoubleComplex)*spx*nz*spy,cudaMemcpyDeviceToHost);
//for (int i=0;i<nz*spy*spx;i++)
//{
//  printf("uu_spec %20.16le %20.16le %d\n",h_res[i].x,h_res[i].y,i+1);
//}
//double *h_res=(double*)malloc(sizeof(double) * spx*nz*spy);
//double *h_res2=(double*)malloc(sizeof(double) * spx*nz*spy);
//
//ok = ok + cudaMemcpy(h_res,d_urs,sizeof(double)*spx*nz*spy,cudaMemcpyDeviceToHost);
//ok = ok + cudaMemcpy(h_res2,d_ucs,sizeof(double)*spx*nz*spy,cudaMemcpyDeviceToHost);
//
//for (int i=0;i<nz*spy*spx;i++)
//{
//  printf("uu_spec %20.16le %20.16le %d\n",h_res[i],h_res2[i],i+1);
//}