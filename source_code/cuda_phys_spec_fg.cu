#include <cuda.h>
#include <stdio.h>
#include <math.h>
#include <cufft.h>
#include <complex.h>

#include "cuda_variables.h"
#include "cuda_spec_phys_fg.h"
#include "cuda_tran.h"


extern "C" void h_fftxfwd_fg(double *in_r, double *out_r, double *out_c, int aliasing)
//-------------------------------------------------------------------------------------
//
//     Subroutine forward fftx transformation, fg style
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, May 2020
//
//-------------------------------------------------------------------------------------
{
  //npsix,fpzpsi,fpypsi
  dblk = 128;

  //copy to device
  ok = ok + cudaMemcpy(psir_d, in_r, sizeof(double)*npsix*fpzpsi*fpypsi, cudaMemcpyHostToDevice);
  if (ok!=0) printf("Error in copying to device!! fg\n");


  //execute
  cufftExecD2Z(plan_x_fwd_psi, psir_d, d_phic);
  ok = cudaGetLastError();
  if (ok!=0) printf("Error in forward fft D2Z!! fg\n");


  //inplace dealiasing
  if(aliasing)
  {
	//ERROR HERE??
    int al_low = floor(2.0/3.0*(npsix/2+1)); //arrays start from zero
	int al_up  = (npsix/2+1)-1;//
	k_cmp_alias<<<((al_up-al_low+1)*fpypsi*fpzpsi+dblk-1)/dblk,dblk>>>(d_phic, al_low, al_up-al_low+1, (npsix/2+1), (npsix/2+1)*fpypsi*fpzpsi);
	ok = cudaGetLastError();
	if (ok!=0) printf("===============>error in call kernel k_cmp_alias not called!! fg\n");
  }


  //copy back to host
  k_sec_separate<<<((npsix/2+1)*fpypsi*fpzpsi+dblk-1)/dblk,dblk>>>(d_phic, psir_d, psic_d, (npsix/2+1)*fpypsi*fpzpsi);
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call kernel k_sec_separate not called!! fg\n");
  ok = ok + cudaMemcpy(out_r, psir_d, sizeof(double)*(npsix/2+1)*fpzpsi*fpypsi, cudaMemcpyDeviceToHost);
  ok = ok + cudaMemcpy(out_c, psic_d, sizeof(double)*(npsix/2+1)*fpzpsi*fpypsi, cudaMemcpyDeviceToHost);

  if (ok!=0) printf("Error in copying to host!! fg\n");


}//end subroutine h_chebback_fg
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
extern "C" void h_chebfwd_fg(double *in_r, double *in_c, double *out_r, double *out_c, int aliasing)
//-------------------------------------------------------------------------------------
//
//     Subroutine forward chebyshev transformation, fg style
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, May 2020
//
//-------------------------------------------------------------------------------------
{
  dblk = 128;

  ok = ok + cudaMemcpy(psir_d, in_r, sizeof(double)*spxpsi*npsiz*spypsi, cudaMemcpyHostToDevice);
  ok = ok + cudaMemcpy(psic_d, in_c, sizeof(double)*spxpsi*npsiz*spypsi, cudaMemcpyHostToDevice);
  if (ok!=0) printf("Error in copying to device!! fg\n");


  //grid,block sizes for xzy to zxy transposition
  int TILE_DIM=32;
  int BLOCK_ROWS=8;
  nbx = spxpsi/TILE_DIM;
  nby = npsiz/TILE_DIM;
  if (spxpsi % TILE_DIM != 0) nbx = nbx + 1;
  if (npsiz  % TILE_DIM != 0) nby = nby + 1;
  dim3 tgrid1(nbx,nby,spypsi);
  dim3 tBlock(TILE_DIM, BLOCK_ROWS,1);

  //transpose to gain coalescence in z
  k_t102<<<tgrid1,tBlock>>>(d_psir, psir_d, spxpsi, npsiz, spypsi);//#OUTPUT,INPUT,SIZES
  k_t102<<<tgrid1,tBlock>>>(d_psic, psic_d, spxpsi, npsiz, spypsi);
  ok = cudaGetLastError();
  if (ok!=0) printf ("===============>error in call kernel k_t102 not called!! fg\n");


  //make input for cufftD2Z even symmetrical ##format OUTPUT,INPUT,sizes
  k_mirr_bigtime<<<(npsiz*spxpsi*spypsi+dblk-1)/dblk,dblk>>>(psir_d, psic_d, d_psir, d_psic, npsiz, 2*(npsiz-1), spxpsi*spypsi*npsiz);
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call kernel k_mirr_bigtime not called!! fg\n");


  //EXECUTE THE DCT
  cufftExecD2Z(plan_z_psi, psir_d, d_phic);//batched
  cufftExecD2Z(plan_z_psi, psic_d, d_phic_c);
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in cufftExecD2Z !! fg\n");


  //Merge the output
  k_sec_copy<<<(spxpsi*spypsi*npsiz+dblk-1)/dblk,dblk>>>(d_phic, d_phic_c, npsiz*spxpsi*spypsi);//#OUTPUT,INPUT,size -spx*spy,32*(nz/32+1)
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call kernel k_sec_copy not called!! fg\n");


  //normalize the result
  int tot_size = spypsi * spxpsi * npsiz;
  double norm_fact = 2.0e0 / (double)(npsiz-1);
  k_norm_cmp<<<(spypsi*spxpsi*npsiz+dblk-1)/dblk,dblk>>>(d_phic, norm_fact, tot_size);


  //manipulate last terms
  k_manip_cmp<<<(npsiz+dblk-1)/dblk,dblk>>>(d_phic, npsiz, tot_size);


  //inplace aliasing
  if(aliasing)
  {
    int al_low  = floor(2.0*(double)(npsiz)/3.0);//first aliased position
	int al_up   = npsiz-1;//last aliased position
    k_cmp_alias<<<((al_up-al_low+1)*spypsi*spxpsi+dblk-1)/dblk,dblk>>>(d_phic, al_low, al_up-al_low+1, npsiz, spxpsi*spypsi*npsiz);
	ok = cudaGetLastError();
    if (ok!=0) printf("===============>error in kernel k_cmp_alias_y not called!! fg \n");
  }

  //transpose back from zxy to xzy
  dim3 tgrid2(nby,nbx,spypsi);
  k_cmp_t102<<<tgrid2,tBlock>>>(d_phic_c, d_phic, npsiz, spxpsi, spypsi);//#OUTPUT,INPUT,SIZES
  ok = cudaGetLastError();
  if (ok!=0) printf("===============>error in kernel k_cmp_t102 not called!!  fg\n");


  //copy back to host
  k_sec_separate<<<(spxpsi*spypsi*npsiz+dblk-1)/dblk,dblk>>>(d_phic_c, psir_d, psic_d, spxpsi*spypsi*npsiz);
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call kernel k_sec_separate not called!! fg\n");

  ok = ok + cudaMemcpy(out_r, psir_d, sizeof(double)*spxpsi*npsiz*spypsi, cudaMemcpyDeviceToHost);
  ok = ok + cudaMemcpy(out_c, psic_d, sizeof(double)*spxpsi*npsiz*spypsi, cudaMemcpyDeviceToHost);

  if (ok!=0) printf("Error in copying to host!! fg\n");

}//end subroutine h_fftxback_fg
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
extern "C" void h_fftymanyfwd_fg(double *in_r, double *in_c, double *out_r, double *out_c, int aliasing)
//-------------------------------------------------------------------------------------
//
//     Subroutine forward FFTy transformation, fg style
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, May 2020
//
//-------------------------------------------------------------------------------------
{
  dblk = 128;

  ok = ok + cudaMemcpy(psir_d, in_r,sizeof(double)*spxpsi*fpzpsi*npsiy, cudaMemcpyHostToDevice);
  ok = ok + cudaMemcpy(psic_d, in_c,sizeof(double)*spxpsi*fpzpsi*npsiy, cudaMemcpyHostToDevice);
  if (ok!=0) printf("Error in copying to device!! fg\n");

  k_merge_cmp<<<(spxpsi*npsiy*fpzpsi+dblk-1)/dblk,dblk>>>(d_phic, psir_d, psic_d, spxpsi*npsiy*fpzpsi);
  if (ok!=0) printf("===============>error in call kernel k_merge_cmp not called!! fg\n");


  //BACKWARDS FFT IN Y
  cufftExecZ2Z(plan_y_many_psi, d_phic, d_phic, CUFFT_FORWARD);//this can be done inplace
  ok = cudaGetLastError();
  if (ok!=0) printf("Error in cufftExecZ2Z forward ffty many!! fg\n");


  //copy back to host
  k_sec_separate<<<(spxpsi*npsiy*fpzpsi+dblk-1)/dblk,dblk>>>(d_phic, psir_d, psic_d, spxpsi*npsiy*fpzpsi);
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call kernel k_sec_separate not called!! psi\n");

  ok = ok + cudaMemcpy(out_r, psir_d, sizeof(double)*spxpsi*fpzpsi*npsiy, cudaMemcpyDeviceToHost);
  ok = ok + cudaMemcpy(out_c, psic_d, sizeof(double)*spxpsi*fpzpsi*npsiy, cudaMemcpyDeviceToHost);
  if (ok!=0) printf("Error in copying to host!!\n");

}//end subroutine h_fftymanyfwd_fg
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
extern "C" void h_ffty_fwd_fg(double *in_r, double *in_c, double *out_r, double *out_c, int aliasing)
{
  dblk = 128;
  int dthr = 8;
  dim3 threadsPerBlock(dthr, dthr, dthr);

  ok = ok + cudaMemcpy(psir_d, in_r,sizeof(double)*spxpsi*fpzpsi*npsiy, cudaMemcpyHostToDevice);
  ok = ok + cudaMemcpy(psic_d, in_c,sizeof(double)*spxpsi*fpzpsi*npsiy, cudaMemcpyHostToDevice);
  if (ok!=0) printf("Error in copying to device!! fg\n");

  k_merge_cmp<<<(spxpsi*npsiy*fpzpsi+dblk-1)/dblk,dblk>>>(d_phic, psir_d, psic_d, spxpsi*npsiy*fpzpsi);
  if (ok!=0) printf("===============>error in call kernel k_merge_cmp not called!! fg\n");



  //BACKWARDS FFT IN Y
  cufftExecZ2Z(plan_y_many_psi, d_phic, d_phic, CUFFT_FORWARD);//this can be done inplace
  ok = cudaGetLastError();
  if (ok!=0) printf("Error in cufftExecZ2Z forward ffty many!! fg\n");


  //perform inplace aliasing
  if (aliasing == 1)
  {
	//3rd dimension aliasing
	int al_low  = floor(2.0/3.0*(npsiy/2+1))-1;//first aliased position
	int al_up   = npsiy-floor(2.0/3.0*(npsiy/2));//last aliased position
    //grid for working on the x pencils
    dim3 grid_aly((spxpsi+dthr-1)/dthr,(al_up-al_low+dthr-1)/dthr,(fpzpsi+dthr-1)/dthr);
	k_alias_3rd<<<grid_aly,threadsPerBlock>>>(d_phic, al_low, al_up-al_low, spxpsi, npsiy, fpzpsi);
    ok = cudaGetLastError();
    if(ok!=0) printf("===============>error in call kernel k_alias_3rd not called!! fg\n");
  }


  //copy back to host
  k_sec_separate<<<(spxpsi*npsiy*fpzpsi+dblk-1)/dblk,dblk>>>(d_phic, psir_d, psic_d, spxpsi*npsiy*fpzpsi);
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call kernel k_sec_separate not called!! psi\n");

  ok = ok + cudaMemcpy(out_r, psir_d, sizeof(double)*spxpsi*fpzpsi*npsiy, cudaMemcpyDeviceToHost);
  ok = ok + cudaMemcpy(out_c, psic_d, sizeof(double)*spxpsi*fpzpsi*npsiy, cudaMemcpyDeviceToHost);
  if (ok!=0) printf("Error in copying to host!!\n");


}//end subroutine h_ffty_fwd_fg
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
