#include <cuda.h>
#include <stdio.h>
#include <math.h>
#include <cufft.h>
#include <complex.h>

#include "cuda_variables.h"
#include "cuda_spec_phys_fg.h"
#include "cuda_tran.h"



extern "C" void h_chebback_fg(double *in_r, double *in_c, double *out_r, double *out_c, int aliasing)
//-------------------------------------------------------------------------------------
//
//     Subroutine inverse chebyshev transformation, fg style
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, May 2020
//
//-------------------------------------------------------------------------------------
{

  dblk = 128;

  ok = ok + cudaMemcpy(psir_d, in_r,sizeof(double)*spxpsi*npsiz*spypsi, cudaMemcpyHostToDevice);
  ok = ok + cudaMemcpy(psic_d, in_c,sizeof(double)*spxpsi*npsiz*spypsi, cudaMemcpyHostToDevice);
  if (ok!=0) printf("Error in copying to device fg!!\n");


  //grid,block sizes for xzy to zxy transposition
  int TILE_DIM=32;
  int BLOCK_ROWS=8;
  nbx = spxpsi/TILE_DIM;
  nby = npsiz/TILE_DIM;
  if (spxpsi % TILE_DIM != 0) nbx = nbx + 1;
  if (npsiz  % TILE_DIM != 0) nby = nby + 1;
  dim3 tgrid1(nbx,nby,spypsi);
  dim3 tBlock(TILE_DIM, BLOCK_ROWS,1);

  //transpose from xzy to zxy to gain coalescence in z
  k_t102<<<tgrid1,tBlock>>>(d_psir, psir_d, spxpsi, npsiz, spypsi);//#OUTPUT,INPUT,SIZES
  k_t102<<<tgrid1,tBlock>>>(d_psic, psic_d, spxpsi, npsiz, spypsi);
  ok = cudaGetLastError();
  if (ok!=0) printf ("===============>error in call kernel k_t102 not called fg!! \n");



  //perform inplace aliasing
  if (aliasing == 1)
  {
    int ali = floor(2.0*double(npsiz)/3.0); //arrays start from zero
    k_alias_1st<<<((npsiz-ali)*spxpsi*spypsi+dblk-1)/dblk,dblk>>>(d_psir, d_psic, ali, npsiz, spxpsi*spypsi*npsiz);
    ok = cudaGetLastError();
    if (ok!=0) printf("===============>error in call kernel k_alias_1st not called fg!! \n");
  }


  //manipulate the input arrays
  k_manip<<<(npsiz+dblk-1)/dblk,dblk>>>(d_psir, npsiz, npsiz*spypsi*spxpsi);
  k_manip<<<(npsiz+dblk-1)/dblk,dblk>>>(d_psic, npsiz, npsiz*spypsi*spxpsi);
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call kernel k_manip not called!! fg \n");



  //make input for cufftD2Z even symmetrical ##format OUTPUT,INPUT,sizes
  k_mirr_bigtime<<<(npsiz*spxpsi*spypsi+dblk-1)/dblk,dblk>>>(psir_d, psic_d, d_psir,  d_psic, npsiz, 2*(npsiz-1), spxpsi*spypsi*npsiz);
  ok = cudaGetLastError();
  if (ok!=0) printf("===============>error in call kernel k_mirr_bigtime not called!! fg \n");



  //execute the inverse Chebyshev DCT
  cufftExecD2Z(plan_z_psi, psir_d, d_phic);
  cufftExecD2Z(plan_z_psi, psic_d, d_phic_c);
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call cufftD2Z not called!! \n");


  //merge the two outputs of the Chebyshev DFT and normalize
  double norm_fac = 0.5e0;
  k_sep_cmp<<<(spxpsi*spypsi*npsiz+dblk-1)/dblk,dblk>>>(d_phic,   psir_d, norm_fac, npsiz*spxpsi*spypsi);
  k_sep_cmp<<<(spxpsi*spypsi*npsiz+dblk-1)/dblk,dblk>>>(d_phic_c, psic_d, norm_fac, npsiz*spxpsi*spypsi);
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call kernel k_sec_copy not called!! fg \n");


  //transpose back from zxy to xzy
  dim3 tgrid2(nby,nbx,spypsi);
  k_t102<<<tgrid2,tBlock>>>(d_psir, psir_d, npsiz, spxpsi, spypsi);
  k_t102<<<tgrid2,tBlock>>>(d_psic, psic_d, npsiz, spxpsi, spypsi);
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call kernel k_transpose_back not called!! fg\n");

  ok = ok + cudaMemcpy(out_r, d_psir, sizeof(double)*spxpsi*npsiz*spypsi, cudaMemcpyDeviceToHost);
  ok = ok + cudaMemcpy(out_c, d_psic, sizeof(double)*spxpsi*npsiz*spypsi, cudaMemcpyDeviceToHost);
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
extern "C" void h_fftxback_fg(double *in_r, double *in_c, double *out_r, int aliasing)
//-------------------------------------------------------------------------------------
//
//     Subroutine inverse fftx transformation, fg style
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, May 2020
//
//-------------------------------------------------------------------------------------
{
	//npsix,fpzpsi,fpypsi


  dblk = 128;

  //copy to device
  ok = ok + cudaMemcpy(psir_d, in_r, sizeof(double)*(npsix/2+1)*fpypsi*fpzpsi, cudaMemcpyHostToDevice);
  ok = ok + cudaMemcpy(psic_d, in_c, sizeof(double)*(npsix/2+1)*fpypsi*fpzpsi, cudaMemcpyHostToDevice);
  if (ok!=0) printf("Error in copying to device!! fftx back\n");
  k_merge_cmp<<<((npsix/2+1)*fpypsi*fpzpsi+dblk-1)/dblk,dblk>>>(d_phic, psir_d, psic_d, (npsix/2+1)*fpypsi*fpzpsi);
  if (ok!=0) printf("===============>error in call kernel k_merge_cmp not called!! fg\n");


  if (aliasing == 1)
  {
    int al_low = floor(2.0/3.0*(npsix/2+1)); //arrays start from zero
    k_alias_1st_cmp<<<(((npsix/2+1)-al_low)*fpypsi*fpzpsi+dblk-1)/dblk,dblk>>>(d_phic, al_low, (npsix/2+1), (npsix/2+1)*fpypsi*fpzpsi);
	ok = cudaGetLastError();
    if (ok!=0) printf("===============>error in second call kernel k_alias_1st not called!! fg\n");
  }



  //BACKWARD FFT to physical space
  cufftExecZ2D(plan_x_psi, d_phic, psir_d);


  //Normalize
  double norm_fact = 1.0e0/(double)npsix;
  int tot_size  = npsix*fpypsi*fpzpsi;
  k_norm_real<<<(fpypsi*npsix*fpzpsi+dblk-1)/dblk,dblk>>>(psir_d, norm_fact, tot_size);


  //copy back to host
  ok = ok + cudaMemcpy(out_r, psir_d, sizeof(double)*npsix*fpzpsi*fpypsi, cudaMemcpyDeviceToHost);
  if (ok!=0) printf("Error in copying back to host!\n");



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
extern "C" void h_fftymanybwd_fg(double *in_r, double *in_c, double *out_r, double *out_c, int aliasing)
//-------------------------------------------------------------------------------------
//
//     Subroutine inverse FFTy transformation, fg style
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, May 2020
//
//-------------------------------------------------------------------------------------
{
	  dblk = 128;

	  ok = ok + cudaMemcpy(psir_d, in_r, sizeof(double)*spxpsi*fpzpsi*npsiy, cudaMemcpyHostToDevice);
	  ok = ok + cudaMemcpy(psic_d, in_c, sizeof(double)*spxpsi*fpzpsi*npsiy, cudaMemcpyHostToDevice);
	  if (ok!=0) printf("Error in copying to device!! fg\n");

	  k_merge_cmp<<<(spxpsi*npsiy*fpzpsi+dblk-1)/dblk,dblk>>>(d_phic,psir_d,psic_d,spxpsi*npsiy*fpzpsi);
	  if (ok!=0) printf("===============>error in call kernel k_merge_cmp not called!! fg\n");


	  //BACKWARDS FFT IN Y //in place
	  cufftExecZ2Z(plan_y_many_psi, d_phic, d_phic, CUFFT_INVERSE);
	  ok = cudaGetLastError();
	  if (ok!=0) printf("Error in cufftExecZ2Z inverse ffty many!! fg\n");


	  //copy back to host
	  k_sec_separate<<<(spxpsi*npsiy*fpzpsi+dblk-1)/dblk,dblk>>>(d_phic, psir_d, psic_d, spxpsi*npsiy*fpzpsi);
	  ok = cudaGetLastError();
	  if(ok!=0) printf("===============>error in call kernel k_sec_separate not called!! fg\n");

	  ok = ok + cudaMemcpy(out_r, psir_d, sizeof(double)*spxpsi*fpzpsi*npsiy, cudaMemcpyDeviceToHost);
	  ok = ok + cudaMemcpy(out_c, psic_d, sizeof(double)*spxpsi*fpzpsi*npsiy, cudaMemcpyDeviceToHost);
	  if (ok!=0) printf("Error in copying to host!!\n");

}//end subroutine h_fftymanybwd_fg
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
