#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
//##include <mpi.h>
#include <string.h>
#include <math.h>
#include <cufft.h>
#include <complex.h>

#include "cuda_variables.h"
#include "cuda_spec_phys.h"
//#include "cuda_surf.h"
#include "cuda_tran.h"

#define MODULO(X,Y) (((X)%(Y))>=0) ? ((X)%(Y)) : ((X)%(Y) + (Y))



extern "C" void h_ffty_many_bwd(double *in_r, double *in_c, double *out_r, double *out_c, int aliasing)
//-------------------------------------------------------------------------------------
//
//     Subroutine h_ffty_many_bwd for spectral to physical FFTy without data transposition
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{
  dblk = 128;

  ok = ok + cudaMemcpy(ur_d,in_r,sizeof(double)*spx*npz*ny,cudaMemcpyHostToDevice);
  ok = ok + cudaMemcpy(uc_d,in_c,sizeof(double)*spx*npz*ny,cudaMemcpyHostToDevice);
  if (ok!=0) printf("Error in copying to device!!\n");
  k_merge_cmp<<<(spx*ny*npz+dblk-1)/dblk,dblk>>>(d_batch,ur_d,uc_d,spx*ny*npz);
  if (ok!=0) printf("===============>error in call kernel k_merge_cmp not called!! ffty_back\n");
//  ok = ok + cudaMemcpy2D(d_batch,       2 * sizeof(d_batch), in_r, sizeof(in_r), sizeof(in_r), spx*nz*spy, cudaMemcpyHostToDevice);
//  ok = ok + cudaMemcpy2D(&d_batch[0].y, 2 * sizeof(d_batch), in_c, sizeof(in_c), sizeof(in_c), spx*nz*spy, cudaMemcpyHostToDevice);
//  if (ok!=0) printf("Error in copying to device!!\n");


  //BACKWARDS FFT IN Y //in place
  cufftExecZ2Z(plan_y_many, d_batch, d_batch, CUFFT_INVERSE);
  ok = cudaGetLastError();
  if (ok!=0) printf("Error in cufftExecZ2Z inverse ffty many!!\n");


  //copy back to host
  k_sec_separate<<<(spx*ny*npz+dblk-1)/dblk,dblk>>>(d_batch, ur_d, uc_d, spx*ny*npz);
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call kernel k_sec_separate not called!! fftx_fwd\n");
  ok = ok + cudaMemcpy(out_r,ur_d,sizeof(double)*spx*npz*ny,cudaMemcpyDeviceToHost);
  ok = ok + cudaMemcpy(out_c,uc_d,sizeof(double)*spx*npz*ny,cudaMemcpyDeviceToHost);

//  ok = ok + cudaMemcpy2D(out_r, sizeof(out_r),d_batch_c,       2 * sizeof(d_batch_c),sizeof(d_batch_c), spx*nz*spy, cudaMemcpyDeviceToHost);
//  ok = ok + cudaMemcpy2D(out_c, sizeof(out_c),&d_batch_c[0].y, 2 * sizeof(d_batch_c),sizeof(d_batch_c), spx*nz*spy, cudaMemcpyDeviceToHost);
  if (ok!=0) printf("Error in copying to host!!\n");
}


extern "C" void h_chebyshev_back(double *in_r, double *in_c, double *out_r, double *out_c, int aliasing)
//-------------------------------------------------------------------------------------
//
//     Subroutine h_chebyshev_back for spectral to physical DCT
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

  ok = ok + cudaMemcpy(ur_d,in_r,sizeof(double)*spx*nz*spy,cudaMemcpyHostToDevice);
  ok = ok + cudaMemcpy(uc_d,in_c,sizeof(double)*spx*nz*spy,cudaMemcpyHostToDevice);
  if (ok!=0) printf("Error in copying to device!!\n");


  //grid,block sizes for xzy to zxy transposition
  int TILE_DIM=32;
  int BLOCK_ROWS=8;
  nbx = spx/TILE_DIM;
  nby = nz/TILE_DIM;
  if (spx % TILE_DIM != 0) nbx = nbx + 1;
  if (nz  % TILE_DIM != 0) nby = nby + 1;
  dim3 tgrid1(nbx,nby,spy);
  dim3 tBlock(TILE_DIM, BLOCK_ROWS,1);

  //transpose from xzy to zxy to gain coalescence in z
  k_t102<<<tgrid1,tBlock>>>(d_uopr, ur_d, spx, nz, spy);//#OUTPUT,INPUT,SIZES
  k_t102<<<tgrid1,tBlock>>>(d_uopc, uc_d, spx, nz, spy);
  ok = cudaGetLastError();
  if (ok!=0) printf ("===============>error in call kernel k_t102 not called!! \n");



  //perform inplace aliasing
  if (aliasing == 1)
  {
	  int ali = floor(2.0*double(nz)/3.0); //arrays start from zero

	  k_alias_1st<<<((nz-ali)*spx*spy+dblk-1)/dblk,dblk>>>(d_uopr, d_uopc, ali, nz, spx*spy*nz);//evaluate a merging of the two kernels

//	  k_alias_small<<<((nz-ali)*spx*spy+dblk-1)/dblk,dblk>>>(d_uopr, ali, nz);
//      k_alias_small<<<((nz-ali)*spx*spy+dblk-1)/dblk,dblk>>>(d_uopc, ali, nz);
      ok = cudaGetLastError();
      if (ok!=0) printf("===============>error in call kernel k_alias_1st not called!! \n");
  }


  //manipulate the input arrays
  //WARNING: this kernel does very little and it should be evaluated if the overhead is big. actually I just need one block with nz threads, NOT  nz,(spx*spy+1)/32*32
  k_manip<<<(nz+dblk-1)/dblk,dblk>>>(d_uopr, nz, nz*spx*spy);
  k_manip<<<(nz+dblk-1)/dblk,dblk>>>(d_uopc, nz, nz*spx*spy);
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call kernel k_manip not called!! \n");




  //make input for cufftD2Z even symmetrical ##format OUTPUT,INPUT,sizes
  k_mirr_bigtime<<<(nz*spx*spy+dblk-1)/dblk,dblk>>>(ur_d, uc_d, d_uopr,  d_uopc, nz, 2*(nz-1), spx*spy*nz);

//  k_mirror<<<spx*spy,32*(nz/32+1),nz*sizeof(double)>>>(ur_d, d_uopr, nz, 2*(nz-1));
//  k_mirror<<<spx*spy,32*(nz/32+1),nz*sizeof(double)>>>(uc_d,  d_uopc, nz, 2*(nz-1));
  ok = cudaGetLastError();
  if (ok!=0) printf("===============>error in call kernel k_mirr_bigtime not called!! \n");




  //execute the inverse Chebyshev DCT
  cufftExecD2Z(plan_z, ur_d, d_batch);
  cufftExecD2Z(plan_z, uc_d,  d_batch_c);
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call cufftD2Z not called!! \n");





  //-------------------------------------------------------------------------------------------------------------------------------------------------------------
  //show memory usage @ this point (should be highest)
  size_t free_byte ;
  size_t total_byte ;
  if (check_mem < 1 && idGPU == 0)//add flag for GPU zero in multi-GPU version
  {
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
  }
  check_mem=1;
  //-------------------------------------------------------------------------------------------------------------------------------------------------------------




  //merge the two outputs of the Chebyshev DFT and normalize
  double norm_fac = 0.5e0;
  k_sep_cmp<<<(spx*spy*nz+dblk-1)/dblk,dblk>>>(d_batch,   ur_d, norm_fac, nz*spx*spy);
  k_sep_cmp<<<(spx*spy*nz+dblk-1)/dblk,dblk>>>(d_batch_c, uc_d, norm_fac, nz*spx*spy);
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call kernel k_sec_copy not called!! \n");



  //transpose back from zxy to xzy
  dim3 tgrid2(nby,nbx,spy);
//  k_cmp_t102<<<tgrid2,tBlock>>>(d_batch_c, d_batch, nz, spx, spy);//#OUTPUT,INPUT,SIZES

  k_t102<<<tgrid2,tBlock>>>(d_uopr,ur_d,nz,spx,spy);
  k_t102<<<tgrid2,tBlock>>>(d_uopc,uc_d,nz,spx,spy);
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call kernel k_transpose_back not called!! cheb_back\n");

  ok = ok + cudaMemcpy(out_r,d_uopr,sizeof(double)*spx*nz*spy,cudaMemcpyDeviceToHost);
  ok = ok + cudaMemcpy(out_c,d_uopc,sizeof(double)*spx*nz*spy,cudaMemcpyDeviceToHost);
//  ok = ok + cudaMemcpy2D(out_r, sizeof(out_r),d_batch_c,       2 * sizeof(d_batch_c),sizeof(d_batch_c), spx*nz*spy, cudaMemcpyDeviceToHost);
//  ok = ok + cudaMemcpy2D(out_c, sizeof(out_c),&d_batch_c[0].y, 2 * sizeof(d_batch_c),sizeof(d_batch_c), spx*nz*spy, cudaMemcpyDeviceToHost);
  if (ok!=0) printf("Error in copying to host!!\n");

}//end subroutine h_chebyshev_back
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
extern "C" void h_ffty_back(double *in_r, double *in_c, double *out_r, double *out_c, int aliasing)
//-------------------------------------------------------------------------------------
//
//     Subroutine h_ffty_back for fft from spectral to physical
//     WARNING: same as h_chebyshev_back
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{
  dblk = 128;

  //copy to device
  ok = ok + cudaMemcpy(ur_d,in_r,sizeof(double)*spx*ny*npz,cudaMemcpyHostToDevice);
  ok = ok + cudaMemcpy(uc_d,in_c,sizeof(double)*spx*ny*npz,cudaMemcpyHostToDevice);
  if (ok!=0) printf("Error in copying to device!!\n");

  //merge to cufftDoubleComplex
  k_merge_cmp<<<(spx*ny*npz+dblk-1)/dblk,dblk>>>(d_batch,ur_d,uc_d,spx*ny*npz);
  if (ok!=0) printf("===============>error in call kernel k_merge_cmp not called!! ffty_back\n");
//  ok = ok + cudaMemcpy2D(d_batch,       2 * sizeof(d_batch), in_r, sizeof(in_r), sizeof(in_r), spx*nz*spy, cudaMemcpyHostToDevice);
//  ok = ok + cudaMemcpy2D(&d_batch[0].y, 2 * sizeof(d_batch), in_c, sizeof(in_c), sizeof(in_c), spx*nz*spy, cudaMemcpyHostToDevice);
//  if (ok!=0) printf("Error in copying to device!!\n");



  //grid, block for xzy to yzx
  int TILE_DIM=32;
  int BLOCK_ROWS=8;
  nbx = spx/TILE_DIM;
  nby = ny/TILE_DIM;
  if (spx   % TILE_DIM != 0) nbx = nbx + 1;
  if (ny  % TILE_DIM != 0) nby = nby + 1;
  dim3 tgrid1(nbx,nby,npz);
  dim3 tBlock(TILE_DIM, BLOCK_ROWS,1);


  //transpose to gain coalescence in y
  k_cmp_t210<<<tgrid1,tBlock>>>(d_batch_c, d_batch, spx, npz, ny);//#OUTPUT,INPUT,sizes
  ok = cudaGetLastError();
  if(ok!=0) printf("error in call kernel k_cmp_t210 not called!! \n");



  //perform inplace aliasing
  if (aliasing == 1)
  {
	int al_low  = floor(2.0/3.0*(ny/2+1));//first aliased position
	int al_up = ny-floor(2.0/3.0*(ny/2));//last aliased position
	k_cmp_alias<<<((al_up-al_low+1)*spx*npz+dblk-1)/dblk,dblk>>>(d_batch_c, al_low, al_up-al_low, spy, spx*ny*npz);
	ok = cudaGetLastError();
    if (ok!=0) printf("===============>error in call kernel k_cmp_alias_y not called!! \n");
  }



  //BACKWARDS FFT IN Y
  cufftExecZ2Z(plan_y, d_batch_c, d_batch, CUFFT_INVERSE);
  ok = cudaGetLastError();
  if (ok!=0) printf("Error in cufftExecZ2Z inverse ffty!!\n");


  //Normalize
  int tot_size = ny * spx * npz;
  double norm_fact = 1.0e0 / (double)ny;
  k_norm_cmp<<<(ny*spx*npz+dblk-1)/dblk,dblk>>>(d_batch, norm_fact, tot_size);
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call kernel k_norm_cmp not called!! \n");


  //transpose back for Fortran
  dim3 tgrid2(nby,nbx,npz);
  k_cmp_t210<<<tgrid2,tBlock>>>(d_batch_c, d_batch, ny, npz, spx);//#OUTPUT,INPUT,sizes
  ok = cudaGetLastError();
  if(ok!=0) printf("error in call kernel k_t210 not called!! \n");

  //copy back to host
  k_sec_separate<<<(spx*ny*npz+dblk-1)/dblk,dblk>>>(d_batch_c,ur_d,uc_d, spx*ny*npz);
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call kernel k_sec_separate not called!! fftx_fwd\n");
  ok = ok + cudaMemcpy(out_r,ur_d,sizeof(double)*spx*npz*ny,cudaMemcpyDeviceToHost);
  ok = ok + cudaMemcpy(out_c,uc_d,sizeof(double)*spx*npz*ny,cudaMemcpyDeviceToHost);

//  ok = ok + cudaMemcpy2D(out_r, sizeof(out_r),d_batch_c,       2 * sizeof(d_batch_c),sizeof(d_batch_c), spx*nz*spy, cudaMemcpyDeviceToHost);
//  ok = ok + cudaMemcpy2D(out_c, sizeof(out_c),&d_batch_c[0].y, 2 * sizeof(d_batch_c),sizeof(d_batch_c), spx*nz*spy, cudaMemcpyDeviceToHost);
  if (ok!=0) printf("Error in copying to host!!\n");

}//end subroutine h_ffty_back
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
extern "C" void h_fftx_back(double *in_r, double *in_c, double *out_r, int aliasing)
//-------------------------------------------------------------------------------------
//
//     Subroutine h_fftx_back for fft from spectral to physical
//     WARNING: same as h_chebyshev_back
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{
  dblk = 128;

  //copy to device
  ok = ok + cudaMemcpy(ur_d,in_r,sizeof(double)*(nx/2+1)*npy*npz,cudaMemcpyHostToDevice);
  ok = ok + cudaMemcpy(uc_d,in_c,sizeof(double)*(nx/2+1)*npy*npz,cudaMemcpyHostToDevice);
  if (ok!=0) printf("Error in copying to device!!\n");
  k_merge_cmp<<<((nx/2+1)*npy*npz+dblk-1)/dblk,dblk>>>(d_batch,ur_d,uc_d,(nx/2+1)*npy*npz);
  if (ok!=0) printf("===============>error in call kernel k_merge_cmp not called!! ffty_back\n");

//  ok = ok + cudaMemcpy2D(d_batch,       2 * sizeof(d_batch), in_r, sizeof(in_r), sizeof(in_r), spx*nz*spy, cudaMemcpyHostToDevice);
//  ok = ok + cudaMemcpy2D(&d_batch[0].y, 2 * sizeof(d_batch), in_c, sizeof(in_c), sizeof(in_c), spx*nz*spy, cudaMemcpyHostToDevice);
//  if (ok!=0) printf("Error in copying to device!!\n");

  if (aliasing == 1)
  {
    int al_low = floor(2.0/3.0*(nx/2+1)); //arrays start from zero
	int al_up  = nx/2;//corresponds to nz-1, alias the last position of the c array to be set to zero
	k_cmp_alias<<<((al_up-al_low+1)*npy*npz+dblk-1)/dblk,dblk>>>(d_batch, al_low, al_up-al_low+1, (nx/2+1), (nx/2+1)*npy*npz);
	ok = cudaGetLastError();
    if (ok!=0) printf("===============>error in second call kernel k_cmp_alias_y not called!! \n");
  }



  //BACKWARD FFT to physical space
  cufftExecZ2D(plan_x, d_batch, ur_d);


  //Normalize
  double norm_fact = 1.0e0/(double)nx;
  int tot_size  = nx*fpy*fpz;
  k_norm_real<<<(fpy*nx*fpz+dblk-1)/dblk,dblk>>>(ur_d, norm_fact, tot_size);



  //copy back to host
  ok = ok + cudaMemcpy(out_r,ur_d,sizeof(double)*nx*fpz*fpy,cudaMemcpyDeviceToHost);
  if (ok!=0) printf("Error in copying back to host!\n");



}//end subroutine h_fftx_back
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
void cuda_spectral_to_phys(double *real_d,double *imm_d,double *out_d)
//-------------------------------------------------------------------------------------
//
//     Subroutine spectral to physical for single GPU
//     WARNING: memory management needs improvement, could consider working on in-place transposition
//     or batched FFTs that need no transposition, despite the current approach grants coalescence for
//     the other kernels
//
//     Next step will be to split the function between the three inverse transformations
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{
  //INVERSE DCT in z

  //grid,block sizes for xzy to zxy transposition
  int TILE_DIM=32;
  int BLOCK_ROWS=8;
  nbx = spx/TILE_DIM;
  nby = nz/TILE_DIM;
  if (spx % TILE_DIM != 0) nbx = nbx + 1;
  if (nz  % TILE_DIM != 0) nby = nby + 1;
  dim3 tgrid1(nbx,nby,spy);
  dim3 tBlock(TILE_DIM, BLOCK_ROWS,1);

  //transpose from xzy to zxy to gain coalescence in z
  k_t102<<<tgrid1,tBlock>>>(d_uopr,  real_d,
		                    spx, nz, spy);//#OUTPUT,INPUT,SIZES
  k_t102<<<tgrid1,tBlock>>>(d_uopc,  imm_d,
		                    spx, nz, spy);
  ok = cudaGetLastError();
  if (ok!=0) printf ("===============>error in call kernel k_t102 not called!! \n");


  //perform inplace aliasing
  if (aliasing == 1)
  {
	  int ali = floor(2.0*double(nz)/3.0); //arrays start from zero

	  k_alias_small<<<((nz-ali)*spx*spy+dblk-1)/dblk,dblk>>>(d_uopr,
			                                                 ali, nz);
      k_alias_small<<<((nz-ali)*spx*spy+dblk-1)/dblk,dblk>>>(d_uopc,
    		                                                 ali, nz);
      ok = cudaGetLastError();
      if (ok!=0) printf("===============>error in call kernel k_alias_small not called!! \n");
  }


  //manipulate the input arrays
  //WARNING: this kernel does very little and it should be evaluated if the overhead is big. actually I just need one block with nz threads, NOT  nz,(spx*spy+1)/32*32
  k_manip<<<(nz+dblk-1)/dblk,dblk>>>(d_uopr,
		                             nz,    nz*spx*spy);
  k_manip<<<(nz+dblk-1)/dblk,dblk>>>(d_uopc,
		                             nz,    nz*spx*spy);
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call kernel k_manip not called!! \n");



//wARNING: tutti i kernel accoppiati possono essere ridotti ad un unico kernel che ha il doppio delle thread ed esegue su due matrici
  //make input for cufftD2Z even symmetrical ##format OUTPUT,INPUT,sizes
  k_mirror<<<spx*spy,32*(nz/32+1),nz*sizeof(double)>>>(real_d, d_uopr,
		                                               nz,     2*(nz-1));
  k_mirror<<<spx*spy,32*(nz/32+1),nz*sizeof(double)>>>(imm_d,  d_uopc,
		                                               nz,     2*(nz-1));
  ok = cudaGetLastError();
  if (ok!=0) printf("===============>error in call kernel k_mirror not called!! \n");




  //execute the inverse Chebyshev DCT
  cufftExecD2Z(plan_z, real_d, d_batch);
  cufftExecD2Z(plan_z, imm_d,  d_batch_c);
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call cufftD2Z not called!! \n");



  //merge the two outputs of the Chebyshev DFT and normalize
  k_sec_copy<<<(spx*spx*nz+dblk-1)/dblk,dblk>>>(d_batch, d_batch_c,
		                                        nz*spx*spy);//#OUTPUT,INPUT,size -spx*spy,32*(nz/32+1)
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call kernel k_sec_copy not called!! \n");






  //INVERSE FFT in y
  //set up grid,blocks for zxy to yxz transposition
  nbx = nz/TILE_DIM;
  nby = spy/TILE_DIM;
  if (nz   % TILE_DIM != 0) nbx = nbx + 1;
  if (spy  % TILE_DIM != 0) nby = nby + 1;
  dim3 tgrid2(nbx,nby,spx);

  //transpose to gain coalescence in y zxy --> yxz
  k_cmp_t210<<<tgrid2,tBlock>>>(d_batch_c, d_batch,
		                        nz, spx, spy);//#OUTPUT,INPUT,sizes
  ok = cudaGetLastError();
  if(ok!=0) printf("error in call kernel k_t210 not called!! \n");




  //perform inplace aliasing
  if (aliasing == 1)
  {
	int al_low  = floor(2.0/3.0*(spy/2+1));//first aliased position
	int al_up = spy-floor(2.0/3.0*(spy/2));//last aliased position
	k_cmp_alias<<<((al_up-al_low+1)*spx*nz+dblk-1)/dblk,dblk>>>(d_batch_c,
			                                                      al_low, al_up-al_low, spy, spx*spy*nz);
	ok = cudaGetLastError();
    if (ok!=0) printf("===============>error in call kernel k_cmp_alias_y not called!! \n");
  }




  //BACKWARDS FFT IN Y
  cufftExecZ2Z(plan_y, d_batch_c, d_batch, CUFFT_INVERSE);
  ok = cudaGetLastError();
  if (ok!=0) printf("Error in cufftExecZ2Z inverse ffty!!\n");

  //Normalize
  int tot_size = spy * spx * nz;
  double norm_fact = 1.0e0 / (double)spy;
  k_norm_cmp<<<(spy*spx*nz+dblk-1)/dblk,dblk>>>(d_batch,
		                                        norm_fact, tot_size);
  ok = cudaGetLastError();
  if(ok!=0) printf("===============>error in call kernel k_norm_cmp not called!! \n");





  //INVERSE FFT in x

  //grid,block sizes for yxz to xyz transposition: WARNING: little repetition (could use the maximum between the two to create the grid??
  nbx = spy/TILE_DIM;
  nby = spx/TILE_DIM;
  if (spy  % TILE_DIM != 0) nbx = nbx + 1;
  if (spx  % TILE_DIM != 0) nby = nby + 1;
  dim3 tgrid3(nbx,nby,nz);


  //transpose to gain coalescence in x:  yxz --> xyz
  k_cmp_t102<<<tgrid3,tBlock>>>(d_batch_c, d_batch,
		                        spy, spx, nz);//##OUTPUT,INPUT,sizes
  ok = cudaGetLastError();
  if (ok!=0) printf("===============>error in call kernel k_cmp_t102 not called!! \n");


  if (aliasing == 1)
  {
    int al_low = floor(2.0/3.0*(nx/2+1)); //arrays start from zero
	int al_up  = nx/2;//corresponds to nz-1, alias the last position of the c array to be set to zero
	k_cmp_alias<<<((al_up-al_low+1)*spy*nz+dblk-1)/dblk,dblk>>>(d_batch_c,
			                                                      al_low, al_up-al_low+1, spx, spx*spy*nz);
	ok = cudaGetLastError();
    if (ok!=0) printf("===============>error in second call kernel k_cmp_alias_y not called!! \n");
  }


  //BACKWARD FFT to physical space
  cufftExecZ2D(plan_x, d_batch_c, out_d);

  //Normalize
  norm_fact = 1.0e0/(double)nx;
  tot_size  = nx*spy*nz;
  k_norm_real<<<(spy*nx*nz+dblk-1)/dblk,dblk>>>(out_d,
		                                        norm_fact, tot_size);




  //OUTPUT is ordered as xyz


}//end subroutine cuda_spectral_to_phys
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
