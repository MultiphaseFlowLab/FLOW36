#include <cuda.h>
#include <stdio.h>
#include <cufft.h>
#include <complex.h>

#include "cuda_variables.h"
#include "free_gpu.h"

extern "C" void h_free_gpu(int spx, int spy, int nz, int nx)
//-------------------------------------------------------------------------------------
//
//     Reset the GPU and free allocated memory
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{

#define phiflag phicompflag
#define psiflag psicompflag
#define expx expansionx
#define expy expansiony
#define expz expansionz


  ok = 0;
  ok = ok + cudaFree(d_batch);
  ok = ok + cudaFree(d_batch_c);

  ok = ok + cudaFree(ur_d);
  ok = ok + cudaFree(uc_d);

  ok = ok + cudaFree(d_uopr);
  ok = ok + cudaFree(d_uopc);

  if (ok!=0) printf("Failure in freeing the arrays!!\n");

#if phiflag == 1 || psiflag == 1 || expx != 1 || expy != 1 || expz != 1
  //fg
  ok = ok + cudaFree(d_phic);
  ok = ok + cudaFree(d_phic_c);
  ok = ok + cudaFree(psir_d);
  ok = ok + cudaFree(psic_d);
  ok = ok + cudaFree(d_psir);
  ok = ok + cudaFree(d_psic);
  if (ok!=0) printf("Failure in freeing the fg arrays!!\n");

  ok = ok + cufftDestroy(plan_x_psi);
  ok = ok + cufftDestroy(plan_z_psi);
  ok = ok + cufftDestroy(plan_x_fwd_psi);
  ok = ok + cufftDestroy(plan_y_many_psi);
  if (ok!=0) printf("Failure in fg plans destruction!!\n");

#endif

  ok = 0;
  ok = ok + cufftDestroy(plan_x);
  ok = ok + cufftDestroy(plan_y);
  ok = ok + cufftDestroy(plan_z);
  ok = ok + cufftDestroy(plan_x_fwd);
  ok = ok + cufftDestroy(plan_y_many);

  if (ok!=0) printf("Error in planes destruction!!\n");

  cudaDeviceReset();
}//end subroutine h_free_gpu
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
