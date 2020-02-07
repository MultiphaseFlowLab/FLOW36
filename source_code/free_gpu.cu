#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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

  ok = 0;
  ok = ok + cudaFree(d_batch);
  ok = ok + cudaFree(d_batch_c);

  ok = ok + cudaFree(ur_d);
  ok = ok + cudaFree(uc_d);
//  ok = ok + cudaFree(d_uout);

  ok = ok + cudaFree(d_uopr);
  ok = ok + cudaFree(d_uopc);
//  ok = ok + cudaFree(d_uext);

  //  ok = ok + cudaFree(vr_d);
//  ok = ok + cudaFree(vc_d);
//  ok = ok + cudaFree(d_vout);

//  ok = ok + cudaFree(wr_d);
//  ok = ok + cudaFree(wc_d);
//  ok = ok + cudaFree(d_wout);

//  ok = ok + cudaFree(d_uu);

//  ok = ok + cudaFree(dz_d);
//  ok = ok + cudaFree(fstart_d);
  if (ok!=0) printf("Failure in freeing the arrays!!\n");

  ok = ok + cufftDestroy(plan_x);
  ok = ok + cufftDestroy(plan_y);
  ok = ok + cufftDestroy(plan_z);
  ok = ok + cufftDestroy(plan_x_fwd);
//  ok = ok + cufftDestroy(plan_y_fwd);
  ok = ok + cufftDestroy(plan_y_many);
//  ok = ok + cufftDestroy(plan_z_fwd);

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
