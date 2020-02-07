#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cufft.h>
#include <complex.h>
#include <mpi.h>

#define VAR_DECLS
#include "cuda_variables.h"
#include "init_gpu.h"

extern "C" void h_initialize_gpu(int spx_f, int spy_f, int nz_f, int nx_f, MPI_Fint *FLOW_COMM)
//-------------------------------------------------------------------------------------
//
//     Initialize the GPU
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{

  //------------------------------------check how many GPU devices are available----------------------------
  MPI_Comm Cflow_comm;
  Cflow_comm = MPI_Comm_f2c(*FLOW_COMM);

  MPI_Comm_size(Cflow_comm, &nGPUs);//use proper communicator here!!
  MPI_Comm_rank(Cflow_comm, &idGPU);

  int deviceCount = 0;
  cudaError_t error_id = cudaGetDeviceCount(&deviceCount);
  if (error_id != cudaSuccess)
  {
	  printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)error_id, cudaGetErrorString(error_id));
	  exit(EXIT_FAILURE);
  }

  // This function call returns 0 if there are no CUDA capable devices.
  if (deviceCount == 0)
  {
	  printf("There are no available device(s) that support CUDA\n");
  }
  else
  {
	  if(idGPU==0) printf("Detected %d CUDA Capable device(s)\n", deviceCount);
  }

  if(idGPU==0) printf("\nSay hello to the GPU!\n");



  if (nGPUs==1)
  {
    ilGPU = 0;
  }
  else
  {
	//openmpi
    ilGPU = atoi(getenv("OMPI_COMM_WORLD_LOCAL_RANK"));
  }


  cudaSetDevice(ilGPU);

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, ilGPU);
  if (idGPU == 0) printf("The device %d used for the MPI process %d is %s at %d GHz: \n", ilGPU, idGPU, deviceProp.name,deviceProp.clockRate/1000);




  spx = spx_f;
  spy = spy_f;
  nz  = nz_f;
  nx  = nx_f;

  int d_dim = nx * nz * spy;
  if (spx*spy*2*(nz-1) > nx * nz * spy)
  {
    d_dim = spx * spy * 2 * (nz-1);
  }


  //allocate on GPU
  //integer arrays
//  ok = ok + cudaMalloc((void **)&dz_d,      sizeof(double)*nz);
//  ok = ok + cudaMalloc((void **)&fstart_d, sizeof(int)*nz);


  //double arrays
  ok = ok + cudaMalloc((void **)&ur_d,    sizeof(double)*d_dim);
  ok = ok + cudaMalloc((void **)&uc_d,    sizeof(double)*d_dim);
//  ok = ok + cudaMalloc((void **)&d_uout,  sizeof(double)*nx*spy*nz);

  ok = ok + cudaMalloc((void **)&d_uopr,  sizeof(double)*d_dim);
  ok = ok + cudaMalloc((void **)&d_uopc,  sizeof(double)*d_dim);
//  ok = ok + cudaMalloc((void **)&d_uext,  sizeof(double)*d_dim);

//  ok = ok + cudaMalloc((void **)&vr_d,    sizeof(double)*spx*nz*spy);
//  ok = ok + cudaMalloc((void **)&vc_d,    sizeof(double)*spx*nz*spy);
//  ok = ok + cudaMalloc((void **)&d_vout,  sizeof(double)*nx*spy*nz);

//  ok = ok + cudaMalloc((void **)&wr_d,    sizeof(double)*spx*nz*spy);
//  ok = ok + cudaMalloc((void **)&wc_d,    sizeof(double)*spx*nz*spy);
//  ok = ok + cudaMalloc((void **)&d_wout,  sizeof(double)*nx*spy*nz);
  if (ok!=0) printf ("Allocation of double arrays failed\n");

  //allocate output arrays for Chebyshev inverse DCT
  ok = ok + cudaMalloc((void **)&d_batch,  sizeof(cufftDoubleComplex)*nz*spx*spy);
  ok = ok + cudaMalloc((void **)&d_batch_c,sizeof(cufftDoubleComplex)*nz*spx*spy);
  if (ok!=0) printf("Output arrays for Chebyshev not allocated correctly!!\n");

  //mixed products
//  ok = ok + cudaMalloc((void **)&d_uu,   sizeof(double)*nx*nz*spy);
//  ok = ok + cudaMalloc((void **)&d_uuc,  sizeof(cufftDoubleComplex)*nz*spx*spy);
//  ok = ok + cudaMalloc((void **)&d_uuc2,  sizeof(cufftDoubleComplex)*nz*spx*spy);

  if (ok!=0) printf("Failed allocation of mixed product arrays!\n");


  //
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);


  cudaEventRecord(start);
  //create plans for cufft transformations
  cufftPlan1d(&plan_z, (nz-1)*2, CUFFT_D2Z, spx*spy);
  cufftPlan1d(&plan_y,      spy, CUFFT_Z2Z, spx*nz);
  cufftPlan1d(&plan_x,       nx, CUFFT_Z2D, spy*nz);
  cufftPlan1d(&plan_x_fwd,  nx, CUFFT_D2Z, spy*nz);
//  cufftPlan1d(&plan_y_fwd, spy, CUFFT_Z2Z, spx*nz);
  ok = cudaGetLastError();
  if (ok!=0) printf("Error in creating the plans !!\n");
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  float milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);
  printf("Time elapsed in cufft planes creation: %f milliseconds\n",milliseconds);



  // CUFFT plan in y for data of order xzy
//  int datasize=spy;//
  int rank = 1; // --- 1D FFTs
  int n[] = { spy }; // --- Size of the Fourier transform
  int istride = spx*nz, ostride = spx*nz; // --- Distance between two successive input/output elements
  int idist = 1, odist = 1; // --- Distance between batches
  int inembed[] = { 0 }; // --- Input size with pitch (ignored for 1D transforms)
  int onembed[] = { 0 }; // --- Output size with pitch (ignored for 1D transforms)
  int batch = spx*nz; // --- Number of batched executions

  ok = ok + cufftPlanMany(&plan_y_many, rank, n,
                          inembed, istride, idist,
                          onembed, ostride, odist, CUFFT_Z2Z, batch);

  if (ok!=0) printf ("Error in creating the batched many plan\n");












}// end subroutine h_initialize_gpu
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


