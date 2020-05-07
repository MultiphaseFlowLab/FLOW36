#include <cuda.h>
#include <stdio.h>
#include <cufft.h>
#include <complex.h>
#include <mpi.h>

#define VAR_DECLS
#include "cuda_variables.h"
#include "init_gpu.h"

extern "C" void h_initialize_gpu(int spx_f, int nx_f, int nsx_f, int npx_f, int fpy_f, int npy_f, int ny_f, int spy_f, int nz_f, int fpz_f, int npz_f,
								 int npsix_f, int fpypsi_f, int fpzpsi_f, int spxpsi_f, int spypsi_f, int npsiz_f, int npsiy_f, MPI_Fint *FLOW_COMM)
//-------------------------------------------------------------------------------------
//
//     Initialize the GPU
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
	 //reference: http://on-demand.gputechconf.com/gtc/2014/presentations/S4236-multi-gpu-programming-mpi.pdf
	//rely on process placement: deviceCount == ranks per node (Piz Daint)
	//should work on marconi100 as well if the processes are distributed equally per socket as there are 2 sockets

	//possible to use environment variables provided by MPI launcher
#ifdef OPENMPI
    ilGPU = atoi(getenv("OMPI_COMM_WORLD_LOCAL_RANK"));
#elif MVAPICH2
    ilGPU = atoi(getenv("MV2_COMM_WORLD_LOCAL_RANK"));
#else
    ilGPU = idGPU % deviceCount;
#endif
  }

  cudaSetDevice(ilGPU);

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, ilGPU);
  if (idGPU == 0) printf("The device %d used for the MPI process %d is %s at %d GHz: \n", ilGPU, idGPU, deviceProp.name,deviceProp.clockRate/1000);




  spx = spx_f;
  nsx = nsx_f;
  nx  = nx_f;
  npx = npx_f;

  spy = spy_f;
  fpy = fpy_f;
  npy = npy_f;
  ny  = ny_f;

  nz  = nz_f;
  fpz = fpz_f;
  npz = npz_f;

  npsix  = npsix_f;
  fpypsi = fpypsi_f;
  fpzpsi = fpzpsi_f;

  spxpsi = spxpsi_f;
  spypsi = spypsi_f;
  npsiz  = npsiz_f;

  npsiy = npsiy_f;



  //calculate max array dim
  int d_dim = nx * fpz * fpy;
  if (spx * spy * 2 * (nz-1) > d_dim) d_dim = spx*spy*2*(nz-1);
  if (spx * fpz * ny         > d_dim) d_dim = spx*fpz*ny;
  if (spx * nz * fpy         > d_dim) d_dim = spx*nz*fpy;

  //allocate on GPU


  //double arrays
  ok = ok + cudaMalloc((void **)&ur_d,    sizeof(double)*d_dim);
  ok = ok + cudaMalloc((void **)&uc_d,    sizeof(double)*d_dim);

  ok = ok + cudaMalloc((void **)&d_uopr,  sizeof(double)*d_dim);
  ok = ok + cudaMalloc((void **)&d_uopc,  sizeof(double)*d_dim);
  if (ok!=0) printf ("Allocation of double arrays failed\n");

  //allocate output arrays for Chebyshev inverse DCT
  ok = ok + cudaMalloc((void **)&d_batch,  sizeof(cufftDoubleComplex)*d_dim);
  ok = ok + cudaMalloc((void **)&d_batch_c,sizeof(cufftDoubleComplex)*d_dim);
  if (ok!=0) printf("Output arrays for Chebyshev not allocated correctly!!\n");

  #if phiflag == 1 || psiflag == 1 || expx != 1 || expy != 1 || expz != 1
  //calculate max array dim
  int dim_psi = npsix * fpzpsi * fpypsi;
  if (spxpsi * spypsi * 2 * (npsiz-1) > dim_psi) dim_psi = spxpsi*spypsi*2*(npsiz-1);
  if (spxpsi * fpzpsi * npsiy         > dim_psi) dim_psi = spxpsi*fpzpsi*npsiy;
  if (spxpsi * npsiz * fpypsi         > dim_psi) dim_psi = spxpsi*npsiz*fpypsi;

  //fg arrays
  //double arrays
  ok = ok + cudaMalloc((void **)&psir_d,    sizeof(double)*dim_psi);
  ok = ok + cudaMalloc((void **)&psic_d,    sizeof(double)*dim_psi);

  ok = ok + cudaMalloc((void **)&d_psir,  sizeof(double)*dim_psi);
  ok = ok + cudaMalloc((void **)&d_psic,  sizeof(double)*dim_psi);
  if (ok!=0) printf ("Allocation of double fg arrays failed\n");

  //allocate output arrays for Chebyshev inverse DCT
  ok = ok + cudaMalloc((void **)&d_phic,  sizeof(cufftDoubleComplex)*dim_psi);
  ok = ok + cudaMalloc((void **)&d_phic_c,sizeof(cufftDoubleComplex)*dim_psi);
  if (ok!=0) printf("Output arrays for Chebyshev fg not allocated correctly!!\n");
  #endif


  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);


  cudaEventRecord(start);
  //create plans for cufft transformations
  cufftPlan1d(&plan_z,      (nz-1)*2,  CUFFT_D2Z,  spx*spy);//OK both ways
  cufftPlan1d(&plan_y,      ny,        CUFFT_Z2Z,  spx*npz);//never tested
  cufftPlan1d(&plan_x,      nx,  CUFFT_Z2D,  fpy*fpz);//OK so npz=fpz and npy=fpy
  cufftPlan1d(&plan_x_fwd,  nx,        CUFFT_D2Z,  fpy*fpz);//OK so ...



  // CUFFT plan in y for data of order xzy
//  int datasize=spy;//
  int rank = 1; // --- 1D FFTs
  int n[] = { ny }; // --- Size of the Fourier transform
  int istride = spx*npz, ostride = spx*npz; // --- Distance between two successive input/output elements
  int idist = 1, odist = 1; // --- Distance between batches
  int inembed[] = { 0 }; // --- Input size with pitch (ignored for 1D transforms)
  int onembed[] = { 0 }; // --- Output size with pitch (ignored for 1D transforms)
  int batch = spx*npz; // --- Number of batched executions

  ok = ok + cufftPlanMany(&plan_y_many, rank, n,
                          inembed, istride, idist,
                          onembed, ostride, odist, CUFFT_Z2Z, batch);

#if phiflag == 1 || psiflag == 1 || expx != 1 || expy != 1 || expz != 1
  cufftPlan1d(&plan_z_psi,  (npsiz-1)*2,  CUFFT_D2Z,  spxpsi*spypsi);//OK both ways
  cufftPlan1d(&plan_x_psi,  npsix,        CUFFT_Z2D,  fpypsi*fpzpsi);//OK so npz=fpz and npy=fpy
  cufftPlan1d(&plan_x_fwd_psi,  npsix,    CUFFT_D2Z,  fpypsi*fpzpsi);//OK so ...

  int n_psi[] = { npsiy }; // --- Size of the Fourier transform
  int istride_psi = spxpsi*fpzpsi, ostride_psi = spxpsi*fpzpsi; // --- Distance between two successive input/output elements
  int batch_psi = spxpsi*fpzpsi; // --- Number of batched executions

  ok = ok + cufftPlanMany(&plan_y_many_psi, rank, n_psi,
                            inembed, istride_psi, idist,
                            onembed, ostride_psi, odist, CUFFT_Z2Z, batch_psi);
#endif

  if (ok!=0) printf ("Error in creating the batched many plan\n");
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  float milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);
  if(idGPU == 0) printf("Time elapsed in cufft planes creation: %f milliseconds\n",milliseconds);


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


