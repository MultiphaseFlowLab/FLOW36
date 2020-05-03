#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
#include <cufft.h>

#define VAR_DECLS
#include "cuda_batch.h"

void h_batch(int nx,int ny,cufftDoubleComplex* h_signal)
{

  int i;
  int signal_size; 
  int ok = 0;


  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  signal_size = nx * ny;
  printf ("totale size %d; %d 1D FFTs of %d elements\n",signal_size,ny,nx);
  long mem_size = sizeof(cufftDoubleComplex) * signal_size;
  cufftDoubleComplex* d_signal;
  cufftDoubleComplex* d_out;
  
  ok = ok + cudaMalloc((void**)&d_signal, mem_size);
  ok = ok + cudaMalloc((void**)&d_out,    mem_size);
  if (ok != 0) printf("Error in allocating the GPU array!\n");
  
  ok = ok + cudaMemcpy(d_signal, h_signal, mem_size,cudaMemcpyHostToDevice);
  if (ok != 0) printf ("Error in coping the array to device\n");

  // CUFFT plan
  cufftHandle plan;
  int datasize=nx;//
  int rank = 1; // --- 1D FFTs
  int n[] = { nx }; // --- Size of the Fourier transform
  int istride = 1, ostride = 1; // --- Distance between two successive input/output elements
  int idist = nx, odist = nx; // --- Distance between batches
  int inembed[] = { 0 }; // --- Input size with pitch (ignored for 1D transforms)
  int onembed[] = { 0 }; // --- Output size with pitch (ignored for 1D transforms)
  int batch = ny; // --- Number of batched executions

  ok = ok + cufftPlanMany(&plan, rank, n,
                          inembed, istride, idist,
                          onembed, ostride, odist, CUFFT_Z2Z, ny);

  if (ok!=0) printf ("Error in creating the batched many plan\n");

  cufftHandle plan2;
  cufftPlan1d(&plan2, nx, CUFFT_Z2Z, ny);
  if (ok!=0) printf ("Error in creating the batched 1d plan\n");



  //warm-up
  ok = ok + cufftExecZ2Z(plan, (cufftDoubleComplex *)d_signal, (cufftDoubleComplex *)d_out, CUFFT_FORWARD);
  float milliseconds;
  float time_all=0.0;
  // Transform signal
  printf("Transforming signal cufftExecC2C\n");
  for (i=0;i<100;i++)
  {
    cudaEventRecord(start);
    ok = ok + cufftExecZ2Z(plan, (cufftDoubleComplex *)d_signal, (cufftDoubleComplex *)d_out, CUFFT_FORWARD); 
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    time_all+=milliseconds;
  }
  milliseconds=time_all/100.0;
  time_all=0.0;
  printf("Effective CUFFT-FORWARD Bandwidth (GB/s): %f\nMeasured time (ms): %f\n", signal_size*sizeof(cufftDoubleComplex)/milliseconds/1e6,milliseconds);

  if (ok != 0) printf ("Error in forward transformation\n");

  // Transform signal back
  printf("Transforming signal back cufftExecC2C\n");
  for (i=0;i<100;i++)
  {
    cudaEventRecord(start);
    ok = ok + (cufftExecZ2Z(plan, (cufftDoubleComplex *)d_out, (cufftDoubleComplex *)d_signal, CUFFT_INVERSE));
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    time_all+=milliseconds;
  }
  milliseconds=time_all/100;
  printf("Effective CUFFT-BACKWARD Bandwidth (GB/s): %f\nMeasured time (ms): %f\n", signal_size*sizeof(cufftDoubleComplex)/milliseconds/1e6,milliseconds);

  if (ok != 0) printf ("Error in backwards transformation\n");

  // Copy device memory to host
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

    printf("\n GPU memory usage: used = %5.2f GB, total (free+used) = %5.2f GB, % = %5.2f\n", used_db/(1073741824), total_db/(1073741824), used_db/total_db*100);
  }


  cufftDoubleComplex* h_inverse_signal = (cufftDoubleComplex*)malloc(sizeof(cufftDoubleComplex) * signal_size);;
  ok = ok + (cudaMemcpy(h_inverse_signal, d_signal, mem_size,cudaMemcpyDeviceToHost));
  if (ok != 0) printf ("Error in coping back data to host\n");
  
  
  for(i=0;i< signal_size;i=i+100)
  {
    h_inverse_signal[i].x= h_inverse_signal[i].x/(double)nx;//signal_size;
    h_inverse_signal[i].y= h_inverse_signal[i].y/(double)nx;//signal_size;

    printf("first : %24.16le %24.16le  after %24.16le %24.16le \n",h_signal[i].x,h_signal[i].y,h_inverse_signal[i].x,h_inverse_signal[i].y);
  }  



  //Destroy CUFFT context
  ok = ok +(cufftDestroy(plan));
  ok = ok + cufftDestroy(plan2);
  if (ok != 0) printf ("destroy plane failed\n");
  // cleanup memory
  free(h_signal);

  free(h_inverse_signal);
  ok = ok +(cudaFree(d_signal));
  if (ok != 0) printf ("cudaFree failed\n");
  cudaDeviceReset();

}

