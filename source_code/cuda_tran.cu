#include <cuda.h>
#include <cufft.h>
#include <cooperative_groups.h>
#include <cuda_runtime.h>
#include <complex.h>
#include "cuda_variables.h"
#include "cuda_tran.h"


#define TILE_DIM   32
#define BLOCK_ROWS 8

__global__ void k_alias_3rd(cufftDoubleComplex *a, int al_low, int al_up, int spx, int ny, int npz)
//-------------------------------------------------------------------------------------
//
//     Perform aliasing in third direction of array:
//	   NEW: 3D grid-block structures, simplifies a lot!

//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, May 2020
//
//-------------------------------------------------------------------------------------
{
  int x_index = blockIdx.x*blockDim.x + threadIdx.x; //absolute x,y,z indexes
  int y_index = blockIdx.y*blockDim.y + threadIdx.y;
  int z_index = blockIdx.z*blockDim.z + threadIdx.z;

  if (x_index < spx)
  {
    if (z_index < npz)
    {
      if (y_index < al_up)
      {
    	int ind = x_index + spx * (z_index + npz * (al_low+y_index));

    	a[ind].x = 0.0e0;
    	a[ind].y = 0.0e0;
      }
    }
  }

}//end kernel k_alias_3rd
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

void __global__ k_alias_1st_cmp(cufftDoubleComplex *a, int al_low, int nx, int dim)
//-------------------------------------------------------------------------------------
//
//     perform aliasing on cufftDoubleComplex velocity array in first dimension
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, March 2020
//
//-------------------------------------------------------------------------------------
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;   //absolute thread index
  int check = al_low + index/(nx-al_low)*nx + index % (nx-al_low);

  //each thread accesses only the values to be modified
  if (check < dim)
  {
    a[check].x = 0.0e0;
    a[check].y = 0.0e0;
  }

}//end kernel k_cmp_alias_small
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
void __global__ k_alias_1st(double *a, double *b, int al_low, int nx, int dim)
//-------------------------------------------------------------------------------------
//
//     perform aliasing on double velocity arrays in first dimension
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, March 2020
//
//-------------------------------------------------------------------------------------
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;   //absolute thread index
  int check = al_low + index/(nx-al_low)*nx + index % (nx-al_low);

  //each thread accesses only the values to be modified
  if (check < dim)
  {
    a[check] = 0.0e0;
    b[check] = 0.0e0;
  }

}//end kernel k_cmp_alias_small
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
__global__ void k_merge_cmp(cufftDoubleComplex *out,
		                    double *re, double *im,
		                    int size)
//-------------------------------------------------------------------------------------
//
//     This kernel merges Re and Im inputs to a cufftDoubleComplex array
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{
  int index = blockIdx.x*blockDim.x + threadIdx.x; //absolute thread index
  if (index<size)
  {
	out[index].x = re[index];
    out[index].y = im[index];
  }
}//end kernel k_merge_cmp
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
__global__ void k_sec_separate(cufftDoubleComplex *in,
		                       double *re, double *im,
		                       int size)
//-------------------------------------------------------------------------------------
//
//     this kernel separates the input R I R I R I into two double arrays R R R and I I I
//     Also, the results are normalised according to the previous FFTy forward
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{
  int index = blockIdx.x*blockDim.x + threadIdx.x; //absolute thread index
  if (index<size)
  {
	re[index] = in[index].x;
    im[index] = in[index].y;
  }
}//end kernel k_sec_separate
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
__global__ void k_mixed_prod(double *a, double *b, double *out, int size)
//-------------------------------------------------------------------------------------
//
//     Multiply two arrays term by term and store the result in out; check for built-in function
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{
  int index = blockIdx.x * blockDim.x + threadIdx.x; // absolute thread index

  if (index < size)
  {
    out[index] = a[index] * b[index];
  }
}//end kernel k_mixed_prod
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
__device__ double atomic_max(double* address, double val)
{
  unsigned long long int* address_as_ull =(unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;

  while(val > __longlong_as_double(old) ) {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val));
  }

  return 0;//__longlong_as_double(old);
}//end kernel atomic max
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
__global__ void k_courant(double *u, double *v, double *w, double *dz_vec, double *result, double dx, double dy, double dt, int size, int fstart, int nx, int spy)
//-------------------------------------------------------------------------------------
//
//     Courant check kernel: very simple reduction, room for large optimisation
//     Also, optimize the way dz_vec is accessed
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{
  extern __shared__ double sdata[];

  double dz = dz_vec[fstart+blockIdx.x];

  int tid = threadIdx.x;
//  if(tid==0) printf("tid %d %lf\n",blockIdx.x,dz);
  int index = blockIdx.x * blockDim.x + threadIdx.x;

  sdata[tid] = dt * (fabs(u[index]/dx)+fabs(v[index]/dy)+fabs(w[index]/dz));
  __syncthreads();

  for (unsigned int s=blockDim.x/2; s>0; s>>=1)//bitwise shift
  {
    if (tid < s)
    {
      sdata[tid] = fmax(sdata[tid],sdata[tid + s]);
    }
    __syncthreads();
  }
  // write result for this block to global mem
  if (tid == 0) atomic_max(&result[0],sdata[0]);
  //if (tid == 0) out[blockIdx.x] = sdata[0];

}//end kernel k_courant
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
__global__ void k_norm_real(double *a,
		                    double n, int size)
//-------------------------------------------------------------------------------------
//
//     Normalize a real array output of a FFT; check the normalizing factor n
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{
  int index = blockIdx.x * blockDim.x + threadIdx.x; // absolute thread index
  if (index < size)
  {
    a[index] = a[index] * n;
  }
}//end kernel k_norm_real
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
__global__ void k_norm_cmp(cufftDoubleComplex *a,
		                    double n, int size)
//-------------------------------------------------------------------------------------
//
//     Normalize a cufftDoubleComplex array output of a FFT; check the normalizing factor n
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index < size)
  {
    a[index].x = a[index].x * n;
    a[index].y = a[index].y * n;
  }
}//end kernel k_norm_cmp
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
__global__ void k_sep_cmp(cufftDoubleComplex *in, double *out, double norm, int size)
//-------------------------------------------------------------------------------------
//
//     this kernel copies the Re of a cufftDoubleComplex array to a double array
//     and normalizes.
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index < size)
  {
    out[index] = in[index].x * norm;
  }
}
__global__ void k_sec_copy(cufftDoubleComplex *out,
		                   cufftDoubleComplex *im,
		                   int size)
//-------------------------------------------------------------------------------------
//
//     this kernel copies back together Re and Im parts of the signal in an unique array of kind R I R I R I;
//     Also, the results are halved
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{
  int index = blockIdx.x*blockDim.x + threadIdx.x;
  if (index<size)
  {
	out[index].x = out[index].x * 0.5e0;
    out[index].y = im[index].x  * 0.5e0;
  }
}//end kernel k_sec_copy
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
__global__ void k_manip(double *a, int nz, int dim)
//-------------------------------------------------------------------------------------
//
//     this kernel performs the multiplication of the last mode
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{
  int index = threadIdx.x + blockIdx.x * blockDim.x + 1; //absolute thread index
  int pos = nz * index - 1;
  if (pos < dim) //redundant
  {
    a[pos] = 2.0e0 * a[pos];
  }
}//end kernel k_manip
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
__global__ void k_manip_cmp (cufftDoubleComplex *a, int nz, int size)
//-------------------------------------------------------------------------------------
//
//     this kernel performs the multiplication of the last mode for complex phys to spec
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{
  int index = threadIdx.x + blockIdx.x * blockDim.x + 1; //absolute thread index plus one
  int pos = nz * index - 1;
  if (pos < size) //redundant check actually
  {
    a[pos].x = 0.5e0 * a[pos].x;
    a[pos].y = 0.5e0 * a[pos].y;
  }
}//end kernel k_manip_cmp
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
__global__ void k_mirr_bigtime(double *out1, double *out2,
		                       double *in1,  double *in2,
                               int nx, int nx_big, int size)
//-------------------------------------------------------------------------------------
//
//	   Mirror array in 1st direction, general form
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{

	extern __shared__ double tile[];

	int index = blockIdx.x * blockDim.x + threadIdx.x; //absolute thread index


	int p  = index / nx;    //period in first array
	int clk = index % nx; //clock in first array

	if (index < size)
	{
	  out1[clk + p * nx_big] = in1[index]; //write copy to output
	  out2[clk + p * nx_big] = in2[index]; //write copy to output

	  if (clk > 0 && clk < nx-1)
	  {
	    out1[nx_big - clk + p * nx_big] = in1[index];
	    out2[nx_big - clk + p * nx_big] = in2[index];
	  }
	}
}//end kernel k_mirr_bigtime
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
__global__ void k_mirror(double *out,
		                 double *in,
                         int nx, int nx_big)
//-------------------------------------------------------------------------------------
//
//     transform the input in an even symmetrical signal along the first direction;
//     shared memory strategy; tested up to 256, exibits good bandwidth comparable to the shared memory copy
//     and takes approximately twice the time since it processes approx twice the data;
//     can be improved by working on the second if, cause for warp divergence
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{

	extern __shared__ double tile[];
	int i = blockIdx.x; //from zero to ny*nz
	//int ind = threadIdx.x + blockDim.x * blockIdx.x;//abs thread index in the array
	int index = threadIdx.x;//index of thread in the block

	int p  = i * nx;
	int pp = i * nx_big;

	if ( index < nx)
	{
	  out[index + pp] = in [index + p]; //write copy to output
	  tile [index] = in[index + p]; //write mirror copy to input
	}

	__syncthreads();//we need to synchronize due to the fact that threads read and write different values

	if ( index  < nx_big - nx)
	{
	  out[index + pp + nx] = tile[nx - 2 - index];//allows to avoid first value
	}
}//end kernel k_mirror
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
__global__ void k_mirror_2n(double *out,
		                    double *in,
                            int nx, int nx_big)
//-------------------------------------------------------------------------------------
//
//     transform the input in a 2N simmetrical signal along the first direction;
//     shared memory strategy; tested up to 256, exibits good bandwidth comparable to the shared memory copy
//     and takes approximately twice the time since it processes approx twice the data;
//     can be improved by working on the second if, cause for warp divergence
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{

	extern __shared__ double tile[];
	int i = blockIdx.x; //from zero to ny*nz
	//int ind = threadIdx.x + blockDim.x * blockIdx.x;//abs thread index in the array
	int index = threadIdx.x;//index of thread in the block

	int p  = i * nx;
	int pp = i * nx_big;

	if ( index < nx)
	{
	  out[index + pp] = in [index + p]; //write copy to output
	  tile [index] = in[index + p]; //write mirror copy to input
	}

	__syncthreads();//we need to synchronize due to the fact that threads read and write different values

	if ( index  < nx_big - nx)
	{
	  out[index + pp + nx] = tile[nx - 1 - index];//allows to avoid first value
	}
}//end kernel k_mirror_2n
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
void __global__ k_alias_small(double *a, int ali, int nx)
//-------------------------------------------------------------------------------------
//
//     perform aliasing on double  array
//	   NOT TESTED!!!
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{
  int index = blockIdx.x*blockDim.x+threadIdx.x;   //absolute thread index
  int check = ali + index/(nx-ali)*nx + index % (nx-ali);

  //each thread accesses only the values to be modified
  a[check] = 0.0e0;

}//end kernel k_alias_small
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
void __global__ k_cmp_alias_small(cufftDoubleComplex *a, int ali, int nx)
//-------------------------------------------------------------------------------------
//
//     perform aliasing on cufftDoubleComplex velocity array
//     access only the values that need to be modified to zero; one per thread
//     this kernel expresses lower bandwidth compared to the cmp_alias and is faster for small arrays
//     at 256**3 the two kernels take the same time as tested on DAVIDE
//     the bigger the array's aliased size, the better for coalescence.
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{
  int index = blockIdx.x*blockDim.x+threadIdx.x;   //absolute thread index
  int check = ali + index/(nx-ali)*nx + index % (nx-ali);

  //each thread accesses only the values to be modified
  a[check].x = 0.0e0;

}//end kernel k_cmp_alias_small
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
void __global__ k_cmp_alias(cufftDoubleComplex *a, int al_low, int al_val, int nx, int tot_size)
//-------------------------------------------------------------------------------------
//
//     perform aliasing on cufftDoubleComplex velocity array
//     access only necessary values, requires positions for first and last aliased values
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;   //absolute thread index
  int check = al_low + index % al_val + index/al_val * nx;

  //each thread accesses only the values to be modified
  if (check < tot_size)
  {
    a[check].x = 0.0e0;
    a[check].y = 0.0e0;
  }
}//end kernel k_cmp_alias
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
void __global__ cmp_alias(cufftDoubleComplex *a, int ali, int nx)
//-------------------------------------------------------------------------------------
//
//     perform aliasing on cufftDoubleComplex velocity array
//     access all length of the array and select the indexes to reset; could be possible to access only the necessary values
//
//     Copyright Multiphase Flow Laboratory, University of Udine
//     authors - D. Di Giusto, Jan 2020
//
//-------------------------------------------------------------------------------------
{
  int index = blockIdx.x*blockDim.x+threadIdx.x;   //absolute thread index
  int check = (index - nx * (index/nx))/ali;

  //each thread has now a value for check
  if (check)
  {
    a[index].x = 0.0e0;
  }
}//end kernel cmp_alias
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
//-------------------------------------------------------------------------------------
//    This file is freely modified from the cuTranspose.
//    Copyright 2016 Ibai Gurrutxaga, Javier Muguerza, Jose L. Jodra.
//    cuTranspose is free software: it can be redistributed it and/or modified
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//-------------------------------------------------------------------------------------
__global__ void k_t102(double *out,
                       double *in,
                       int nx, int ny, int nz)
//-------------------------------------------------------------------------------------
//    Shared memory transposition from xyz to yxz for double data type;
//    Each thread handles one double at loop;
//
//-------------------------------------------------------------------------------------
{
	__shared__ double tile[TILE_DIM][TILE_DIM + 1];

	int x_in, y_in, z,
	    x_out, y_out,
	    ind_in,
	    ind_out;


	x_in = threadIdx.x + TILE_DIM * blockIdx.x;
	y_in = threadIdx.y + TILE_DIM * blockIdx.y;

	z = blockIdx.z;

	x_out = threadIdx.y + TILE_DIM * blockIdx.x;
	y_out = threadIdx.x + TILE_DIM * blockIdx.y;


	ind_in = x_in + (y_in + z * ny) * nx;
	ind_out = y_out + (x_out + z * nx) * ny;

  	for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
	{
	  if( x_in < nx && y_in + j < ny )
	  {
		tile[threadIdx.x][threadIdx.y+j] = in[ind_in+j*nx];
	  }
	}

	__syncthreads();

  	for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
	{
	  if( x_out + j < nx && y_out < ny )
	  {
		out[ind_out+j*ny] = tile[threadIdx.y+j][threadIdx.x];
	  }
	}
}//end kernel k_t102
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
__global__ void k_cmp_t102(cufftDoubleComplex *out,
                           cufftDoubleComplex *in,
                           int nx, int ny, int nz)
//-------------------------------------------------------------------------------------
//    Shared memory transposition from xyz to yxz for cufftDoubleComplex data type;
//    Each thread handles one couple of R I;
//
//-------------------------------------------------------------------------------------
{
	__shared__ cufftDoubleComplex tile[TILE_DIM][TILE_DIM + 1];

	int x_in, y_in, z,
	    x_out, y_out,
	    ind_in,
	    ind_out;


	x_in = threadIdx.x + TILE_DIM * blockIdx.x;
	y_in = threadIdx.y + TILE_DIM * blockIdx.y;

	z = blockIdx.z;

	x_out = threadIdx.y + TILE_DIM * blockIdx.x;
	y_out = threadIdx.x + TILE_DIM * blockIdx.y;


	ind_in = x_in + (y_in + z * ny) * nx;
	ind_out = y_out + (x_out + z * nx) * ny;

  	for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
	{
	  if( x_in < nx && y_in + j < ny )
	  {
		tile[threadIdx.x][threadIdx.y+j].x = in[ind_in+j*nx].x;
		tile[threadIdx.x][threadIdx.y+j].y = in[ind_in+j*nx].y;
	  }
	}

	__syncthreads();

  	for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
	{
	  if( x_out + j < nx && y_out < ny )
	  {
		out[ind_out+j*ny].x = tile[threadIdx.y+j][threadIdx.x].x;
		out[ind_out+j*ny].y = tile[threadIdx.y+j][threadIdx.x].y;
	  }
	}
}//end kernel k_cmp_t102
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
__global__ void k_cmp_t210(cufftDoubleComplex *out,
		                   cufftDoubleComplex *in,
		                   int np0, int np1, int np2)
//-------------------------------------------------------------------------------------
//    Shared memory transposition from xyz to zyx for cufftDoubleComplex data type;
//    Each thread handles one couple of R I;
//    Warning-GPU: could be optimised to make each thread handle one value only?
//-------------------------------------------------------------------------------------
{
	__shared__ cufftDoubleComplex tile[TILE_DIM][TILE_DIM + 1];

	int x_in, y, z_in,
	    x_out, z_out,
	    ind_in,
	    ind_out;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x_in = lx + TILE_DIM * bx;
	z_in = ly + TILE_DIM * by;

	y = blockIdx.z;

	x_out = ly + TILE_DIM * bx;
	z_out = lx + TILE_DIM * by;


	ind_in = x_in + (y + z_in * np1) * np0;
	ind_out = z_out + (y + x_out * np1) * np2;


  	for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
	{
  	  if( x_in < np0 && z_in + j < np2 )
	  {
			tile[lx][ly+j].x = in[ind_in+j*np0*np1].x;
			tile[lx][ly+j].y = in[ind_in+j*np0*np1].y;
	  }
	}

	__syncthreads();

  	for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
	{
	  if( z_out < np2 && x_out + j < np0 )
	  {
		out[ind_out+j*np2*np1].x = tile[ly + j][lx].x;
		out[ind_out+j*np2*np1].y = tile[ly + j][lx].y;
	  }
	}
}//end kernel k_cmp_t210
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
__global__ void k_cmp_t021(cufftDoubleComplex *out,
                           cufftDoubleComplex *in,
                           int nx, int ny, int nz)
//-------------------------------------------------------------------------------------
//    Shared memory transposition from xyz to xzy for cufftDoubleComplex data type;
//    Each thread handles one couple of R I;
//    Warning-GPU: could be optimised to make each thread handle one value only?
//
//-------------------------------------------------------------------------------------
{
	__shared__ cufftDoubleComplex tile[TILE_DIM][TILE_DIM + 1];

	int x, y, z,ind;

	x = threadIdx.x + TILE_DIM * blockIdx.x;
	y = threadIdx.y + TILE_DIM * blockIdx.y;
	z = blockIdx.z;

	int lx = threadIdx.x,
	    ly = threadIdx.y;

	ind = x + (y + z * ny) * nx;

  	for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
	{
	  if( x < nx && y < ny )
	  {
	  	tile[lx][ly+j].x = in[ind+j*nx].x;
		tile[lx][ly+j].y = in[ind+j*nx].y;
	  }
	}
	__syncthreads();

	ind = x + (z + y * nz) * nx;

	for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
	{
	  if( x < nx && y < ny)
	  {
	  	out[ind+j*nx*nz].x = tile[lx][ly+j].x;
		out[ind+j*nx*nz].y = tile[lx][ly+j].y;
	  }
	}
}//end kernel k_cmp_t021
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
