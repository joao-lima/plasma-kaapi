/*
 -- Innovative Computing Laboratory
 -- Electrical Engineering and Computer Science Department
 -- University of Tennessee
 -- (C) Copyright 2009
 
 Redistribution  and  use  in  source and binary forms, with or without
 modification,  are  permitted  provided  that the following conditions
 are met:
 
 * Redistributions  of  source  code  must  retain  the above copyright
 notice,  this  list  of  conditions  and  the  following  disclaimer.
 * Redistributions  in  binary  form must reproduce the above copyright
 notice,  this list of conditions and the following disclaimer in the
 documentation  and/or other materials provided with the distribution.
 * Neither  the  name of the University of Tennessee, Knoxville nor the
 names of its contributors may be used to endorse or promote products
 derived from this software without specific prior written permission.
 
 THIS  SOFTWARE  IS  PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 ``AS IS''  AND  ANY  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A  PARTICULAR  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL,  EXEMPLARY,  OR  CONSEQUENTIAL  DAMAGES  (INCLUDING,  BUT NOT
 LIMITED  TO,  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA,  OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 THEORY  OF  LIABILITY,  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF  THIS  SOFTWARE,  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/**
 *
 * @file magmablas.cu
 *
 *  PLASMA core_blas kernel for CUBLAS
 * (c) INRIA
 *
 * @author Joao Lima
 *
 * Based on MAGMA (version 1.1.0) http://icl.cs.utk.edu/magma
 *
 **/

#ifndef BLOCK_SIZE
#define BLOCK_SIZE  64
#endif

#define PRECISION_d

#include <cuda.h>
#include "magmablas.h"

typedef struct {
  double *A;
  int n, lda, j0, npivots;
  short ipiv[BLOCK_SIZE];
} kaapixx_dlaswp_params_t2;

typedef struct {
  float *A;
  int n, lda, j0, npivots;
  short ipiv[BLOCK_SIZE];
} kaapixx_slaswp_params_t2;

__global__ void mydlaswp2_kaapixx( kaapixx_dlaswp_params_t2 params )
{
  unsigned int tid = threadIdx.x + __mul24(blockDim.x, blockIdx.x);
  if( tid < params.n )
  {
    int lda = params.lda;
    double *A = params.A + tid + lda * params.j0;
    
    for( int i = 0; i < params.npivots; i++ )
    {
      int j = params.ipiv[i];
      double *p1 = A + i*lda;
      double *p2 = A + j*lda;
      double temp = *p1;
      *p1 = *p2;
      *p2 = temp;
    }
  }
}

__global__ void myslaswp2_kaapixx( kaapixx_slaswp_params_t2 params )
{
  unsigned int tid = threadIdx.x + __mul24(blockDim.x, blockIdx.x);
  if( tid < params.n )
  {
    int lda = params.lda;
    float *A = params.A + tid + lda * params.j0;
    
    for( int i = 0; i < params.npivots; i++ )
    {
      int j = params.ipiv[i];
      float *p1 = A + i*lda;
      float *p2 = A + j*lda;
      float temp = *p1;
      *p1 = *p2;
      *p2 = temp;
    }
  }
}

#define	DSIZE_1SHARED	32

__global__ void dtranspose_32_kaapixx( double *B, int ldb, double *A, int lda )
{
  __shared__ double a[32][DSIZE_1SHARED+1];
  
  int inx = threadIdx.x;
  int iny = threadIdx.y;
  int ibx = blockIdx.x*32;
  int iby = blockIdx.y*32;
  
  A += ibx + inx + __mul24( iby + iny, lda );
  B += iby + inx + __mul24( ibx + iny, ldb );
  
  a[iny+0][inx] = A[0*lda];
  a[iny+8][inx] = A[8*lda];
  a[iny+16][inx] = A[16*lda];
  a[iny+24][inx] = A[24*lda];
  
  __syncthreads();
  
#if defined(PRECISION_s) || defined(PRECISION_d) || defined(PRECISION_c)
  B[0*ldb] = a[inx][iny+0];
  B[8*ldb] = a[inx][iny+8];
  B[16*ldb] = a[inx][iny+16];
  B[24*ldb] = a[inx][iny+24];
#else /* defined(PRECISION_z) */
  B[0*ldb]    = a[inx][iny+0];
  B[8*ldb]    = a[inx][iny+8];
  B[0*ldb+16] = a[inx+16][iny+0];
  B[8*ldb+16] = a[inx+16][iny+8];
  
  __syncthreads();
  A += DSIZE_1SHARED;
  B += __mul24( 16, ldb);
  
  a[iny+0][inx] = A[0*lda];
  a[iny+8][inx] = A[8*lda];
  a[iny+16][inx] = A[16*lda];
  a[iny+24][inx] = A[24*lda];
  
  __syncthreads();
  
  B[0*ldb] = a[inx][iny+0];
  B[8*ldb] = a[inx][iny+8];
  B[0*ldb+16] = a[inx+16][iny+0];
  B[8*ldb+16] = a[inx+16][iny+8];
#endif
}

typedef struct {
  double *A1;
  double *A2;
  int n, lda1, lda2, npivots;
  short ipiv[BLOCK_SIZE];
} kaapixx_magmagpu_dswapblk_params_t;

typedef struct {
  float *A1;
  float *A2;
  int n, lda1, lda2, npivots;
  short ipiv[BLOCK_SIZE];
} kaapixx_magmagpu_sswapblk_params_t;

__global__ void magmagpu_dswapblkrm_kaapixx( kaapixx_magmagpu_dswapblk_params_t params )
{
  unsigned int y = threadIdx.x + blockDim.x*blockIdx.x;
  if( y < params.n )
  {
    double *A1 = params.A1 + y - params.lda1;
    double *A2 = params.A2 + y;
    
    for( int i = 0; i < params.npivots; i++ )
    {
      A1 += params.lda1;
      if ( params.ipiv[i] == -1 )
        continue;
      double tmp1  = *A1;
      double *tmp2 = A2 + params.ipiv[i]*params.lda2;
      *A1   = *tmp2;
      *tmp2 = tmp1;
    }
  }
}

__global__ void magmagpu_sswapblkrm_kaapixx( kaapixx_magmagpu_sswapblk_params_t params )
{
  unsigned int y = threadIdx.x + blockDim.x*blockIdx.x;
  if( y < params.n )
  {
    float *A1 = params.A1 + y - params.lda1;
    float *A2 = params.A2 + y;
    
    for( int i = 0; i < params.npivots; i++ )
    {
      A1 += params.lda1;
      if ( params.ipiv[i] == -1 )
        continue;
      float tmp1  = *A1;
      float *tmp2 = A2 + params.ipiv[i]*params.lda2;
      *A1   = *tmp2;
      *tmp2 = tmp1;
    }
  }
}

__global__ void magmagpu_dswapblkcm_kaapixx( kaapixx_magmagpu_dswapblk_params_t params )
{
  unsigned int y = threadIdx.x + blockDim.x*blockIdx.x;
  unsigned int offset1 = __mul24( y, params.lda1);
  unsigned int offset2 = __mul24( y, params.lda2);
  if( y < params.n )
  {
    double *A1 = params.A1 + offset1 - 1;
    double *A2 = params.A2 + offset2;
    
    for( int i = 0; i < params.npivots; i++ )
    {
      A1++;
      if ( params.ipiv[i] == -1 )
        continue;
      double tmp1  = *A1;
      double *tmp2 = A2 + params.ipiv[i];
      *A1   = *tmp2;
      *tmp2 = tmp1;
    }
  }
  __syncthreads();
}

__global__ void magmagpu_sswapblkcm_kaapixx( kaapixx_magmagpu_sswapblk_params_t params )
{
  unsigned int y = threadIdx.x + blockDim.x*blockIdx.x;
  unsigned int offset1 = __mul24( y, params.lda1);
  unsigned int offset2 = __mul24( y, params.lda2);
  if( y < params.n )
  {
    float *A1 = params.A1 + offset1 - 1;
    float *A2 = params.A2 + offset2;
    
    for( int i = 0; i < params.npivots; i++ )
    {
      A1++;
      if ( params.ipiv[i] == -1 )
        continue;
      float tmp1  = *A1;
      float *tmp2 = A2 + params.ipiv[i];
      *A1   = *tmp2;
      *tmp2 = tmp1;
    }
  }
  __syncthreads();
}

__global__ void magmagpu_dlacpy_kaapixx_kernel_LU(int m, int n,
                                                  double *a, int lda,
                                                  double *b, int ldb )
{
  const int i = blockIdx.x*64 + threadIdx.x;
  int j;
  
  if( i >= n )
    return;
  
  for (j = 0; j < m; j++) {
    b[i + j * ldb] = a[i + j * lda];
  }
}

__global__ void magmagpu_dlacpy_kaapixx_kernel_U(int m, int n,
                                                 double *a, int lda,
                                                 double *b, int ldb )
{
  const int i = blockIdx.x*64 + threadIdx.x;
  int j;
  
  if( i >= m )
    return;
  
  const int i__ = min(i, n);
  
  for (j = 0; j < i__; j++) {
    b[j + i * ldb] = a[j + i * lda];
  }
}

__global__ void magmagpu_dlacpy_kaapixx_kernel_L(int m, int n,
                                                 double *a, int lda,
                                                 double *b, int ldb )
{
  const int i = blockIdx.x*64 + threadIdx.x;
  int j;
  
  if( i >= m )
    return;
  
  for (j = i; j < n; j++) {
    b[j + i * ldb] = a[j + i * lda];
  }
}

#if defined(__cplusplus)
extern "C" {
#endif
  
#include <stdio.h>
  
void dlaswp3_kaapixx( cudaStream_t stream, kaapixx_dlaswp_params_t2 &params )
{
  int blocksize = 64;
  dim3 blocks = (params.n+blocksize-1) / blocksize;
  mydlaswp2_kaapixx<<< blocks, blocksize, 0, stream >>>( params );
#if defined(CONFIG_DEBUG)
  cudaError_t res = cudaGetLastError();
  if(res != cudaSuccess){
    fprintf(stdout, "ERROR in %s: %s\n", __FUNCTION__, cudaGetErrorString(res));
    fflush(stdout);
  }
#endif
}

void slaswp3_kaapixx( cudaStream_t stream, kaapixx_slaswp_params_t2 &params )
{
  int blocksize = 64;
  dim3 blocks = (params.n+blocksize-1) / blocksize;
  myslaswp2_kaapixx<<< blocks, blocksize, 0, stream >>>( params );
#if defined(CONFIG_DEBUG)
  cudaError_t res = cudaGetLastError();
  if(res != cudaSuccess){
    fprintf(stdout, "ERROR in %s: %s\n", __FUNCTION__, cudaGetErrorString(res));
    fflush(stdout);
  }
#endif
}

void
magmablas_dlaswp_kaapixx( cudaStream_t stream, int n, double *dAT, int lda,
                 int i1, int i2, int *ipiv, int inci )
{
  int k;
  
  for( k=(i1-1); k<i2; k+=BLOCK_SIZE )
  {
    int sb = min(BLOCK_SIZE, i2-k);
    //dlaswp_params_t params = { dAT, lda, lda, ind + k };
    kaapixx_dlaswp_params_t2 params = { dAT+k*lda, n, lda, 0, sb };
    for( int j = 0; j < sb; j++ )
    {
      params.ipiv[j] = ipiv[(k+j)*inci] - k - 1;
    }
    dlaswp3_kaapixx( stream, params );
  }
}

void
magmablas_slaswp_kaapixx( cudaStream_t stream, int n, float *dAT, int lda,
                         int i1, int i2, int *ipiv, int inci )
{
  int k;
  
  for( k=(i1-1); k<i2; k+=BLOCK_SIZE )
  {
    int sb = min(BLOCK_SIZE, i2-k);
    //dlaswp_params_t params = { dAT, lda, lda, ind + k };
    kaapixx_slaswp_params_t2 params = { dAT+k*lda, n, lda, 0, sb };
    for( int j = 0; j < sb; j++ )
    {
      params.ipiv[j] = ipiv[(k+j)*inci] - k - 1;
    }
    slaswp3_kaapixx( stream, params );
  }
}

void magmablas_dtranspose_kaapixx(cudaStream_t stream, double *odata, int ldo,
                          double *idata, int ldi,
                          int m, int n )
{
  dim3 threads( DSIZE_1SHARED, 8, 1 );
  dim3 grid( m/32, n/32, 1 );
  dtranspose_32_kaapixx<<<grid, threads, 0, stream>>>( odata, ldo, idata, ldi );
}

void magmablas_dswapblk_kaapixx(cudaStream_t stream, int n,
                        double *dA1T, int lda1,
                        double *dA2T, int lda2,
                        int i1, int i2, int *ipiv, int inci, int offset)
{
  int  blocksize = 64;
  dim3 blocks( (n+blocksize-1) / blocksize, 1, 1);
  int  k, im;
  
  for( k=(i1-1); k<i2; k+=BLOCK_SIZE )
  {
    int sb = min(BLOCK_SIZE, i2-k);
    kaapixx_magmagpu_dswapblk_params_t params = { dA1T+k*lda1, dA2T, n, lda1, lda2, sb };
    for( int j = 0; j < sb; j++ )
    {
      im = ipiv[(k+j)*inci] - 1;
      if ( (k+j) == im)
        params.ipiv[j] = -1;
      else
        params.ipiv[j] = im - offset;
    }
    magmagpu_dswapblkrm_kaapixx<<< blocks, blocksize, 0, stream >>>( params );
  }
}

void magmablas_sswapblk_kaapixx(cudaStream_t stream, int n,
                                float *dA1T, int lda1,
                                float *dA2T, int lda2,
                                int i1, int i2, int *ipiv, int inci, int offset)
{
  int  blocksize = 64;
  dim3 blocks( (n+blocksize-1) / blocksize, 1, 1);
  int  k, im;
  
  for( k=(i1-1); k<i2; k+=BLOCK_SIZE )
  {
    int sb = min(BLOCK_SIZE, i2-k);
    kaapixx_magmagpu_sswapblk_params_t params = { dA1T+k*lda1, dA2T, n, lda1, lda2, sb };
    for( int j = 0; j < sb; j++ )
    {
      im = ipiv[(k+j)*inci] - 1;
      if ( (k+j) == im)
        params.ipiv[j] = -1;
      else
        params.ipiv[j] = im - offset;
    }
    magmagpu_sswapblkrm_kaapixx<<< blocks, blocksize, 0, stream >>>( params );
  }
}

void magmablas_dlacpy_kaapixx(cudaStream_t stream, char uplo, int m, int n,
                              double *a, int lda,
                              double *b, int ldb )
{
  dim3 threads( 64 );
  dim3 grid( m/64 + (m%64 != 0) );
  
  if ( m == 0 || n == 0 )
    return;

  if ( (uplo == 'U') || (uplo == 'u') ) {
    magmagpu_dlacpy_kaapixx_kernel_U<<< grid, threads, 0, stream >>> ( m, n, a, lda, b, ldb );
  } else if ( (uplo == 'L') || (uplo == 'l') ) {
    magmagpu_dlacpy_kaapixx_kernel_L<<< grid, threads, 0, stream >>> ( m, n, a, lda, b, ldb );
  } else {
    magmagpu_dlacpy_kaapixx_kernel_LU<<< grid, threads, 0, stream >>> ( m, n, a, lda, b, ldb );
  }
}
  
  
#if defined(__cplusplus)
}
#endif
