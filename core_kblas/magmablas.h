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
 * @file magmablas.h
 *
 *  PLASMA core_blas kernel for CUBLAS
 * (c) INRIA
 *
 * @author Joao Lima
 *
 * Based on MAGMA (version 1.1.0) http://icl.cs.utk.edu/magma
 *
 **/

#ifndef _MAGMABLAS_XKAAPI_H_
#define _MAGMABLAS_XKAAPI_H_

#include <cuda_runtime_api.h>

#if defined(__cplusplus)
extern "C" {
#endif

void magmablas_dlaswp_kaapixx(cudaStream_t stream, int n, double *dAT, int lda,
                  int i1, int i2, int *ipiv, int inci);

void magmablas_dtranspose_kaapixx(cudaStream_t stream, double *odata, int ldo,
                     double *idata, int ldi, 
                     int m, int n );

void magmablas_dswapblk_kaapixx(cudaStream_t stream, int n,
                    double *dA1T, int lda1,
                    double *dA2T, int lda2,
                    int i1, int i2, int *ipiv, int inci, int offset);
  
void magmablas_dlacpy_kaapixx(cudaStream_t stream, char uplo, int m, int n,
                                double *a, int lda,
                              double *b, int ldb );
  
void magmablas_slaswp_kaapixx(cudaStream_t stream, int n, float *dAT, int lda,
                              int i1, int i2, int *ipiv, int inci);

void magmablas_stranspose_kaapixx(cudaStream_t stream, float *odata, int ldo,
                                  float *idata, int ldi,
                                  int m, int n );

void magmablas_sswapblk_kaapixx(cudaStream_t stream, int n,
                                float *dA1T, int lda1,
                                float *dA2T, int lda2,
                                int i1, int i2, int *ipiv, int inci, int offset);
  
#if defined(__cplusplus)
}
#endif

#endif /* _MAGMABLAS_XKAAPI_H_ */