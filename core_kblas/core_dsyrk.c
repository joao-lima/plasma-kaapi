/*
 -- Innovative Computing Laboratory
 -- Electrical Engineering and Computer Science Department
 -- University of Tennessee
 -- (C) Copyright 2008-2010
 
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
 * @file cublas_core_dsyrk.c
 *
 *  PLASMA core_blas kernel for CUBLAS
 * (c) INRIA
 *
 * @version 2.4.6
 * @author Thierry Gautier
 *
 **/

#include "core_kblas.h"

void CORE_dsyrk_cublas(int uplo, int trans,
                int n, int k,
                double alpha, const double *A, int LDA,
                double beta, double *C, int LDC)
{
  cublasStatus_t status = cublasDsyrk_v2(
                                         kaapi_cuda_cublas_handle(),
                                         uplo,
                                         trans,
                                         n, k,
                                         &alpha, A, LDA,
                                         &beta, C, LDC);
  PLASMA_CUBLAS_ASSERT(status);
#if CONFIG_VERBOSE
  fprintf(stdout,"%s: a=%p c=%p n=%d k=%d\n", __FUNCTION__,
          A, C, n, k);
  fflush(stdout);
#endif
}


/***************************************************************************
 *
 **/
void CORE_dsyrk_quark_cublas(Quark *quark)
{
    PLASMA_enum uplo;
    PLASMA_enum trans;
    int n;
    int k;
    double alpha;
    double *A;
    int lda;
    double beta;
    double *C;
    int ldc;

    quark_unpack_args_10(quark, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
    CORE_dsyrk_cublas(
                     PLASMA_CUBLAS_convertToFillMode(uplo),
                     PLASMA_CUBLAS_convertToOp(trans),
                     n, k,
                     alpha, A, lda,
                     beta, C, ldc);
}
