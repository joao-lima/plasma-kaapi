
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
 * @file cublas_core_dgemm.c
 *
 *  PLASMA core_blas kernel for CUBLAS
 * (c) INRIA
 *
 * @author Joao Lima
 *
 **/

#include "core_kblas.h"
#include "magmablas.h"

extern cudaStream_t kaapi_cuda_kernel_stream(void);

int CORE_dgessm_cublas(int M, int N, int K, int IB,
                int *IPIV,
                double *L, int LDL,
                double *A, int LDA)
{
  static double zone  =  1.0;
  static double mzone = -1.0;
  static int    ione  =  1;
  
  int i, sb;
  int tmp, tmp2;
  
  int* piv;
  cublasStatus_t status;
  
  /* Check input arguments */
  if (M < 0) {
    coreblas_error(1, "Illegal value of M");
    return -1;
  }
  if (N < 0) {
    coreblas_error(2, "Illegal value of N");
    return -2;
  }
  if (K < 0) {
    coreblas_error(3, "Illegal value of K");
    return -3;
  }
  if (IB < 0) {
    coreblas_error(4, "Illegal value of IB");
    return -4;
  }
  if ((LDL < max(1,M)) && (M > 0)) {
    coreblas_error(7, "Illegal value of LDL");
    return -7;
  }
  if ((LDA < max(1,M)) && (M > 0)) {
    coreblas_error(9, "Illegal value of LDA");
    return -9;
  }
  
  /* Quick return */
  if ((M == 0) || (N == 0) || (K == 0) || (IB == 0))
    return PLASMA_SUCCESS;
  
  piv = kaapi_memory_get_host_pointer(IPIV);
  
  for(i = 0; i < K; i += IB) {
    sb = min(IB, K-i);
    /*
     * Apply interchanges to columns I*IB+1:IB*( I+1 )+1.
     */
    tmp  = i+1;
    tmp2 = i+sb;
    magmablas_dlaswp_kaapixx(kaapi_cuda_kernel_stream(), N, A, LDA, tmp, tmp2, piv, ione);
    /*
     * Compute block row of U.
     */
    status = cublasDtrsm_v2( kaapi_cuda_cublas_handle(),
                CUBLAS_SIDE_LEFT,
                CUBLAS_FILL_MODE_LOWER,
                CUBLAS_OP_N,
                CUBLAS_DIAG_UNIT,
                sb, N, &(zone),
                &L[LDL*i+i], LDL,
                &A[i], LDA );
    PLASMA_CUBLAS_ASSERT(status);
    
    if (i+sb < M) {
      /*
       * Update trailing submatrix.
       */
      status = cublasDgemm_v2( kaapi_cuda_cublas_handle(),
                  CUBLAS_OP_N,
                  CUBLAS_OP_N,
                  M-(i+sb), N, sb,
                  &(mzone), &L[LDL*i+(i+sb)], LDL,
                  &A[i], LDA,
                  &(zone), &A[i+sb], LDA );
      PLASMA_CUBLAS_ASSERT(status);
    }
  }
  return PLASMA_SUCCESS;
}

void CORE_dgessm_quark_cublas(Quark *quark)
{
  int m;
  int n;
  int k;
  int ib;
  int *IPIV;
  double *L;
  int ldl;
  double *A;
  int lda;
  
  quark_unpack_args_9(quark, m, n, k, ib, IPIV, L, ldl, A, lda);
  CORE_dgessm_cublas(m, n, k, ib, IPIV, L, ldl, A, lda);
}
