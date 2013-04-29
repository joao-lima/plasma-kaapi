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

#include "core_kblas.h"

int CORE_dssssm_cublas(int M1, int N1, int M2, int N2, int K, int IB,
                double *A1, int LDA1,
                double *A2, int LDA2,
                double *L1, int LDL1,
                double *L2, int LDL2,
                int *IPIV)
{
  static double zone  = 1.0;
  static double mzone =-1.0;
  cublasStatus_t status;
  int* piv;
  
  int i, ii, sb;
  int im, ip;
  
  /* Check input arguments */
  if (M1 < 0) {
    coreblas_error(1, "Illegal value of M1");
    return -1;
  }
  if (N1 < 0) {
    coreblas_error(2, "Illegal value of N1");
    return -2;
  }
  if (M2 < 0) {
    coreblas_error(3, "Illegal value of M2");
    return -3;
  }
  if (N2 < 0) {
    coreblas_error(4, "Illegal value of N2");
    return -4;
  }
  if (K < 0) {
    coreblas_error(5, "Illegal value of K");
    return -5;
  }
  if (IB < 0) {
    coreblas_error(6, "Illegal value of IB");
    return -6;
  }
  if (LDA1 < max(1,M1)) {
    coreblas_error(8, "Illegal value of LDA1");
    return -8;
  }
  if (LDA2 < max(1,M2)) {
    coreblas_error(10, "Illegal value of LDA2");
    return -10;
  }
  if (LDL1 < max(1,IB)) {
    coreblas_error(12, "Illegal value of LDL1");
    return -12;
  }
  if (LDL2 < max(1,M2)) {
    coreblas_error(14, "Illegal value of LDL2");
    return -14;
  }
  
  /* Quick return */
  if ((M1 == 0) || (N1 == 0) || (M2 == 0) || (N2 == 0) || (K == 0) || (IB == 0))
    return PLASMA_SUCCESS;
  
#if CONFIG_VERBOSE
  fprintf(stdout, "%s: A1=%p A2=%p L1=%p L2=%p M1=%d N1=%d M2=%d N2=%d K=%d IB=%d LDA1=%d LDA2=%d LDL1=%d LDL2=%d\n",
          __FUNCTION__, A1, A2, L1, L2, M1, N1, M2, N2, K, IB, LDA1, LDA2, LDL1, LDL2);
  fflush(stdout);
#endif
  
  ip = 0;
  
  piv = kaapi_memory_get_host_pointer(IPIV);
  
  for(ii = 0; ii < K; ii += IB) {
    sb = min(K-ii, IB);
    
    for(i = 0; i < sb; i++) {
      im = piv[ip]-1;
      
      if (im != (ii+i)) {
        im = im - M1;
        status = cublasDswap_v2
        (
           kaapi_cuda_cublas_handle(),
           N1,
             &A1[ii+i], LDA1,
             &A2[im], LDA2
         );
        PLASMA_CUBLAS_ASSERT(status);
      }
      ip = ip + 1;
    }
    
    status = cublasDtrsm_v2
    (
       kaapi_cuda_cublas_handle(),
       PLASMA_CUBLAS_convertToSideMode(PlasmaLeft),
       PLASMA_CUBLAS_convertToFillMode(PlasmaLower),
       PLASMA_CUBLAS_convertToOp(PlasmaNoTrans),
       PLASMA_CUBLAS_convertToDiagType(PlasmaUnit),
       sb, N1, &zone,
       &L1[LDL1*ii], LDL1,
       &A1[ii], LDA1
     );
    PLASMA_CUBLAS_ASSERT(status);
    
    status = cublasDgemm_v2
    (
       kaapi_cuda_cublas_handle(),
       PLASMA_CUBLAS_convertToOp(PlasmaNoTrans),
       PLASMA_CUBLAS_convertToOp(PlasmaNoTrans),
       M2, N2, sb,
       &mzone,
       &L2[LDL2*ii], LDL2,
       &A1[ii], LDA1,
       &zone, A2, LDA2
     );
    PLASMA_CUBLAS_ASSERT(status);
  }
  return PLASMA_SUCCESS;
}


void CORE_dssssm_quark_cublas(Quark *quark)
{
  int m1;
  int n1;
  int m2;
  int n2;
  int k;
  int ib;
  double *A1;
  int lda1;
  double *A2;
  int lda2;
  double *L1;
  int ldl1;
  double *L2;
  int ldl2;
  int *IPIV;
  
  quark_unpack_args_15(quark, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, L1, ldl1, L2, ldl2, IPIV);
  CORE_dssssm_cublas(m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, L1, ldl1, L2, ldl2, IPIV);
}

int CORE_dssssm_cublas_v2(int M1, int N1, int M2, int N2, int K, int IB,
                       double *A1, int LDA1,
                       double *A2, int LDA2,
                       double *L1, int LDL1,
                       double *L2, int LDL2,
                       int *piv)
{
  static double zone  = 1.0;
  static double mzone =-1.0;
  cublasStatus_t status;
  
  int i, ii, sb;
  int im, ip;
  
  /* Check input arguments */
  if (M1 < 0) {
    coreblas_error(1, "Illegal value of M1");
    return -1;
  }
  if (N1 < 0) {
    coreblas_error(2, "Illegal value of N1");
    return -2;
  }
  if (M2 < 0) {
    coreblas_error(3, "Illegal value of M2");
    return -3;
  }
  if (N2 < 0) {
    coreblas_error(4, "Illegal value of N2");
    return -4;
  }
  if (K < 0) {
    coreblas_error(5, "Illegal value of K");
    return -5;
  }
  if (IB < 0) {
    coreblas_error(6, "Illegal value of IB");
    return -6;
  }
  if (LDA1 < max(1,M1)) {
    coreblas_error(8, "Illegal value of LDA1");
    return -8;
  }
  if (LDA2 < max(1,M2)) {
    coreblas_error(10, "Illegal value of LDA2");
    return -10;
  }
  if (LDL1 < max(1,IB)) {
    coreblas_error(12, "Illegal value of LDL1");
    return -12;
  }
  if (LDL2 < max(1,M2)) {
    coreblas_error(14, "Illegal value of LDL2");
    return -14;
  }
  
  /* Quick return */
  if ((M1 == 0) || (N1 == 0) || (M2 == 0) || (N2 == 0) || (K == 0) || (IB == 0))
    return PLASMA_SUCCESS;
  
#if CONFIG_VERBOSE
  fprintf(stdout, "%s: A1=%p A2=%p L1=%p L2=%p M1=%d N1=%d M2=%d N2=%d K=%d IB=%d LDA1=%d LDA2=%d LDL1=%d LDL2=%d\n",
          __FUNCTION__, A1, A2, L1, L2, M1, N1, M2, N2, K, IB, LDA1, LDA2, LDL1, LDL2);
  fflush(stdout);
#endif
  
  ip = 0;
  
  for(ii = 0; ii < K; ii += IB) {
    sb = min(K-ii, IB);
    
    for(i = 0; i < sb; i++) {
      im = piv[ip]-1;
      
      if (im != (ii+i)) {
        im = im - M1;
        status = cublasDswap_v2
        (
         kaapi_cuda_cublas_handle(),
         N1,
         &A1[ii+i], LDA1,
         //           &a1[ii+i], 1,
         //           &a2[im], 1
         &A2[im], LDA2
         );
        PLASMA_CUBLAS_ASSERT(status);
      }
      ip = ip + 1;
    }
    
    status = cublasDtrsm_v2
    (
     kaapi_cuda_cublas_handle(),
     PLASMA_CUBLAS_convertToSideMode(PlasmaLeft),
     PLASMA_CUBLAS_convertToFillMode(PlasmaLower),
     PLASMA_CUBLAS_convertToOp(PlasmaNoTrans),
     PLASMA_CUBLAS_convertToDiagType(PlasmaUnit),
     sb, N1, &zone,
     &L1[LDL1*ii], LDL1,
     &A1[ii], LDA1
     );
    PLASMA_CUBLAS_ASSERT(status);
    
    status = cublasDgemm_v2
    (
     kaapi_cuda_cublas_handle(),
     PLASMA_CUBLAS_convertToOp(PlasmaNoTrans),
     PLASMA_CUBLAS_convertToOp(PlasmaNoTrans),
     M2, N2, sb,
     &mzone,
     &L2[LDL2*ii], LDL2,
     &A1[ii], LDA1,
     &zone, A2, LDA2
     );
    PLASMA_CUBLAS_ASSERT(status);
  }
  return PLASMA_SUCCESS;
}