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

extern cudaStream_t kaapi_cuda_kernel_stream(void);

/***************************************************************************//**
*
* @ingroup CORE_double
*
*  CORE_dtsrfb applies a complex block reflector H or its transpose H' to a
*  complex rectangular matrix formed by coupling two tiles A1 and A2,
*  from either the left or the right.
*
*******************************************************************************
*
* @param[in] side
*         @arg PlasmaLeft  : apply Q or Q**T from the Left;
*         @arg PlasmaRight : apply Q or Q**T from the Right.
*
* @param[in] trans
*         @arg PlasmaNoTrans   : No transpose, apply Q;
*         @arg PlasmaTrans : ConjTranspose, apply Q**T.
*
* @param[in] direct
*         Indicates how H is formed from a product of elementary
*         reflectors
*         @arg PlasmaForward  : H = H(1) H(2) . . . H(k) (Forward)
*         @arg PlasmaBackward : H = H(k) . . . H(2) H(1) (Backward)
*
* @param[in] storev
*         Indicates how the vectors which define the elementary
*         reflectors are stored:
*         @arg PlasmaColumnwise
*         @arg PlasmaRowwise
*
* @param[in] M1
*         The number of columns of the tile A1. M1 >= 0.
*
* @param[in] N1
*         The number of rows of the tile A1. N1 >= 0.
*
* @param[in] M2
*         The number of columns of the tile A2. M2 >= 0.
*
* @param[in] N2
*         The number of rows of the tile A2. N2 >= 0.
*
* @param[in] K
*          The order of the matrix T (= the number of elementary
*          reflectors whose product defines the block reflector).
*
* @param[in,out] A1
*         On entry, the M1-by-N1 tile A1.
*         On exit, A1 is overwritten by the application of Q.
*
* @param[in] LDA1
*         The leading dimension of the array A1. LDA1 >= max(1,N1).
*
* @param[in,out] A2
*         On entry, the M2-by-N2 tile A2.
*         On exit, A2 is overwritten by the application of Q.
*
* @param[in] LDA2
*         The leading dimension of the tile A2. LDA2 >= max(1,N2).
*
* @param[in] V
*         (LDV,K) if STOREV = 'C'
*         (LDV,M2) if STOREV = 'R' and SIDE = 'L'
*         (LDV,N2) if STOREV = 'R' and SIDE = 'R'
*         The matrix V.
*
* @param[in] LDV
*         The leading dimension of the array V.
*         If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M2);
*         if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N2);
*         if STOREV = 'R', LDV >= K.
*
* @param[out] T
*         The triangular K-by-K matrix T in the representation of the
*         block reflector.
*         T is upper triangular by block (economic storage);
*         The rest of the array is not referenced.
*
* @param[in] LDT
*         The leading dimension of the array T. LDT >= K.
*
* @param[in,out] WORK
*
* @param[in] LDWORK
*         The dimension of the array WORK.
*
*******************************************************************************
*
* @return
*          \retval PLASMA_SUCCESS successful exit
*          \retval <0 if -i, the i-th argument had an illegal value
*
******************************************************************************/

int CORE_dtsrfb_cublas(int side, int trans, int direct, int storev,
                int M1, int N1, int M2, int N2, int K,
                double *A1, int LDA1,
                double *A2, int LDA2,
                double *V, int LDV,
                double *T, int LDT,
                double *WORK, int LDWORK)
{
  static double zone  =  1.0;
  static double mzone = -1.0;
  
  int j;
  
#if defined(CONFIG_VERBOSE)
  fprintf(stdout, "%s: M1=%d N1=%d M2=%d N2=%d K=%d A1=%p LDA1=%d A2=%p LDA2=%d V=%p LDV=%d T=%p LDT=%d WORK=%p LDWORK=%d\n", __FUNCTION__,
          M1, N1, M2, N2, K, A1, LDA1, A2, LDA2, V, LDV, T, LDT, WORK, LDWORK);
  fflush(stdout);
#endif
  /* Check input arguments */
  if (M1 < 0) {
    coreblas_error(5, "Illegal value of M1");
    return -5;
  }
  if (N1 < 0) {
    coreblas_error(6, "Illegal value of N1");
    return -6;
  }
  if ( (M2 < 0) || 
      ( (M2 != M1) && (side == PlasmaRight) ) ){
    coreblas_error(7, "Illegal value of M2");
    return -7;
  }
  if ( (N2 < 0) || 
      ( (N2 != N1) && (side == PlasmaLeft) ) ){
    coreblas_error(8, "Illegal value of N2");
    return -8;
  }
  if (K < 0) {
    coreblas_error(9, "Illegal value of K");
    return -9;
  }
  
  /* Quick return */
  if ((M1 == 0) || (N1 == 0) || (M2 == 0) || (N2 == 0) || (K == 0))
    return PLASMA_SUCCESS;
  
  cublasStatus_t status;

  if (storev == PlasmaColumnwise) {
    if (direct == PlasmaForward) {
      if (side == PlasmaLeft) {
        /*
         * B = A1 + V' * A2
         */
        magmablas_dlacpy_kaapixx(kaapi_cuda_kernel_stream(),
                            'A',
//                            lapack_const(PlasmaUpperLower),
                            K, N1,
                            A1, LDA1, WORK, LDWORK);

        status = cublasDgemm_v2( kaapi_cuda_cublas_handle(),
                                 PLASMA_CUBLAS_convertToOp(PlasmaTrans),
                                 PLASMA_CUBLAS_convertToOp(PlasmaNoTrans),
                                 K, N2, M2,
                                 &zone, V, LDV,
                                 A2, LDA2,
                                 &zone, WORK, LDWORK);
        PLASMA_CUBLAS_ASSERT(status);
        
        /*
         * A2 = A2 - V*T*B -> B = T*B, A2 = A2 - V*B
         */
        status = cublasDtrmm_v2( kaapi_cuda_cublas_handle(),
                                PLASMA_CUBLAS_convertToSideMode(PlasmaLeft),
                                PLASMA_CUBLAS_convertToFillMode(PlasmaUpper),
                                PLASMA_CUBLAS_convertToOp(trans),
                                PLASMA_CUBLAS_convertToDiagType(PlasmaNonUnit),
                                K, N2, &zone, T, LDT, WORK, LDWORK, WORK, LDWORK
                                );
        PLASMA_CUBLAS_ASSERT(status);
        
        status = cublasDgemm_v2( kaapi_cuda_cublas_handle(),
                                PLASMA_CUBLAS_convertToOp(PlasmaNoTrans),
                                PLASMA_CUBLAS_convertToOp(PlasmaNoTrans),
                                M2, N2, K,
                                &mzone, V, LDV,
                                WORK, LDWORK,
                                &zone, A2, LDA2);
        PLASMA_CUBLAS_ASSERT(status);
        
        /*
         * A1 = A1 - B
         */
        for(j = 0; j < N1; j++) {
          status = cublasDaxpy_v2( kaapi_cuda_cublas_handle(),
                      K, &mzone,
                      &WORK[LDWORK*j], 1,
                      &A1[LDA1*j], 1);
          PLASMA_CUBLAS_ASSERT(status);          
        }
      }
      /*
       * Columnwise / Forward / Right
       */
      else {
        /*
         * B = A1 + A2 * V
         */
        magmablas_dlacpy_kaapixx(kaapi_cuda_kernel_stream(),
                                 'A',
//                                 lapack_const(PlasmaUpperLower),
                                 M1, K,
                                 A1, LDA1, WORK, LDWORK);

        status = cublasDgemm_v2( kaapi_cuda_cublas_handle(),
                                PLASMA_CUBLAS_convertToOp(PlasmaNoTrans),
                                PLASMA_CUBLAS_convertToOp(PlasmaNoTrans),
                                M2, K, N2,
                                &zone, A2, LDA2,
                                V, LDV,
                                &zone, WORK, LDWORK);
        PLASMA_CUBLAS_ASSERT(status);
        
        /*
         * A2 = A2 - B*T*V' -> B = B*T, A2 = A2 - B*V'
         */
        status = cublasDtrmm_v2( kaapi_cuda_cublas_handle(),
                                PLASMA_CUBLAS_convertToSideMode(PlasmaRight),
                                PLASMA_CUBLAS_convertToFillMode(PlasmaUpper),
                                PLASMA_CUBLAS_convertToOp(trans),
                                PLASMA_CUBLAS_convertToDiagType(PlasmaNonUnit),
                                M1, K, &zone, T, LDT, WORK, LDWORK, WORK, LDWORK
                                );
        PLASMA_CUBLAS_ASSERT(status);
        
        status = cublasDgemm_v2( kaapi_cuda_cublas_handle(),
                                PLASMA_CUBLAS_convertToOp(PlasmaNoTrans),
                                PLASMA_CUBLAS_convertToOp(PlasmaTrans),
                                M2, N2, K,
                                &mzone, WORK, LDWORK,
                                V, LDV,
                                &zone, A2, LDA2);
        PLASMA_CUBLAS_ASSERT(status);

        /*
         * A1 = A1 - B
         */
        for(j = 0; j < K; j++) {
          status = cublasDaxpy_v2( kaapi_cuda_cublas_handle(),
                                  M1, &mzone,
                                  &WORK[LDWORK*j], 1,
                                  &A1[LDA1*j], 1);
          PLASMA_CUBLAS_ASSERT(status);          
        }
      }
    }
    else {
      coreblas_error(3, "Not implemented (ColMajor / Backward / Left or Right)");
      return PLASMA_ERR_NOT_SUPPORTED;
    }
  }
  else {
    /* TODO */
    if (direct == PlasmaForward) {
      /*
       * Rowwise / Forward / Left
       */
      if (side == PlasmaLeft) {
        /*
         * B = A1 + V * A2
         */
        magmablas_dlacpy_kaapixx(kaapi_cuda_kernel_stream(),
                                 'A',
//                                 lapack_const(PlasmaUpperLower),
                            K, N1,
                            A1, LDA1, WORK, LDWORK);
        
        status = cublasDgemm_v2( kaapi_cuda_cublas_handle(),
                                PLASMA_CUBLAS_convertToOp(PlasmaNoTrans),
                                PLASMA_CUBLAS_convertToOp(PlasmaNoTrans),
                                K, N2, M2,
                                &zone, V, LDV,
                                A2, LDA2,
                                &zone, WORK, LDWORK);
        PLASMA_CUBLAS_ASSERT(status);

        /*
         * A2 = A2 - V'*T*B -> B = T*B, A2 = A2 - V'*B
         */
        status = cublasDtrmm_v2( kaapi_cuda_cublas_handle(),
                                PLASMA_CUBLAS_convertToSideMode(PlasmaLeft),
                                PLASMA_CUBLAS_convertToFillMode(PlasmaUpper),
                                PLASMA_CUBLAS_convertToOp(trans),
                                PLASMA_CUBLAS_convertToDiagType(PlasmaNonUnit),
                                K, N2, &zone, T, LDT, WORK, LDWORK, WORK, LDWORK
                                );
        PLASMA_CUBLAS_ASSERT(status);

        status = cublasDgemm_v2( kaapi_cuda_cublas_handle(),
                                PLASMA_CUBLAS_convertToOp(PlasmaTrans),
                                PLASMA_CUBLAS_convertToOp(PlasmaNoTrans),
                                M2, N2, K,
                                &mzone, V, LDV,
                                WORK, LDWORK,
                                &zone, A2, LDA2
                                );
        PLASMA_CUBLAS_ASSERT(status);

        /*
         * A1 = A1 - B
         */
        for(j=0; j<N1; j++) {
          status = cublasDaxpy_v2( kaapi_cuda_cublas_handle(),
                                  K, &mzone,
                                  &WORK[LDWORK*j], 1,
                                  &A1[LDA1*j], 1);
          PLASMA_CUBLAS_ASSERT(status);
        }
      }
      /*
       * Rowwise / Forward / Right
       */
      else {
        /*
         * B = A1 + A2 * V'
         */
        magmablas_dlacpy_kaapixx(kaapi_cuda_kernel_stream(),
                                 'A',
//                                 lapack_const(PlasmaUpperLower),
                            M1, K,
                            A1, LDA1, WORK, LDWORK);
        
        status = cublasDgemm_v2( kaapi_cuda_cublas_handle(),
                                PLASMA_CUBLAS_convertToOp(PlasmaNoTrans),
                                PLASMA_CUBLAS_convertToOp(PlasmaTrans),
                                M2, K, N2,
                                &zone, A2, LDA2,
                                V, LDV,
                                &zone, WORK, LDWORK);
        PLASMA_CUBLAS_ASSERT(status);
        
        /*
         * A2 = A2 - B*T*V -> B = B*T, A2 = A2 - B*V'
         */
        status = cublasDtrmm_v2( kaapi_cuda_cublas_handle(),
                                PLASMA_CUBLAS_convertToSideMode(PlasmaRight),
                                PLASMA_CUBLAS_convertToFillMode(PlasmaUpper),
                                PLASMA_CUBLAS_convertToOp(trans),
                                PLASMA_CUBLAS_convertToDiagType(PlasmaNonUnit),
                                M1, K, &zone, T, LDT, WORK, LDWORK, WORK, LDWORK
                                );
        PLASMA_CUBLAS_ASSERT(status);

        status = cublasDgemm_v2( kaapi_cuda_cublas_handle(),
                                PLASMA_CUBLAS_convertToOp(PlasmaNoTrans),
                                PLASMA_CUBLAS_convertToOp(PlasmaNoTrans),
                                M2, N2, K,
                                &mzone, WORK, LDWORK,
                                V, LDV,
                                &zone, A2, LDA2
                                );
        PLASMA_CUBLAS_ASSERT(status);

        /*
         * A1 = A1 - B
         */
        for(j = 0; j < K; j++) {
          status = cublasDaxpy_v2( kaapi_cuda_cublas_handle(),
                                  M1, &mzone,
                                  &WORK[LDWORK*j], 1,
                                  &A1[LDA1*j], 1);
          PLASMA_CUBLAS_ASSERT(status);
        }
      }
    }
    else {
      coreblas_error(3, "Not implemented (RowMajor / Backward / Left or Right)");
      return PLASMA_ERR_NOT_SUPPORTED;
    }
  }
  return PLASMA_SUCCESS;
}
