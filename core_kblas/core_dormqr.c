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

/***************************************************************************//**
*
* @ingroup CORE_double
*
*  CORE_dormqr overwrites the general complex M-by-N tile C with
*
*                    SIDE = 'L'     SIDE = 'R'
*    TRANS = 'N':      Q * C          C * Q
*    TRANS = 'C':      Q**T * C       C * Q**T
*
*  where Q is a complex unitary matrix defined as the product of k
*  elementary reflectors
*
*    Q = H(1) H(2) . . . H(k)
*
*  as returned by CORE_dgeqrt. Q is of order M if SIDE = 'L' and of order N
*  if SIDE = 'R'.
*
*******************************************************************************
*
* @param[in] side
*         @arg PlasmaLeft  : apply Q or Q**T from the Left;
*         @arg PlasmaRight : apply Q or Q**T from the Right.
*
* @param[in] trans
*         @arg PlasmaNoTrans   :  No transpose, apply Q;
*         @arg PlasmaTrans :  Transpose, apply Q**T.
*
* @param[in] M
*         The number of rows of the tile C.  M >= 0.
*
* @param[in] N
*         The number of columns of the tile C.  N >= 0.
*
* @param[in] K
*         The number of elementary reflectors whose product defines
*         the matrix Q.
*         If SIDE = PlasmaLeft,  M >= K >= 0;
*         if SIDE = PlasmaRight, N >= K >= 0.
*
* @param[in] IB
*         The inner-blocking size.  IB >= 0.
*
* @param[in] A
*         Dimension:  (LDA,K)
*         The i-th column must contain the vector which defines the
*         elementary reflector H(i), for i = 1,2,...,k, as returned by
*         CORE_dgeqrt in the first k columns of its array argument A.
*
* @param[in] LDA
*         The leading dimension of the array A.
*         If SIDE = PlasmaLeft,  LDA >= max(1,M);
*         if SIDE = PlasmaRight, LDA >= max(1,N).
*
* @param[out] T
*         The IB-by-K triangular factor T of the block reflector.
*         T is upper triangular by block (economic storage);
*         The rest of the array is not referenced.
*
* @param[in] LDT
*         The leading dimension of the array T. LDT >= IB.
*
* @param[in,out] C
*         On entry, the M-by-N tile C.
*         On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
*
* @param[in] LDC
*         The leading dimension of the array C. LDC >= max(1,M).
*
* @param[in,out] WORK
*         On exit, if INFO = 0, WORK(1) returns the optimal LDWORK.
*
* @param[in] LDWORK
*         The dimension of the array WORK.
*         If SIDE = PlasmaLeft,  LDWORK >= max(1,N);
*         if SIDE = PlasmaRight, LDWORK >= max(1,M).
*
*******************************************************************************
*
* @return
*          \retval PLASMA_SUCCESS successful exit
*          \retval <0 if -i, the i-th argument had an illegal value
*
******************************************************************************/

int CORE_dormqr_cublas(int side, int trans,
                int M, int N, int K, int IB,
                double *A, int LDA,
                double *T, int LDT,
                double *C, int LDC,
                double *WORK, int LDWORK)
{
  int i, kb;
  int i1, i3;
  int nq, nw;
  int ic = 0;
  int jc = 0;
  int ni = N;
  int mi = M;
  
  /* Check input arguments */
  if ((side != PlasmaLeft) && (side != PlasmaRight)) {
    coreblas_error(1, "Illegal value of side");
    return -1;
  }
  /*
   * NQ is the order of Q and NW is the minimum dimension of WORK
   */
  if (side == PlasmaLeft) {
    nq = M;
    nw = N;
  }
  else {
    nq = N;
    nw = M;
  }
  
  if ((trans != PlasmaNoTrans) && (trans != PlasmaTrans)) {
    coreblas_error(2, "Illegal value of trans");
    return -2;
  }
  if (M < 0) {
    coreblas_error(3, "Illegal value of M");
    return -3;
  }
  if (N < 0) {
    coreblas_error(4, "Illegal value of N");
    return -4;
  }
  if ((K < 0) || (K > nq)) {
    coreblas_error(5, "Illegal value of K");
    return -5;
  }
  if ((IB < 0) || ( (IB == 0) && ((M > 0) && (N > 0)) )) {
    coreblas_error(6, "Illegal value of IB");
    return -6;
  }
  if ((LDA < max(1,nq)) && (nq > 0)) {
    coreblas_error(8, "Illegal value of LDA");
    return -8;
  }
  if ((LDC < max(1,M)) && (M > 0)) {
    coreblas_error(12, "Illegal value of LDC");
    return -12;
  }
  if ((LDWORK < max(1,nw)) && (nw > 0)) {
    coreblas_error(14, "Illegal value of LDWORK");
    return -14;
  }
  
  /* Quick return */
  if ((M == 0) || (N == 0) || (K == 0))
    return PLASMA_SUCCESS;
  
  if (((side == PlasmaLeft) && (trans != PlasmaNoTrans))
      || ((side == PlasmaRight) && (trans == PlasmaNoTrans))) {
    i1 = 0;
    i3 = IB;
  }
  else {
    i1 = ( ( K-1 ) / IB )*IB;
    i3 = -IB;
  }
  
  for(i = i1; (i >- 1) && (i < K); i+=i3 ) {
    kb = min(IB, K-i);
    
    if (side == PlasmaLeft) {
      /*
       * H or H' is applied to C(i:m,1:n)
       */
      mi = M - i;
      ic = i;
    }
    else {
      /*
       * H or H' is applied to C(1:m,i:n)
       */
      ni = N - i;
      jc = i;
    }
    /*
     * Apply H or H'
     */
    CORE_dlarfb_gemm_cublas(
                        side, trans,
                        PlasmaForward,
                        PlasmaColumnwise,
                        mi, ni, kb,
                        &A[LDA*i+i], LDA,
                        &T[LDT*i], LDT,
                        &C[LDC*jc+ic], LDC,
                        WORK, LDWORK);
  }
  return PLASMA_SUCCESS;
}

void CORE_dormqr_quark_cublas(Quark *quark)
{
  int side;
  int trans;
  int m;
  int n;
  int k;
  int ib;
  double *A;
  int lda;
  double *T;
  int ldt;
  double *C;
  int ldc;
  double *WORK;
  int ldwork;
  
  quark_unpack_args_14(quark, side, trans, m, n, k, ib,
                       A, lda, T, ldt, C, ldc, WORK, ldwork);
  CORE_dormqr_cublas(side, trans, m, n, k, ib,
              A, lda, T, ldt, C, ldc, WORK, ldwork);
}
