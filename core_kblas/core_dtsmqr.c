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
*  CORE_dtsmqr overwrites the general complex M1-by-N1 tile A1 and
*  M2-by-N2 tile A2 with
*
*                        SIDE = 'L'        SIDE = 'R'
*    TRANS = 'N':         Q * | A1 |     | A1 A2 | * Q
*                             | A2 |
*
*    TRANS = 'C':      Q**T * | A1 |     | A1 A2 | * Q**T
*                             | A2 |
*
*  where Q is a complex unitary matrix defined as the product of k
*  elementary reflectors
*
*    Q = H(1) H(2) . . . H(k)
*
*  as returned by CORE_DTSQRT.
*
*******************************************************************************
*
* @param[in] side
*         @arg PlasmaLeft  : apply Q or Q**T from the Left;
*         @arg PlasmaRight : apply Q or Q**T from the Right.
*
* @param[in] trans
*         @arg PlasmaNoTrans   :  No transpose, apply Q;
*         @arg PlasmaTrans :  ConjTranspose, apply Q**T.
*
* @param[in] M1
*         The number of rows of the tile A1. M1 >= 0.
*
* @param[in] N1
*         The number of columns of the tile A1. N1 >= 0.
*
* @param[in] M2
*         The number of rows of the tile A2. M2 >= 0.
*         M2 = M1 if side == PlasmaRight.
*
* @param[in] N2
*         The number of columns of the tile A2. N2 >= 0.
*         N2 = N1 if side == PlasmaLeft.
*
* @param[in] K
*         The number of elementary reflectors whose product defines
*         the matrix Q.
*
* @param[in] IB
*         The inner-blocking size.  IB >= 0.
*
* @param[in,out] A1
*         On entry, the M1-by-N1 tile A1.
*         On exit, A1 is overwritten by the application of Q.
*
* @param[in] LDA1
*         The leading dimension of the array A1. LDA1 >= max(1,M1).
*
* @param[in,out] A2
*         On entry, the M2-by-N2 tile A2.
*         On exit, A2 is overwritten by the application of Q.
*
* @param[in] LDA2
*         The leading dimension of the tile A2. LDA2 >= max(1,M2).
*
* @param[in] V
*         The i-th row must contain the vector which defines the
*         elementary reflector H(i), for i = 1,2,...,k, as returned by
*         CORE_DTSQRT in the first k columns of its array argument V.
*
* @param[in] LDV
*         The leading dimension of the array V. LDV >= max(1,K).
*
* @param[out] T
*         The IB-by-N1 triangular factor T of the block reflector.
*         T is upper triangular by block (economic storage);
*         The rest of the array is not referenced.
*
* @param[in] LDT
*         The leading dimension of the array T. LDT >= IB.
*
* @param[out] WORK
*         Workspace array of size 
*             LDWORK-by-N1 if side == PlasmaLeft
*             LDWORK-by-IB if side == PlasmaRight
*
* @param[in] LDWORK
*         The leading dimension of the array WORK. 
*             LDWORK >= max(1,IB) if side == PlasmaLeft
*             LDWORK >= max(1,M1) if side == PlasmaRight
*
*******************************************************************************
*
* @return
*          \retval PLASMA_SUCCESS successful exit
*          \retval <0 if -i, the i-th argument had an illegal value
*
******************************************************************************/

int CORE_dtsmqr_cublas(int side, int trans,
                int M1, int N1, int M2, int N2, int K, int IB,
                double *A1, int LDA1,
                double *A2, int LDA2,
                double *V, int LDV,
                double *T, int LDT,
                double *WORK, int LDWORK)
{
  int i, i1, i3;
  int NQ, NW;
  int kb;
  int ic = 0;
  int jc = 0;
  int mi = M1;
  int ni = N1;
  
  /* Check input arguments */
  if ((side != PlasmaLeft) && (side != PlasmaRight)) {
    coreblas_error(1, "Illegal value of side");
    return -1;
  }
  
  /* NQ is the order of Q */
  if (side == PlasmaLeft) {
    NQ = M2;
    NW = IB;
  }
  else {
    NQ = N2;
    NW = M1;
  }
  
  if ((trans != PlasmaNoTrans) && (trans != PlasmaTrans)) {
    coreblas_error(2, "Illegal value of trans");
    return -2;
  }
  if (M1 < 0) {
    coreblas_error(3, "Illegal value of M1");
    return -3;
  }
  if (N1 < 0) {
    coreblas_error(4, "Illegal value of N1");
    return -4;
  }
  if ( (M2 < 0) || 
      ( (M2 != M1) && (side == PlasmaRight) ) ){
    coreblas_error(5, "Illegal value of M2");
    return -5;
  }
  if ( (N2 < 0) || 
      ( (N2 != N1) && (side == PlasmaLeft) ) ){
    coreblas_error(6, "Illegal value of N2");
    return -6;
  }
  if ((K < 0) || 
      ( (side == PlasmaLeft)  && (K > M1) ) ||
      ( (side == PlasmaRight) && (K > N1) ) ) {
    coreblas_error(7, "Illegal value of K");
    return -7;
  }
  if (IB < 0) {
    coreblas_error(8, "Illegal value of IB");
    return -8;
  }
  if (LDA1 < max(1,M1)){
    coreblas_error(10, "Illegal value of LDA1");
    return -10;
  }
  if (LDA2 < max(1,M2)){
    coreblas_error(12, "Illegal value of LDA2");
    return -12;
  }
  if (LDV < max(1,NQ)){
    coreblas_error(14, "Illegal value of LDV");
    return -14;
  }
  if (LDT < max(1,IB)){
    coreblas_error(16, "Illegal value of LDT");
    return -16;
  }
  if (LDWORK < max(1,NW)){
    coreblas_error(18, "Illegal value of LDWORK");
    return -18;
  }
  
  /* Quick return */
  if ((M1 == 0) || (N1 == 0) || (M2 == 0) || (N2 == 0) || (K == 0) || (IB == 0))
    return PLASMA_SUCCESS;
  
  if (((side == PlasmaLeft)  && (trans != PlasmaNoTrans))
      || ((side == PlasmaRight) && (trans == PlasmaNoTrans))) {
    i1 = 0;
    i3 = IB;
  }
  else {
    i1 = ((K-1) / IB)*IB;
    i3 = -IB;
  }
  
  for(i = i1; (i > -1) && (i < K); i += i3) {
    kb = min(IB, K-i);
    
    if (side == PlasmaLeft) {
      /*
       * H or H' is applied to C(i:m,1:n)
       */
      mi = M1 - i;
      ic = i;
    }
    else {
      /*
       * H or H' is applied to C(1:m,i:n)
       */
      ni = N1 - i;
      jc = i;
    }
    
    /*
     * Apply H or H'
     */
    CORE_dtsrfb_cublas(
                side, trans, PlasmaForward, PlasmaColumnwise,
                mi, ni, M2, N2, kb,
                &A1[LDA1*jc+ic], LDA1,
                A2, LDA2,
                &V[LDV*i], LDV,
                &T[LDT*i], LDT,
                WORK, LDWORK);
  }
  return PLASMA_SUCCESS;
}

void CORE_dtsmqr_quark_cublas(Quark *quark)
{
  int side;
  int trans;
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
  double *V;
  int ldv;
  double *T;
  int ldt;
  double *WORK;
  int ldwork;
  
  quark_unpack_args_18(quark, side, trans, m1, n1, m2, n2, k, ib,
                       A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
  CORE_dtsmqr_cublas(side, trans, m1, n1, m2, n2, k, ib,
              A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
}