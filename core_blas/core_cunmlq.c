/**
 *
 * @file core_cunmlq.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @author Dulceneia Becker
 * @date 2010-11-15
 * @generated c Thu Sep 15 12:09:00 2011
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex32_t
 *
 *  CORE_cunmlq overwrites the general complex M-by-N tile C with
 *
 *                    SIDE = 'L'     SIDE = 'R'
 *    TRANS = 'N':      Q * C          C * Q
 *    TRANS = 'C':      Q**H * C       C * Q**H
 *
 *  where Q is a complex unitary matrix defined as the product of k
 *  elementary reflectors
 *
 *    Q = H(k) . . . H(2) H(1)
 *
 *  as returned by CORE_cgelqt. Q is of order M if SIDE = 'L' and of order N
 *  if SIDE = 'R'.
 *
 *******************************************************************************
 *
 * @param[in] side
 *         @arg PlasmaLeft  : apply Q or Q**H from the Left;
 *         @arg PlasmaRight : apply Q or Q**H from the Right.
 *
 * @param[in] trans
 *         @arg PlasmaNoTrans   :  No transpose, apply Q;
 *         @arg PlasmaConjTrans :  Transpose, apply Q**H.
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
 *         Dimension:  (LDA,M) if SIDE = PlasmaLeft,
 *                     (LDA,N) if SIDE = PlasmaRight,
 *         The i-th row must contain the vector which defines the
 *         elementary reflector H(i), for i = 1,2,...,k, as returned by
 *         CORE_cgelqt in the first k rows of its array argument A.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,K).
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
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cunmlq = PCORE_cunmlq
#define CORE_cunmlq PCORE_cunmlq
#endif
int CORE_cunmlq(int side, int trans,
                int M, int N, int K, int IB,
                PLASMA_Complex32_t *A, int LDA,
                PLASMA_Complex32_t *T, int LDT,
                PLASMA_Complex32_t *C, int LDC,
                PLASMA_Complex32_t *WORK, int LDWORK)
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

    if ((trans != PlasmaNoTrans) && (trans != PlasmaConjTrans)) {
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
    if ((LDA < max(1,K)) && (K > 0)) {
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

    if (((side == PlasmaLeft) && (trans == PlasmaNoTrans))
        || ((side == PlasmaRight) && (trans != PlasmaNoTrans))) {
        i1 = 0;
        i3 = IB;
    }
    else {
        i1 = ( ( K-1 ) / IB )*IB;
        i3 = -IB;
    }

    if( trans == PlasmaNoTrans) {
        trans = PlasmaConjTrans;
    }
    else {
        trans = PlasmaNoTrans;
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
        LAPACKE_clarfb_work(LAPACK_COL_MAJOR,
            lapack_const(side),
            lapack_const(trans),
            lapack_const(PlasmaForward),
            lapack_const(PlasmaRowwise),
            mi, ni, kb,
            &A[LDA*i+i], LDA,
            &T[LDT*i], LDT,
            &C[LDC*jc+ic], LDC,
            WORK, LDWORK);
    }
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_cunmlq(Quark *quark, Quark_Task_Flags *task_flags,
                       int side, int trans,
                       int m, int n, int k, int ib, int nb,
                       PLASMA_Complex32_t *A, int lda,
                       PLASMA_Complex32_t *T, int ldt,
                       PLASMA_Complex32_t *C, int ldc)
{
    DAG_CORE_UNMLQ;
    QUARK_Insert_Task(quark, CORE_cunmlq_quark, task_flags,
        sizeof(PLASMA_enum),                &side,  VALUE,
        sizeof(PLASMA_enum),                &trans, VALUE,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &k,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A,             INPUT | QUARK_REGION_U,
        sizeof(int),                        &lda,   VALUE,
        sizeof(PLASMA_Complex32_t)*ib*nb,    T,             INPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    C,             INOUT,
        sizeof(int),                        &ldc,   VALUE,
        sizeof(PLASMA_Complex32_t)*ib*nb,    NULL,          SCRATCH,
        sizeof(int),                        &nb,    VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cunmlq_quark = PCORE_cunmlq_quark
#define CORE_cunmlq_quark PCORE_cunmlq_quark
#endif
void CORE_cunmlq_quark(Quark *quark)
{
    int side;
    int trans;
    int m;
    int n;
    int k;
    int ib;
    PLASMA_Complex32_t *A;
    int lda;
    PLASMA_Complex32_t *T;
    int ldt;
    PLASMA_Complex32_t *C;
    int ldc;
    PLASMA_Complex32_t *WORK;
    int ldwork;

    quark_unpack_args_14(quark, side, trans, m, n, k, ib,
                         A, lda, T, ldt, C, ldc, WORK, ldwork);
    CORE_cunmlq(side, trans, m, n, k, ib, 
                A, lda, T, ldt, C, ldc, WORK, ldwork);
}
