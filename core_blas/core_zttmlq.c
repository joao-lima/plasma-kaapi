/**
 *
 * @file core_zttmlq.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Dulceneia Becker
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 *  CORE_zttmlq overwrites the general complex M1-by-N1 tile A1 and
 *  M2-by-N2 tile A2 (N1 == N2) with
 *
 *                        SIDE = 'L'        SIDE = 'R'
 *    TRANS = 'N':         Q * | A1 |       | A1 | * Q
 *                             | A2 |       | A2 |
 *
 *    TRANS = 'C':      Q**H * | A1 |       | A1 | * Q**H
 *                             | A2 |       | A2 |
 *
 *  where Q is a complex unitary matrix defined as the product of k
 *  elementary reflectors
 *
 *    Q = H(1) H(2) . . . H(k)
 *
 *  as returned by CORE_zttqrt.
 *
 *  Only LEFT is supported !!!
 *
 *******************************************************************************
 *
 * @param[in] side
 *         @arg PlasmaLeft  : apply Q or Q**H from the Left;
 *         @arg PlasmaRight : apply Q or Q**H from the Right.
 *
 * @param[in] trans
 *         @arg PlasmaNoTrans   :  No transpose, apply Q;
 *         @arg PlasmaConjTrans :  ConjTranspose, apply Q**H.
 *
 * @param[in] M1
 *         The number of rows of the tile A1. M1 >= 0.
 *
 * @param[in] N1
 *         The number of columns of the tile A1. N1 >= 0.
 *
 * @param[in] M2
 *         The number of rows of the tile A2. M2 >= 0.
 *
 * @param[in] N2
 *         The number of columns of the tile A2. N2 >= 0.
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
 *         CORE_ZTTQRT in the first k rows of its array argument V.
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
 *         Workspace array of size LDWORK-by-N1.
 *
 * @param[in] LDWORK
 *         The dimension of the array WORK. LDWORK >= max(1,IB).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zttmlq = PCORE_zttmlq
#define CORE_zttmlq PCORE_zttmlq
#define CORE_zttrfb PCORE_zttrfb
int  CORE_zttrfb(int side, int trans, int direct, int storev,
                 int M1, int N1, int M2, int N2, int IB,
                 PLASMA_Complex64_t *A1, int LDA1,
                 PLASMA_Complex64_t *A2, int LDA2,
                 PLASMA_Complex64_t *V, int LDV,
                 PLASMA_Complex64_t *T, int LDT,
                 PLASMA_Complex64_t *WORK, int LDWORK);
#endif
int CORE_zttmlq(int side, int trans,
                int M1, int N1, int M2, int N2, int K, int IB,
                PLASMA_Complex64_t *A1, int LDA1,
                PLASMA_Complex64_t *A2, int LDA2,
                PLASMA_Complex64_t *V, int LDV,
                PLASMA_Complex64_t *T, int LDT,
                PLASMA_Complex64_t *WORK, int LDWORK)
{
    int i,  i1, i3;
    int NQ, NW;
    int kb;
    int ic = 0;
    int jc = 0;
    int mi1 = M1;
    int mi2 = M2;
    int ni1 = N1;
    int ni2 = N2;

    /* Check input arguments */
    if ((side != PlasmaLeft) && (side != PlasmaRight)) {
        coreblas_error(1, "Illegal value of side");
        return -1;
    }

    /* NQ is the order of Q */
    if (side == PlasmaLeft) {
        NQ = N2;
        NW = IB;
    }
    else {
        NQ = M2;
        NW = N1;
    }

    if ((trans != PlasmaNoTrans) && (trans != PlasmaConjTrans)) {
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
    if ((M2 < 0) ||
        ( (side == PlasmaRight) && (M1 != M2) ) ) {
        coreblas_error(5, "Illegal value of M2");
        return -5;
    }
    if ((N2 < 0) ||
        ( (side == PlasmaLeft) && (N1 != N2) ) ) {
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

    if (((side == PlasmaLeft) && (trans == PlasmaNoTrans))
        || ((side == PlasmaRight) && (trans != PlasmaNoTrans))) {
        i1 = 0;
        i3 = IB;
    }
    else {
        i1 = ( ( K-1 ) / IB )*IB;
        i3 = -IB;
    }

    /* Transpose */
    if (trans == PlasmaNoTrans) {
        trans = PlasmaConjTrans;
    }
    else {
        trans = PlasmaNoTrans;
    }

    for (i = i1; (i > -1) && (i < K); i+=i3) {
        kb = min(IB, K-i);

        if (side == PlasmaLeft) {
            mi1 = M1 - i;
            mi2 = i + kb;
            ic = i;
        }
        else {
            /* ni2 must be the number of rows of V */
            ni1 = kb;
            ni2 = i + kb;
            /*if ((K == IB) && (NeqIB > 0)) { ni2 = N2 - i; }*/
            jc = i;
        }

        /*
         * H or H' is applied to C(i:m,1:n)
         */
        CORE_zttrfb(
            side, trans, PlasmaForward, PlasmaRowwise,
            mi1, ni1, mi2, ni2, kb,
            &A1[LDA1*jc+ic], LDA1,
            A2, LDA2,
            &V[i], LDV,
            &T[LDT*i], LDT,
            WORK, LDWORK);
    }
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zttmlq(Quark *quark, Quark_Task_Flags *task_flags,
                       int side, int trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       PLASMA_Complex64_t *A1, int lda1,
                       PLASMA_Complex64_t *A2, int lda2,
                       PLASMA_Complex64_t *V, int ldv,
                       PLASMA_Complex64_t *T, int ldt)
{
    int ldwork = side == PlasmaLeft ? ib : nb;

    DAG_CORE_TTMLQ;
    QUARK_Insert_Task(quark, CORE_zttmlq_quark, task_flags,
        sizeof(PLASMA_enum),                &side,  VALUE,
        sizeof(PLASMA_enum),                &trans, VALUE,
        sizeof(int),                        &m1,    VALUE,
        sizeof(int),                        &n1,    VALUE,
        sizeof(int),                        &m2,    VALUE,
        sizeof(int),                        &n2,    VALUE,
        sizeof(int),                        &k,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A1,            INOUT,
        sizeof(int),                        &lda1,  VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A2,            INOUT,
        sizeof(int),                        &lda2,  VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    V,             INPUT|QUARK_REGION_D|QUARK_REGION_L,
        sizeof(int),                        &ldv,   VALUE,
        sizeof(PLASMA_Complex64_t)*ib*nb,    T,             INPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(PLASMA_Complex64_t)*ib*nb,    NULL,          SCRATCH,
        sizeof(int),                        &ldwork,    VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zttmlq_quark = PCORE_zttmlq_quark
#define CORE_zttmlq_quark PCORE_zttmlq_quark
#endif
void CORE_zttmlq_quark(Quark *quark)
{
    int side;
    int trans;
    int m1;
    int n1;
    int m2;
    int n2;
    int k;
    int ib;
    PLASMA_Complex64_t *A1;
    int lda1;
    PLASMA_Complex64_t *A2;
    int lda2;
    PLASMA_Complex64_t *V;
    int ldv;
    PLASMA_Complex64_t *T;
    int ldt;
    PLASMA_Complex64_t *WORK;
    int ldwork;

    quark_unpack_args_18(quark, side, trans, m1, n1, m2, n2, k, ib, 
                         A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
    CORE_zttmlq(side, trans, m1, n1, m2, n2, k, ib, A1, lda1, 
                A2, lda2, V, ldv, T, ldt, WORK, ldwork);
}
