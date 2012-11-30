/**
 *
 * @file core_sgelqt.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @date 2010-11-15
 * @generated s Thu Sep 15 12:08:57 2011
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_float
 *
 *  CORE_sgelqt - computes a LQ factorization of a complex M-by-N tile A: A = L * Q.
 *
 *  The tile Q is represented as a product of elementary reflectors
 *
 *    Q = H(k)' . . . H(2)' H(1)', where k = min(M,N).
 *
 *  Each H(i) has the form
 *
 *    H(i) = I - tau * v * v'
 *
 *  where tau is a complex scalar, and v is a complex vector with
 *  v(1:i-1) = 0 and v(i) = 1; g(v(i+1:n)) is stored on exit in
 *  A(i,i+1:n), and tau in TAU(i).
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the tile A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A.  N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile A.
 *         On exit, the elements on and below the diagonal of the array
 *         contain the M-by-min(M,N) lower trapezoidal tile L (L is
 *         lower triangular if M <= N); the elements above the diagonal,
 *         with the array TAU, represent the unitary tile Q as a
 *         product of elementary reflectors (see Further Details).
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 * @param[out] T
 *         The IB-by-N triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[out] TAU
 *         The scalar factors of the elementary reflectors (see Further
 *         Details).
 *
 * @param[out] WORK
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sgelqt = PCORE_sgelqt
#define CORE_sgelqt PCORE_sgelqt
#endif
int CORE_sgelqt(int M, int N, int IB,
                float *A, int LDA,
                float *T, int LDT,
                float *TAU,
                float *WORK)
{
    int i, k, sb;
    int iinfo;

    /* Check input arguments */
    if (M < 0) {
        coreblas_error(1, "Illegal value of M");
        return -1;
    }
    if (N < 0) {
        coreblas_error(2, "Illegal value of N");
        return -2;
    }
    if ((IB < 0) || ( (IB == 0) && ((M > 0) && (N > 0)) )) {
        coreblas_error(3, "Illegal value of IB");
        return -3;
    }
    if ((LDA < max(1,M)) && (M > 0)) {
        coreblas_error(5, "Illegal value of LDA");
        return -5;
    }
    if ((LDT < max(1,IB)) && (IB > 0)) {
        coreblas_error(7, "Illegal value of LDT");
        return -7;
    }

    /* Quick return */
    if ((M == 0) || (N == 0) || (IB == 0))
        return PLASMA_SUCCESS;

    k = min(M, N);

    for(i = 0; i < k; i += IB) {
        sb = min(IB, k-i);

        iinfo = LAPACKE_sgelq2_work(LAPACK_COL_MAJOR, sb, N-i, 
                                    &A[LDA*i+i], LDA, &TAU[i], WORK);

        LAPACKE_slarft_work(LAPACK_COL_MAJOR,
            lapack_const(PlasmaForward),
            lapack_const(PlasmaRowwise),
            N-i, sb,
            &A[LDA*i+i], LDA, &TAU[i],
            &T[LDT*i], LDT);
        
        if (M > i+sb) {
            LAPACKE_slarfb_work(
                LAPACK_COL_MAJOR,
                lapack_const(PlasmaRight),
                lapack_const(PlasmaNoTrans),
                lapack_const(PlasmaForward),
                lapack_const(PlasmaRowwise),
                M-i-sb, N-i, sb,
                &A[LDA*i+i],      LDA,
                &T[LDT*i],        LDT,
                &A[LDA*i+(i+sb)], LDA,
                WORK, M-i-sb);
        }
    }
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sgelqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       float *A, int lda,
                       float *T, int ldt)
{
    DAG_CORE_GELQT;
    QUARK_Insert_Task(quark, CORE_sgelqt_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(float)*nb*nb,    A,             INOUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(float)*ib*nb,    T,             OUTPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(float)*nb,       NULL,          SCRATCH,
        sizeof(float)*ib*nb,    NULL,          SCRATCH,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sgelqt_quark = PCORE_sgelqt_quark
#define CORE_sgelqt_quark PCORE_sgelqt_quark
#endif
void CORE_sgelqt_quark(Quark *quark)
{
    int m;
    int n;
    int ib;
    float *A;
    int lda;
    float *T;
    int ldt;
    float *TAU;
    float *WORK;

    quark_unpack_args_9(quark, m, n, ib, A, lda, T, ldt, TAU, WORK);
    CORE_sgelqt(m, n, ib, A, lda, T, ldt, TAU, WORK);
}
