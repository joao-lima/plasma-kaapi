/**
 *
 * @file core_dgelqt.c
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
 * @generated d Thu Sep 15 12:08:56 2011
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 *  CORE_dgelqt - computes a LQ factorization of a complex M-by-N tile A: A = L * Q.
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
#pragma weak CORE_dgelqt = PCORE_dgelqt
#define CORE_dgelqt PCORE_dgelqt
#endif
int CORE_dgelqt(int M, int N, int IB,
                double *A, int LDA,
                double *T, int LDT,
                double *TAU,
                double *WORK)
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

        iinfo = LAPACKE_dgelq2_work(LAPACK_COL_MAJOR, sb, N-i, 
                                    &A[LDA*i+i], LDA, &TAU[i], WORK);

        LAPACKE_dlarft_work(LAPACK_COL_MAJOR,
            lapack_const(PlasmaForward),
            lapack_const(PlasmaRowwise),
            N-i, sb,
            &A[LDA*i+i], LDA, &TAU[i],
            &T[LDT*i], LDT);
        
        if (M > i+sb) {
            LAPACKE_dlarfb_work(
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
void QUARK_CORE_dgelqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       double *A, int lda,
                       double *T, int ldt)
{
    DAG_CORE_GELQT;
    QUARK_Insert_Task(quark, CORE_dgelqt_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(double)*nb*nb,    A,             INOUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(double)*ib*nb,    T,             OUTPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(double)*nb,       NULL,          SCRATCH,
        sizeof(double)*ib*nb,    NULL,          SCRATCH,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dgelqt_quark = PCORE_dgelqt_quark
#define CORE_dgelqt_quark PCORE_dgelqt_quark
#endif
void CORE_dgelqt_quark(Quark *quark)
{
    int m;
    int n;
    int ib;
    double *A;
    int lda;
    double *T;
    int ldt;
    double *TAU;
    double *WORK;

    quark_unpack_args_9(quark, m, n, ib, A, lda, T, ldt, TAU, WORK);
    CORE_dgelqt(m, n, ib, A, lda, T, ldt, TAU, WORK);
}
