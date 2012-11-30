/**
 *
 * @file core_dgessm.c
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
 * @generated d Thu Sep 15 12:08:58 2011
 *
 **/
#include <cblas.h>
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 *  CORE_dgessm applies the factor L computed by CORE_dgetrf_incpiv to
 *  a complex M-by-N tile A.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the tile A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A.  N >= 0.
 *
 * @param[in] K
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in] IPIV
 *         as returned by CORE_dgetrf_incpiv.
 *
 * @param[in] L
 *         The NB-by-NB lower triangular tile.
 *
 * @param[in] LDL
 *         The leading dimension of the array L.  LDL >= max(1,NB).
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile A.
 *         On exit, updated by the application of L.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *         \retval PLASMA_SUCCESS successful exit
 *         \retval <0 if INFO = -k, the k-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dgessm = PCORE_dgessm
#define CORE_dgessm PCORE_dgessm
#endif
int CORE_dgessm(int M, int N, int K, int IB,
                int *IPIV,
                double *L, int LDL,
                double *A, int LDA)
{
    static double zone  =  1.0;
    static double mzone = -1.0;
    static int                ione  =  1;

    int i, sb;
    int tmp, tmp2;

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

    for(i = 0; i < K; i += IB) {
        sb = min(IB, K-i);
        /*
         * Apply interchanges to columns I*IB+1:IB*( I+1 )+1.
         */
        tmp  = i+1;
        tmp2 = i+sb;
        LAPACKE_dlaswp_work(LAPACK_COL_MAJOR, N, A, LDA, tmp, tmp2, IPIV, ione);
        /*
         * Compute block row of U.
         */
        cblas_dtrsm(
            CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
            sb, N, (zone),
            &L[LDL*i+i], LDL,
            &A[i], LDA );

        if (i+sb < M) {
        /*
        * Update trailing submatrix.
        */
        cblas_dgemm(
            CblasColMajor, CblasNoTrans, CblasNoTrans,
            M-(i+sb), N, sb,
            (mzone), &L[LDL*i+(i+sb)], LDL,
            &A[i], LDA,
            (zone), &A[i+sb], LDA );
        }
    }
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dgessm(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int k, int ib, int nb,
                       int *IPIV,
                       double *L, int ldl,
                       double *A, int lda)
{
    DAG_CORE_GESSM;
    QUARK_Insert_Task(quark, CORE_dgessm_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &k,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(int)*nb,                      IPIV,          INPUT,
        sizeof(double)*nb*nb,    L,             INPUT | QUARK_REGION_L,
        sizeof(int),                        &ldl,   VALUE,
        sizeof(double)*nb*nb,    A,             INOUT,
        sizeof(int),                        &lda,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dgessm_quark = PCORE_dgessm_quark
#define CORE_dgessm_quark PCORE_dgessm_quark
#endif
void CORE_dgessm_quark(Quark *quark)
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
    CORE_dgessm(m, n, k, ib, IPIV, L, ldl, A, lda);
}
