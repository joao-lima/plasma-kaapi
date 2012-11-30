/**
 *
 * @file core_dgetrf_incpiv.c
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
 * @generated d Thu Sep 15 12:08:57 2011
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 *  CORE_dgetrf_incpiv computes an LU factorization of a general M-by-N tile A
 *  using partial pivoting with row interchanges.
 *
 *  The factorization has the form
 *
 *    A = P * L * U
 *
 *  where P is a permutation matrix, L is lower triangular with unit
 *  diagonal elements (lower trapezoidal if m > n), and U is upper
 *  triangular (upper trapezoidal if m < n).
 *
 *  This is the right-looking Level 2.5 BLAS version of the algorithm.
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
 *         On entry, the M-by-N tile to be factored.
 *         On exit, the factors L and U from the factorization
 *         A = P*L*U; the unit diagonal elements of L are not stored.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 * @param[out] IPIV
 *         The pivot indices; for 1 <= i <= min(M,N), row i of the
 *         tile was interchanged with row IPIV(i).
 *
 * @param[out] INFO
 *         See returned value.
 *
 *******************************************************************************
 *
 * @return
 *         \retval PLASMA_SUCCESS successful exit
 *         \retval <0 if INFO = -k, the k-th argument had an illegal value
 *         \retval >0 if INFO = k, U(k,k) is exactly zero. The factorization
 *              has been completed, but the factor U is exactly
 *              singular, and division by zero will occur if it is used
 *              to solve a system of equations.
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dgetrf_incpiv = PCORE_dgetrf_incpiv
#define CORE_dgetrf_incpiv PCORE_dgetrf_incpiv
#define CORE_dgessm PCORE_dgessm
int  CORE_dgessm(int M, int N, int K, int IB,
                 int *IPIV,
                 double *L, int LDL,
                 double *A, int LDA);
#endif
int CORE_dgetrf_incpiv(int M, int N, int IB,
                double *A, int LDA,
                int *IPIV, int *INFO)
{
    int i, j, k, sb;
    int iinfo;

    /* Check input arguments */
    *INFO = 0;
    if (M < 0) {
        coreblas_error(1, "Illegal value of M");
        return -1;
    }
    if (N < 0) {
        coreblas_error(2, "Illegal value of N");
        return -2;
    }
    if (IB < 0) {
        coreblas_error(3, "Illegal value of IB");
        return -3;
    }
    if ((LDA < max(1,M)) && (M > 0)) {
        coreblas_error(5, "Illegal value of LDA");
        return -5;
    }

    /* Quick return */
    if ((M == 0) || (N == 0) || (IB == 0))
        return PLASMA_SUCCESS;

    k = min(M, N);

    for(i =0 ; i < k; i += IB) {
        sb = min(IB, k-i);
        /*
         * Factor diagonal and subdiagonal blocks and test for exact singularity.
         */
        iinfo = LAPACKE_dgetf2_work(LAPACK_COL_MAJOR, M-i, sb, &A[LDA*i+i], LDA, &IPIV[i]);
        /*
         * Adjust INFO and the pivot indices.
         */
        if((*INFO == 0) && (iinfo > 0))
            *INFO = iinfo + i;

        if (i+sb < N) {
            CORE_dgessm(
                M-i, N-(i+sb), sb, sb,
                &IPIV[i],
                &A[LDA*i+i], LDA,
                &A[LDA*(i+sb)+i], LDA);
        }

        for(j = i; j < i+sb; j++) {
            IPIV[j] = i + IPIV[j];
        }
    }
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dgetrf_incpiv(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       double *A, int lda,
                       int *IPIV,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       PLASMA_bool check_info, int iinfo)
{
    DAG_CORE_GETRF;
    QUARK_Insert_Task(quark, CORE_dgetrf_incpiv_quark, task_flags,
        sizeof(int),                        &m,             VALUE,
        sizeof(int),                        &n,             VALUE,
        sizeof(int),                        &ib,            VALUE,
        sizeof(double)*nb*nb,    A,                     INOUT,
        sizeof(int),                        &lda,           VALUE,
        sizeof(int)*nb,                      IPIV,                  OUTPUT,
        sizeof(PLASMA_sequence*),           &sequence,      VALUE,
        sizeof(PLASMA_request*),            &request,       VALUE,
        sizeof(PLASMA_bool),                &check_info,    VALUE,
        sizeof(int),                        &iinfo,         VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dgetrf_incpiv_quark = PCORE_dgetrf_incpiv_quark
#define CORE_dgetrf_incpiv_quark PCORE_dgetrf_incpiv_quark
#endif
void CORE_dgetrf_incpiv_quark(Quark *quark)
{
    int m;
    int n;
    int ib;
    double *A;
    int lda;
    int *IPIV;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    PLASMA_bool check_info;
    int iinfo;

    int info;

    quark_unpack_args_10(quark, m, n, ib, A, lda, IPIV, sequence, request, check_info, iinfo);
    CORE_dgetrf_incpiv(m, n, ib, A, lda, IPIV, &info);
    if (info != PLASMA_SUCCESS && check_info)
        plasma_sequence_flush(quark, sequence, request, iinfo+info);
}
