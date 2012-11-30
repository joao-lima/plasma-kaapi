/**
 *
 * @file core_daxpy.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Thu Sep 15 12:08:59 2011
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_daxpy = PCORE_daxpy
#define CORE_daxpy PCORE_daxpy
#endif
void CORE_daxpy(int M, int N,  double alpha,
                double *A, int LDA,
                double *B, int LDB)
{
    int j;

    if (M == LDA)
        cblas_daxpy(M*N, (alpha), A, 1, B, 1);
    else {
        for (j = 0; j < N; j++)
            cblas_daxpy(M, (alpha), &A[j*LDA], 1, &B[j*LDA], 1);
    }
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_daxpy(Quark *quark, Quark_Task_Flags *task_flags,
                      int m, int n, int nb, double alpha,
                      double *A, int lda,
                      double *B, int ldb)
{
    DAG_CORE_AXPY;
    QUARK_Insert_Task(quark, CORE_daxpy_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(double),         &alpha, VALUE,
        sizeof(double)*nb*nb,    A,             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(double)*nb*nb,    B,             INOUT,
        sizeof(int),                        &ldb,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_daxpy_quark = PCORE_daxpy_quark
#define CORE_daxpy_quark PCORE_daxpy_quark
#endif
void CORE_daxpy_quark(Quark *quark)
{
    int M;
    int N;
    double alpha;
    double *A;
    int LDA;
    double *B;
    int LDB;

    int j;

    quark_unpack_args_7(quark, M, N, alpha, A, LDA, B, LDB);
    if (M == LDA)
        cblas_daxpy(M*N, (alpha), A, 1, B, 1);
    else {
        for (j = 0; j < N; j++)
            cblas_daxpy(M, (alpha), &A[j*LDA], 1, &B[j*LDA], 1);
    }
}

