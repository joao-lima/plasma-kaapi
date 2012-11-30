/**
 *
 * @file core_zaxpy.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zaxpy = PCORE_zaxpy
#define CORE_zaxpy PCORE_zaxpy
#endif
void CORE_zaxpy(int M, int N,  PLASMA_Complex64_t alpha,
                PLASMA_Complex64_t *A, int LDA,
                PLASMA_Complex64_t *B, int LDB)
{
    int j;

    if (M == LDA)
        cblas_zaxpy(M*N, CBLAS_SADDR(alpha), A, 1, B, 1);
    else {
        for (j = 0; j < N; j++)
            cblas_zaxpy(M, CBLAS_SADDR(alpha), &A[j*LDA], 1, &B[j*LDA], 1);
    }
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zaxpy(Quark *quark, Quark_Task_Flags *task_flags,
                      int m, int n, int nb, PLASMA_Complex64_t alpha,
                      PLASMA_Complex64_t *A, int lda,
                      PLASMA_Complex64_t *B, int ldb)
{
    DAG_CORE_AXPY;
    QUARK_Insert_Task(quark, CORE_zaxpy_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(PLASMA_Complex64_t),         &alpha, VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A,             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    B,             INOUT,
        sizeof(int),                        &ldb,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zaxpy_quark = PCORE_zaxpy_quark
#define CORE_zaxpy_quark PCORE_zaxpy_quark
#endif
void CORE_zaxpy_quark(Quark *quark)
{
    int M;
    int N;
    PLASMA_Complex64_t alpha;
    PLASMA_Complex64_t *A;
    int LDA;
    PLASMA_Complex64_t *B;
    int LDB;

    int j;

    quark_unpack_args_7(quark, M, N, alpha, A, LDA, B, LDB);
    if (M == LDA)
        cblas_zaxpy(M*N, CBLAS_SADDR(alpha), A, 1, B, 1);
    else {
        for (j = 0; j < N; j++)
            cblas_zaxpy(M, CBLAS_SADDR(alpha), &A[j*LDA], 1, &B[j*LDA], 1);
    }
}

