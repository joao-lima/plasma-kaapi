/**
 *
 * @file core_caxpy.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated c Thu Sep 15 12:08:59 2011
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex32_t
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_caxpy = PCORE_caxpy
#define CORE_caxpy PCORE_caxpy
#endif
void CORE_caxpy(int M, int N,  PLASMA_Complex32_t alpha,
                PLASMA_Complex32_t *A, int LDA,
                PLASMA_Complex32_t *B, int LDB)
{
    int j;

    if (M == LDA)
        cblas_caxpy(M*N, CBLAS_SADDR(alpha), A, 1, B, 1);
    else {
        for (j = 0; j < N; j++)
            cblas_caxpy(M, CBLAS_SADDR(alpha), &A[j*LDA], 1, &B[j*LDA], 1);
    }
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_caxpy(Quark *quark, Quark_Task_Flags *task_flags,
                      int m, int n, int nb, PLASMA_Complex32_t alpha,
                      PLASMA_Complex32_t *A, int lda,
                      PLASMA_Complex32_t *B, int ldb)
{
    DAG_CORE_AXPY;
    QUARK_Insert_Task(quark, CORE_caxpy_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(PLASMA_Complex32_t),         &alpha, VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A,             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    B,             INOUT,
        sizeof(int),                        &ldb,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_caxpy_quark = PCORE_caxpy_quark
#define CORE_caxpy_quark PCORE_caxpy_quark
#endif
void CORE_caxpy_quark(Quark *quark)
{
    int M;
    int N;
    PLASMA_Complex32_t alpha;
    PLASMA_Complex32_t *A;
    int LDA;
    PLASMA_Complex32_t *B;
    int LDB;

    int j;

    quark_unpack_args_7(quark, M, N, alpha, A, LDA, B, LDB);
    if (M == LDA)
        cblas_caxpy(M*N, CBLAS_SADDR(alpha), A, 1, B, 1);
    else {
        for (j = 0; j < N; j++)
            cblas_caxpy(M, CBLAS_SADDR(alpha), &A[j*LDA], 1, &B[j*LDA], 1);
    }
}

