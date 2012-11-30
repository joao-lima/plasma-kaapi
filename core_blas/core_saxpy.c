/**
 *
 * @file core_saxpy.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Thu Sep 15 12:08:59 2011
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_float
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_saxpy = PCORE_saxpy
#define CORE_saxpy PCORE_saxpy
#endif
void CORE_saxpy(int M, int N,  float alpha,
                float *A, int LDA,
                float *B, int LDB)
{
    int j;

    if (M == LDA)
        cblas_saxpy(M*N, (alpha), A, 1, B, 1);
    else {
        for (j = 0; j < N; j++)
            cblas_saxpy(M, (alpha), &A[j*LDA], 1, &B[j*LDA], 1);
    }
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_saxpy(Quark *quark, Quark_Task_Flags *task_flags,
                      int m, int n, int nb, float alpha,
                      float *A, int lda,
                      float *B, int ldb)
{
    DAG_CORE_AXPY;
    QUARK_Insert_Task(quark, CORE_saxpy_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(float),         &alpha, VALUE,
        sizeof(float)*nb*nb,    A,             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(float)*nb*nb,    B,             INOUT,
        sizeof(int),                        &ldb,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_saxpy_quark = PCORE_saxpy_quark
#define CORE_saxpy_quark PCORE_saxpy_quark
#endif
void CORE_saxpy_quark(Quark *quark)
{
    int M;
    int N;
    float alpha;
    float *A;
    int LDA;
    float *B;
    int LDB;

    int j;

    quark_unpack_args_7(quark, M, N, alpha, A, LDA, B, LDB);
    if (M == LDA)
        cblas_saxpy(M*N, (alpha), A, 1, B, 1);
    else {
        for (j = 0; j < N; j++)
            cblas_saxpy(M, (alpha), &A[j*LDA], 1, &B[j*LDA], 1);
    }
}

