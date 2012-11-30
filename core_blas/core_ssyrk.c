/**
 *
 * @file core_ssyrk.c
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
 * @generated s Thu Sep 15 12:08:58 2011
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_float
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ssyrk = PCORE_ssyrk
#define CORE_ssyrk PCORE_ssyrk
#endif
void CORE_ssyrk(int uplo, int trans,
                int N, int K,
                float alpha, float *A, int LDA,
                float beta, float *C, int LDC)
{
    cblas_ssyrk(
        CblasColMajor,
        (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
        N, K,
        (alpha), A, LDA,
        (beta), C, LDC);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_ssyrk(Quark *quark, Quark_Task_Flags *task_flags,
                      int uplo, int trans,
                      int n, int k, int nb,
                      float alpha, float *A, int lda,
                      float beta, float *C, int ldc)
{
    DAG_CORE_SYRK;
    QUARK_Insert_Task(quark, CORE_ssyrk_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &trans,     VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(float),         &alpha,     VALUE,
        sizeof(float)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(float),         &beta,      VALUE,
        sizeof(float)*nb*nb,    C,                 INOUT,
        sizeof(int),                        &ldc,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ssyrk_quark = PCORE_ssyrk_quark
#define CORE_ssyrk_quark PCORE_ssyrk_quark
#endif
void CORE_ssyrk_quark(Quark *quark)
{
    int uplo;
    int trans;
    int n;
    int k;
    float alpha;
    float *A;
    int lda;
    float beta;
    float *C;
    int ldc;

    quark_unpack_args_10(quark, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
    cblas_ssyrk(
        CblasColMajor,
        (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
        n, k,
        (alpha), A, lda,
        (beta), C, ldc);
}
