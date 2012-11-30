/**
 *
 * @file core_cherk.c
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
 * @generated c Thu Sep 15 12:09:00 2011
 *
 **/
#include "common.h"

#undef REAL
#define COMPLEX
#ifdef COMPLEX
/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex32_t
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cherk = PCORE_cherk
#define CORE_cherk PCORE_cherk
#endif
void CORE_cherk(int uplo, int trans,
                int N, int K,
                float alpha, PLASMA_Complex32_t *A, int LDA,
                float beta, PLASMA_Complex32_t *C, int LDC)
{
    cblas_cherk(
        CblasColMajor,
        (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
        N, K,
        alpha, A, LDA,
        beta, C, LDC);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_cherk(Quark *quark, Quark_Task_Flags *task_flags,
                      int uplo, int trans,
                      int n, int k, int nb,
                      float alpha, PLASMA_Complex32_t *A, int lda,
                      float beta, PLASMA_Complex32_t *C, int ldc)
{
    DAG_CORE_HERK;
    QUARK_Insert_Task(quark, CORE_cherk_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &trans,     VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(float),                     &alpha,     VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(float),                     &beta,      VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    C,                 INOUT,
        sizeof(int),                        &ldc,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cherk_quark = PCORE_cherk_quark
#define CORE_cherk_quark PCORE_cherk_quark
#endif
void CORE_cherk_quark(Quark *quark)
{
    int uplo;
    int trans;
    int n;
    int k;
    float alpha;
    PLASMA_Complex32_t *A;
    int lda;
    float beta;
    PLASMA_Complex32_t *C;
    int ldc;

    quark_unpack_args_10(quark, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
    cblas_cherk(
        CblasColMajor,
        (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
        n, k,
        alpha, A, lda,
        beta, C, ldc);
}
#endif
