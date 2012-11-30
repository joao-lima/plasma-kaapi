/**
 *
 * @file core_cher2k.c
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
#pragma weak CORE_cher2k = PCORE_cher2k
#define CORE_cher2k PCORE_cher2k
#endif
void CORE_cher2k(int uplo, int trans,
                 int N, int K,
                 PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int LDA,
                 PLASMA_Complex32_t *B, int LDB,
                 float beta, PLASMA_Complex32_t *C, int LDC)
{
    cblas_cher2k(
        CblasColMajor,
        (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
        N, K,
        CBLAS_SADDR(alpha), A, LDA, B, LDB,
        beta, C, LDC);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_cher2k(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int trans,
                       int n, int k, int nb,
                       PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int lda,
                       PLASMA_Complex32_t *B, int ldb,
                       float beta, PLASMA_Complex32_t *C, int ldc)
{
    DAG_CORE_HER2K;
    QUARK_Insert_Task(quark, CORE_cher2k_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &trans,     VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(PLASMA_Complex32_t),         &alpha,     VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    B,                 INPUT,
        sizeof(int),                        &ldb,       VALUE,
        sizeof(float),                     &beta,      VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    C,                 INOUT,
        sizeof(int),                        &ldc,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cher2k_quark = PCORE_cher2k_quark
#define CORE_cher2k_quark PCORE_cher2k_quark
#endif
void CORE_cher2k_quark(Quark *quark)
{
    int uplo;
    int trans;
    int n;
    int k;
    PLASMA_Complex32_t alpha;
    PLASMA_Complex32_t *A;
    int lda;
    PLASMA_Complex32_t *B;
    int ldb;
    float beta;
    PLASMA_Complex32_t *C;
    int ldc;

    quark_unpack_args_12(quark, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    CORE_cher2k(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
#endif
