/**
 *
 * @file core_csyrk.c
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
 * @generated c Thu Sep 15 12:08:58 2011
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex32_t
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_csyrk = PCORE_csyrk
#define CORE_csyrk PCORE_csyrk
#endif
void CORE_csyrk(int uplo, int trans,
                int N, int K,
                PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int LDA,
                PLASMA_Complex32_t beta, PLASMA_Complex32_t *C, int LDC)
{
    cblas_csyrk(
        CblasColMajor,
        (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
        N, K,
        CBLAS_SADDR(alpha), A, LDA,
        CBLAS_SADDR(beta), C, LDC);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_csyrk(Quark *quark, Quark_Task_Flags *task_flags,
                      int uplo, int trans,
                      int n, int k, int nb,
                      PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int lda,
                      PLASMA_Complex32_t beta, PLASMA_Complex32_t *C, int ldc)
{
    DAG_CORE_SYRK;
    QUARK_Insert_Task(quark, CORE_csyrk_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &trans,     VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(PLASMA_Complex32_t),         &alpha,     VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(PLASMA_Complex32_t),         &beta,      VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    C,                 INOUT,
        sizeof(int),                        &ldc,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_csyrk_quark = PCORE_csyrk_quark
#define CORE_csyrk_quark PCORE_csyrk_quark
#endif
void CORE_csyrk_quark(Quark *quark)
{
    int uplo;
    int trans;
    int n;
    int k;
    PLASMA_Complex32_t alpha;
    PLASMA_Complex32_t *A;
    int lda;
    PLASMA_Complex32_t beta;
    PLASMA_Complex32_t *C;
    int ldc;

    quark_unpack_args_10(quark, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
    cblas_csyrk(
        CblasColMajor,
        (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
        n, k,
        CBLAS_SADDR(alpha), A, lda,
        CBLAS_SADDR(beta), C, ldc);
}
