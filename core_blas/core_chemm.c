/**
 *
 * @file core_chemm.c
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
#pragma weak CORE_chemm = PCORE_chemm
#define CORE_chemm PCORE_chemm
#endif
void CORE_chemm(int side, int uplo,
                int M, int N,
                PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int LDA,
                PLASMA_Complex32_t *B, int LDB,
                PLASMA_Complex32_t beta, PLASMA_Complex32_t *C, int LDC)
{
    cblas_chemm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        M, N,
        CBLAS_SADDR(alpha), A, LDA,
        B, LDB,
        CBLAS_SADDR(beta), C, LDC);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_chemm(Quark *quark, Quark_Task_Flags *task_flags,
                      int side, int uplo,
                      int m, int n, int nb,
                      PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int lda,
                      PLASMA_Complex32_t *B, int ldb,
                      PLASMA_Complex32_t beta, PLASMA_Complex32_t *C, int ldc)
{
    DAG_CORE_HEMM;
    QUARK_Insert_Task(quark, CORE_chemm_quark, task_flags,
        sizeof(PLASMA_enum),                &side,    VALUE,
        sizeof(PLASMA_enum),                &uplo,    VALUE,
        sizeof(int),                        &m,       VALUE,
        sizeof(int),                        &n,       VALUE,
        sizeof(PLASMA_Complex32_t),         &alpha,   VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A,               INPUT,
        sizeof(int),                        &lda,     VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    B,               INPUT,
        sizeof(int),                        &ldb,     VALUE,
        sizeof(PLASMA_Complex32_t),         &beta,    VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    C,               INOUT,
        sizeof(int),                        &ldc,     VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_chemm_quark = PCORE_chemm_quark
#define CORE_chemm_quark PCORE_chemm_quark
#endif
void CORE_chemm_quark(Quark *quark)
{
    int side;
    int uplo;
    int M;
    int N;
    PLASMA_Complex32_t alpha;
    PLASMA_Complex32_t *A;
    int LDA;
    PLASMA_Complex32_t *B;
    int LDB;
    PLASMA_Complex32_t beta;
    PLASMA_Complex32_t *C;
    int LDC;

    quark_unpack_args_12(quark, side, uplo, M, N, alpha, A, LDA, B, LDB, beta, C, LDC);
    cblas_chemm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        M, N,
        CBLAS_SADDR(alpha), A, LDA,
        B, LDB,
        CBLAS_SADDR(beta), C, LDC);
}
#endif
