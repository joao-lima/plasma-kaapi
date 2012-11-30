/**
 *
 * @file core_ctrsm.c
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
#pragma weak CORE_ctrsm = PCORE_ctrsm
#define CORE_ctrsm PCORE_ctrsm
#endif
void CORE_ctrsm(int side, int uplo,
                int transA, int diag,
                int M, int N,
                PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int LDA,
                PLASMA_Complex32_t *B, int LDB)
{
    cblas_ctrsm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        (CBLAS_TRANSPOSE)transA, (CBLAS_DIAG)diag,
        M, N,
        CBLAS_SADDR(alpha), A, LDA,
        B, LDB);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_ctrsm(Quark *quark, Quark_Task_Flags *task_flags,
                      int side, int uplo, int transA, int diag,
                      int m, int n, int nb,
                      PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int lda,
                      PLASMA_Complex32_t *B, int ldb)
{
    DAG_CORE_TRSM;
    QUARK_Insert_Task(quark, CORE_ctrsm_quark, task_flags,
        sizeof(PLASMA_enum),                &side,      VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &diag,      VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(PLASMA_Complex32_t),         &alpha,     VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    B,                 INOUT | LOCALITY,
        sizeof(int),                        &ldb,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ctrsm_quark = PCORE_ctrsm_quark
#define CORE_ctrsm_quark PCORE_ctrsm_quark
#endif
void CORE_ctrsm_quark(Quark *quark)
{
    int side;
    int uplo;
    int transA;
    int diag;
    int m;
    int n;
    PLASMA_Complex32_t alpha;
    PLASMA_Complex32_t *A;
    int lda;
    PLASMA_Complex32_t *B;
    int ldb;

    quark_unpack_args_11(quark, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
    cblas_ctrsm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        (CBLAS_TRANSPOSE)transA, (CBLAS_DIAG)diag,
        m, n,
        CBLAS_SADDR(alpha), A, lda,
        B, ldb);
}
