/**
 *
 * @file core_ctrmm.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Julien Langou
 * @author Henricus Bouwmeester
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
#pragma weak CORE_ctrmm = PCORE_ctrmm
#define CORE_ctrmm PCORE_ctrmm
#endif
void CORE_ctrmm(int side, int uplo,
                int transA, int diag,
                int M, int N, 
                PLASMA_Complex32_t alpha, 
                PLASMA_Complex32_t *A, int LDA,
                PLASMA_Complex32_t *B, int LDB)
{
    cblas_ctrmm(
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
void QUARK_CORE_ctrmm(Quark *quark, Quark_Task_Flags *task_flags,
                      int side, int uplo, int transA, int diag,
                      int m, int n, int nb,
                      PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int lda,
                      PLASMA_Complex32_t *B, int ldb)
{
    DAG_CORE_TRMM;
    QUARK_Insert_Task(quark, CORE_ctrmm_quark, task_flags,
        sizeof(PLASMA_enum),                &side,      VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &diag,      VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(PLASMA_Complex32_t),         &alpha,     VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    B,                 INOUT,
        sizeof(int),                        &ldb,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ctrmm_quark = PCORE_ctrmm_quark
#define CORE_ctrmm_quark PCORE_ctrmm_quark
#endif
void CORE_ctrmm_quark(Quark *quark)
{
    int side;
    int uplo;
    int transA;
    int diag;
    int M;
    int N;
    PLASMA_Complex32_t alpha;
    PLASMA_Complex32_t *A;
    int LDA;
    PLASMA_Complex32_t *B;
    int LDB;

    quark_unpack_args_11(quark, side, uplo, transA, diag, M, N, alpha, A, LDA, B, LDB);
    cblas_ctrmm(
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
void QUARK_CORE_ctrmm_p2(Quark *quark, Quark_Task_Flags *task_flags,
                         int side, int uplo, int transA, int diag,
                         int m, int n, int nb,
                         PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int lda,
                         PLASMA_Complex32_t **B, int ldb)
{
    DAG_CORE_TRMM;
    QUARK_Insert_Task(quark, CORE_ctrmm_p2_quark, task_flags,
        sizeof(PLASMA_enum),                &side,      VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &diag,      VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(PLASMA_Complex32_t),         &alpha,     VALUE,
        sizeof(PLASMA_Complex32_t)*lda*nb,   A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(PLASMA_Complex32_t*),         B,                 INOUT,
        sizeof(int),                        &ldb,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ctrmm_p2_quark = PCORE_ctrmm_p2_quark
#define CORE_ctrmm_p2_quark PCORE_ctrmm_p2_quark
#endif
void CORE_ctrmm_p2_quark(Quark *quark)
{
    int side;
    int uplo;
    int transA;
    int diag;
    int M;
    int N;
    PLASMA_Complex32_t alpha;
    PLASMA_Complex32_t *A;
    int LDA;
    PLASMA_Complex32_t **B;
    int LDB;

    quark_unpack_args_11(quark, side, uplo, transA, diag, M, N, alpha, A, LDA, B, LDB);
    cblas_ctrmm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        (CBLAS_TRANSPOSE)transA, (CBLAS_DIAG)diag,
        M, N,
        CBLAS_SADDR(alpha), A, LDA,
        *B, LDB);
}
