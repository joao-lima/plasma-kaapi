/**
 *
 * @file core_strsm.c
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
#pragma weak CORE_strsm = PCORE_strsm
#define CORE_strsm PCORE_strsm
#endif
void CORE_strsm(int side, int uplo,
                int transA, int diag,
                int M, int N,
                float alpha, float *A, int LDA,
                float *B, int LDB)
{
    cblas_strsm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        (CBLAS_TRANSPOSE)transA, (CBLAS_DIAG)diag,
        M, N,
        (alpha), A, LDA,
        B, LDB);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_strsm(Quark *quark, Quark_Task_Flags *task_flags,
                      int side, int uplo, int transA, int diag,
                      int m, int n, int nb,
                      float alpha, float *A, int lda,
                      float *B, int ldb)
{
    DAG_CORE_TRSM;
    QUARK_Insert_Task(quark, CORE_strsm_quark, task_flags,
        sizeof(PLASMA_enum),                &side,      VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &diag,      VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(float),         &alpha,     VALUE,
        sizeof(float)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(float)*nb*nb,    B,                 INOUT | LOCALITY,
        sizeof(int),                        &ldb,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_strsm_quark = PCORE_strsm_quark
#define CORE_strsm_quark PCORE_strsm_quark
#endif
void CORE_strsm_quark(Quark *quark)
{
    int side;
    int uplo;
    int transA;
    int diag;
    int m;
    int n;
    float alpha;
    float *A;
    int lda;
    float *B;
    int ldb;

    quark_unpack_args_11(quark, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
    cblas_strsm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        (CBLAS_TRANSPOSE)transA, (CBLAS_DIAG)diag,
        m, n,
        (alpha), A, lda,
        B, ldb);
}
