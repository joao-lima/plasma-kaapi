/**
 *
 * @file core_dtrmm.c
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
 * @generated d Thu Sep 15 12:08:59 2011
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dtrmm = PCORE_dtrmm
#define CORE_dtrmm PCORE_dtrmm
#endif
void CORE_dtrmm(int side, int uplo,
                int transA, int diag,
                int M, int N, 
                double alpha, 
                double *A, int LDA,
                double *B, int LDB)
{
    cblas_dtrmm(
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
void QUARK_CORE_dtrmm(Quark *quark, Quark_Task_Flags *task_flags,
                      int side, int uplo, int transA, int diag,
                      int m, int n, int nb,
                      double alpha, double *A, int lda,
                      double *B, int ldb)
{
    DAG_CORE_TRMM;
    QUARK_Insert_Task(quark, CORE_dtrmm_quark, task_flags,
        sizeof(PLASMA_enum),                &side,      VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &diag,      VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(double),         &alpha,     VALUE,
        sizeof(double)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(double)*nb*nb,    B,                 INOUT,
        sizeof(int),                        &ldb,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dtrmm_quark = PCORE_dtrmm_quark
#define CORE_dtrmm_quark PCORE_dtrmm_quark
#endif
void CORE_dtrmm_quark(Quark *quark)
{
    int side;
    int uplo;
    int transA;
    int diag;
    int M;
    int N;
    double alpha;
    double *A;
    int LDA;
    double *B;
    int LDB;

    quark_unpack_args_11(quark, side, uplo, transA, diag, M, N, alpha, A, LDA, B, LDB);
    cblas_dtrmm(
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
void QUARK_CORE_dtrmm_p2(Quark *quark, Quark_Task_Flags *task_flags,
                         int side, int uplo, int transA, int diag,
                         int m, int n, int nb,
                         double alpha, double *A, int lda,
                         double **B, int ldb)
{
    DAG_CORE_TRMM;
    QUARK_Insert_Task(quark, CORE_dtrmm_p2_quark, task_flags,
        sizeof(PLASMA_enum),                &side,      VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &diag,      VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(double),         &alpha,     VALUE,
        sizeof(double)*lda*nb,   A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(double*),         B,                 INOUT,
        sizeof(int),                        &ldb,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dtrmm_p2_quark = PCORE_dtrmm_p2_quark
#define CORE_dtrmm_p2_quark PCORE_dtrmm_p2_quark
#endif
void CORE_dtrmm_p2_quark(Quark *quark)
{
    int side;
    int uplo;
    int transA;
    int diag;
    int M;
    int N;
    double alpha;
    double *A;
    int LDA;
    double **B;
    int LDB;

    quark_unpack_args_11(quark, side, uplo, transA, diag, M, N, alpha, A, LDA, B, LDB);
    cblas_dtrmm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        (CBLAS_TRANSPOSE)transA, (CBLAS_DIAG)diag,
        M, N,
        (alpha), A, LDA,
        *B, LDB);
}
