/**
 *
 * @file core_dsyr2k.c
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
#pragma weak CORE_dsyr2k = PCORE_dsyr2k
#define CORE_dsyr2k PCORE_dsyr2k
#endif
void CORE_dsyr2k(int uplo, int trans,
                 int N, int K,
                 double alpha, double *A, int LDA,
                 double *B, int LDB,
                 double beta, double *C, int LDC)
{
    cblas_dsyr2k(
        CblasColMajor,
        (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
        N, K,
        (alpha), A, LDA, B, LDB,
        (beta), C, LDC);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dsyr2k(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int trans,
                       int n, int k, int nb,
                       double alpha, double *A, int lda,
                       double *B, int ldb,
                       double beta, double *C, int ldc)
{
    DAG_CORE_SYR2K;
    QUARK_Insert_Task(quark, CORE_dsyr2k_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &trans,     VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(double),         &alpha,     VALUE,
        sizeof(double)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(double)*nb*nb,    B,                 INPUT,
        sizeof(int),                        &ldb,       VALUE,
        sizeof(double),         &beta,      VALUE,
        sizeof(double)*nb*nb,    C,                 INOUT,
        sizeof(int),                        &ldc,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dsyr2k_quark = PCORE_dsyr2k_quark
#define CORE_dsyr2k_quark PCORE_dsyr2k_quark
#endif
void CORE_dsyr2k_quark(Quark *quark)
{
    int uplo;
    int trans;
    int n;
    int k;
    double alpha;
    double *A;
    int lda;
    double *B;
    int ldb;
    double beta;
    double *C;
    int ldc;

    quark_unpack_args_12(quark, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    CORE_dsyr2k(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
