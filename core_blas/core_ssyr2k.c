/**
 *
 * @file core_ssyr2k.c
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
#pragma weak CORE_ssyr2k = PCORE_ssyr2k
#define CORE_ssyr2k PCORE_ssyr2k
#endif
void CORE_ssyr2k(int uplo, int trans,
                 int N, int K,
                 float alpha, float *A, int LDA,
                 float *B, int LDB,
                 float beta, float *C, int LDC)
{
    cblas_ssyr2k(
        CblasColMajor,
        (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
        N, K,
        (alpha), A, LDA, B, LDB,
        (beta), C, LDC);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_ssyr2k(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int trans,
                       int n, int k, int nb,
                       float alpha, float *A, int lda,
                       float *B, int ldb,
                       float beta, float *C, int ldc)
{
    DAG_CORE_SYR2K;
    QUARK_Insert_Task(quark, CORE_ssyr2k_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &trans,     VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(float),         &alpha,     VALUE,
        sizeof(float)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(float)*nb*nb,    B,                 INPUT,
        sizeof(int),                        &ldb,       VALUE,
        sizeof(float),         &beta,      VALUE,
        sizeof(float)*nb*nb,    C,                 INOUT,
        sizeof(int),                        &ldc,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ssyr2k_quark = PCORE_ssyr2k_quark
#define CORE_ssyr2k_quark PCORE_ssyr2k_quark
#endif
void CORE_ssyr2k_quark(Quark *quark)
{
    int uplo;
    int trans;
    int n;
    int k;
    float alpha;
    float *A;
    int lda;
    float *B;
    int ldb;
    float beta;
    float *C;
    int ldc;

    quark_unpack_args_12(quark, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    CORE_ssyr2k(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
