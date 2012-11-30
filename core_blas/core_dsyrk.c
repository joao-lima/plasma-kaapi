/**
 *
 * @file core_dsyrk.c
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
 * @generated d Thu Sep 15 12:08:58 2011
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dsyrk = PCORE_dsyrk
#define CORE_dsyrk PCORE_dsyrk
#endif
void CORE_dsyrk(int uplo, int trans,
                int N, int K,
                double alpha, double *A, int LDA,
                double beta, double *C, int LDC)
{
    cblas_dsyrk(
        CblasColMajor,
        (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
        N, K,
        (alpha), A, LDA,
        (beta), C, LDC);
#if 1
  fprintf(stdout,"%s: a=%p c=%p n=%d k=%d\n", __FUNCTION__,
          A, C, N, K);
  fflush(stdout);
#endif
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dsyrk(Quark *quark, Quark_Task_Flags *task_flags,
                      int uplo, int trans,
                      int n, int k, int nb,
                      double alpha, double *A, int lda,
                      double beta, double *C, int ldc)
{
    DAG_CORE_SYRK;
    QUARK_Insert_Task(quark, CORE_dsyrk_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &trans,     VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(double),         &alpha,     VALUE,
        sizeof(double)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(double),         &beta,      VALUE,
        sizeof(double)*nb*nb,    C,                 INOUT,
        sizeof(int),                        &ldc,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dsyrk_quark = PCORE_dsyrk_quark
#define CORE_dsyrk_quark PCORE_dsyrk_quark
#endif
void CORE_dsyrk_quark(Quark *quark)
{
    int uplo;
    int trans;
    int n;
    int k;
    double alpha;
    double *A;
    int lda;
    double beta;
    double *C;
    int ldc;

    quark_unpack_args_10(quark, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
    cblas_dsyrk(
        CblasColMajor,
        (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
        n, k,
        (alpha), A, lda,
        (beta), C, ldc);
}
