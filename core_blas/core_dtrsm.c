/**
 *
 * @file core_dtrsm.c
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
#pragma weak CORE_dtrsm = PCORE_dtrsm
#define CORE_dtrsm PCORE_dtrsm
#endif
void CORE_dtrsm(int side, int uplo,
                int transA, int diag,
                int M, int N,
                double alpha, double *A, int LDA,
                double *B, int LDB)
{
    cblas_dtrsm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        (CBLAS_TRANSPOSE)transA, (CBLAS_DIAG)diag,
        M, N,
        (alpha), A, LDA,
        B, LDB);
#if 1
  fprintf(stdout,"%s: a=%p b=%p m=%d n=%d\n", __FUNCTION__,
          A, B, M, N);
  fflush(stdout);
#endif
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dtrsm(Quark *quark, Quark_Task_Flags *task_flags,
                      int side, int uplo, int transA, int diag,
                      int m, int n, int nb,
                      double alpha, double *A, int lda,
                      double *B, int ldb)
{
    DAG_CORE_TRSM;
    QUARK_Insert_Task(quark, CORE_dtrsm_quark, task_flags,
        sizeof(PLASMA_enum),                &side,      VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &diag,      VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(double),         &alpha,     VALUE,
        sizeof(double)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(double)*nb*nb,    B,                 INOUT | LOCALITY,
        sizeof(int),                        &ldb,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dtrsm_quark = PCORE_dtrsm_quark
#define CORE_dtrsm_quark PCORE_dtrsm_quark
#endif
void CORE_dtrsm_quark(Quark *quark)
{
    int side;
    int uplo;
    int transA;
    int diag;
    int m;
    int n;
    double alpha;
    double *A;
    int lda;
    double *B;
    int ldb;

    quark_unpack_args_11(quark, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
    cblas_dtrsm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        (CBLAS_TRANSPOSE)transA, (CBLAS_DIAG)diag,
        m, n,
        (alpha), A, lda,
        B, ldb);
}
