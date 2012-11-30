/**
 *
 * @file core_dsymm.c
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
#pragma weak CORE_dsymm = PCORE_dsymm
#define CORE_dsymm PCORE_dsymm
#endif
void CORE_dsymm(int side, int uplo,
                int M, int N,
                double alpha, double *A, int LDA,
                double *B, int LDB,
                double beta, double *C, int LDC)
{
    cblas_dsymm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        M, N,
        (alpha), A, LDA,
        B, LDB,
        (beta), C, LDC);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dsymm(Quark *quark, Quark_Task_Flags *task_flags,
                      int side, int uplo,
                      int m, int n, int nb,
                      double alpha, double *A, int lda,
                      double *B, int ldb,
                      double beta, double *C, int ldc)
{
    DAG_CORE_SYMM;
    QUARK_Insert_Task(quark, CORE_dsymm_quark, task_flags,
        sizeof(PLASMA_enum),                &side,    VALUE,
        sizeof(PLASMA_enum),                &uplo,    VALUE,
        sizeof(int),                        &m,       VALUE,
        sizeof(int),                        &n,       VALUE,
        sizeof(double),         &alpha,   VALUE,
        sizeof(double)*nb*nb,    A,               INPUT,
        sizeof(int),                        &lda,     VALUE,
        sizeof(double)*nb*nb,    B,               INPUT,
        sizeof(int),                        &ldb,     VALUE,
        sizeof(double),         &beta,    VALUE,
        sizeof(double)*nb*nb,    C,               INOUT,
        sizeof(int),                        &ldc,     VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dsymm_quark = PCORE_dsymm_quark
#define CORE_dsymm_quark PCORE_dsymm_quark
#endif
void CORE_dsymm_quark(Quark *quark)
{
    int side;
    int uplo;
    int M;
    int N;
    double alpha;
    double *A;
    int LDA;
    double *B;
    int LDB;
    double beta;
    double *C;
    int LDC;

    quark_unpack_args_12(quark, side, uplo, M, N, alpha, A, LDA, B, LDB, beta, C, LDC);
    cblas_dsymm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        M, N,
        (alpha), A, LDA,
        B, LDB,
        (beta), C, LDC);
}
