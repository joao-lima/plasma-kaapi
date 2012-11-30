/**
 *
 * @file core_ssymm.c
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
 * @generated s Thu Sep 15 12:08:58 2011
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_float
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ssymm = PCORE_ssymm
#define CORE_ssymm PCORE_ssymm
#endif
void CORE_ssymm(int side, int uplo,
                int M, int N,
                float alpha, float *A, int LDA,
                float *B, int LDB,
                float beta, float *C, int LDC)
{
    cblas_ssymm(
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
void QUARK_CORE_ssymm(Quark *quark, Quark_Task_Flags *task_flags,
                      int side, int uplo,
                      int m, int n, int nb,
                      float alpha, float *A, int lda,
                      float *B, int ldb,
                      float beta, float *C, int ldc)
{
    DAG_CORE_SYMM;
    QUARK_Insert_Task(quark, CORE_ssymm_quark, task_flags,
        sizeof(PLASMA_enum),                &side,    VALUE,
        sizeof(PLASMA_enum),                &uplo,    VALUE,
        sizeof(int),                        &m,       VALUE,
        sizeof(int),                        &n,       VALUE,
        sizeof(float),         &alpha,   VALUE,
        sizeof(float)*nb*nb,    A,               INPUT,
        sizeof(int),                        &lda,     VALUE,
        sizeof(float)*nb*nb,    B,               INPUT,
        sizeof(int),                        &ldb,     VALUE,
        sizeof(float),         &beta,    VALUE,
        sizeof(float)*nb*nb,    C,               INOUT,
        sizeof(int),                        &ldc,     VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ssymm_quark = PCORE_ssymm_quark
#define CORE_ssymm_quark PCORE_ssymm_quark
#endif
void CORE_ssymm_quark(Quark *quark)
{
    int side;
    int uplo;
    int M;
    int N;
    float alpha;
    float *A;
    int LDA;
    float *B;
    int LDB;
    float beta;
    float *C;
    int LDC;

    quark_unpack_args_12(quark, side, uplo, M, N, alpha, A, LDA, B, LDB, beta, C, LDC);
    cblas_ssymm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        M, N,
        (alpha), A, LDA,
        B, LDB,
        (beta), C, LDC);
}
