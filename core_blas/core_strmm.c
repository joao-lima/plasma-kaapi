/**
 *
 * @file core_strmm.c
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
#pragma weak CORE_strmm = PCORE_strmm
#define CORE_strmm PCORE_strmm
#endif
void CORE_strmm(int side, int uplo,
                int transA, int diag,
                int M, int N, 
                float alpha, 
                float *A, int LDA,
                float *B, int LDB)
{
    cblas_strmm(
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
void QUARK_CORE_strmm(Quark *quark, Quark_Task_Flags *task_flags,
                      int side, int uplo, int transA, int diag,
                      int m, int n, int nb,
                      float alpha, float *A, int lda,
                      float *B, int ldb)
{
    DAG_CORE_TRMM;
    QUARK_Insert_Task(quark, CORE_strmm_quark, task_flags,
        sizeof(PLASMA_enum),                &side,      VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &diag,      VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(float),         &alpha,     VALUE,
        sizeof(float)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(float)*nb*nb,    B,                 INOUT,
        sizeof(int),                        &ldb,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_strmm_quark = PCORE_strmm_quark
#define CORE_strmm_quark PCORE_strmm_quark
#endif
void CORE_strmm_quark(Quark *quark)
{
    int side;
    int uplo;
    int transA;
    int diag;
    int M;
    int N;
    float alpha;
    float *A;
    int LDA;
    float *B;
    int LDB;

    quark_unpack_args_11(quark, side, uplo, transA, diag, M, N, alpha, A, LDA, B, LDB);
    cblas_strmm(
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
void QUARK_CORE_strmm_p2(Quark *quark, Quark_Task_Flags *task_flags,
                         int side, int uplo, int transA, int diag,
                         int m, int n, int nb,
                         float alpha, float *A, int lda,
                         float **B, int ldb)
{
    DAG_CORE_TRMM;
    QUARK_Insert_Task(quark, CORE_strmm_p2_quark, task_flags,
        sizeof(PLASMA_enum),                &side,      VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &diag,      VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(float),         &alpha,     VALUE,
        sizeof(float)*lda*nb,   A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(float*),         B,                 INOUT,
        sizeof(int),                        &ldb,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_strmm_p2_quark = PCORE_strmm_p2_quark
#define CORE_strmm_p2_quark PCORE_strmm_p2_quark
#endif
void CORE_strmm_p2_quark(Quark *quark)
{
    int side;
    int uplo;
    int transA;
    int diag;
    int M;
    int N;
    float alpha;
    float *A;
    int LDA;
    float **B;
    int LDB;

    quark_unpack_args_11(quark, side, uplo, transA, diag, M, N, alpha, A, LDA, B, LDB);
    cblas_strmm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        (CBLAS_TRANSPOSE)transA, (CBLAS_DIAG)diag,
        M, N,
        (alpha), A, LDA,
        *B, LDB);
}
