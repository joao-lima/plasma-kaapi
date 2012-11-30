/**
 *
 * @file core_dgemm.c
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
#pragma weak CORE_dgemm = PCORE_dgemm
#define CORE_dgemm PCORE_dgemm
#endif
void CORE_dgemm(int transA, int transB,
                int M, int N, int K,
                double alpha, double *A, int LDA,
                double *B, int LDB,
                double beta, double *C, int LDC)
{
    cblas_dgemm(
        CblasColMajor,
        (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
        M, N, K,
        (alpha), A, LDA,
        B, LDB,
        (beta), C, LDC);
#if 1
  fprintf(stdout,"%s: a=%p b=%p c=%p m=%d n=%d k=%d\n", __FUNCTION__,
          A, B, C, M, N, K);
  fflush(stdout);
#endif
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dgemm(Quark *quark, Quark_Task_Flags *task_flags,
                      int transA, int transB,
                      int m, int n, int k, int nb,
                      double alpha, double *A, int lda,
                      double *B, int ldb,
                      double beta, double *C, int ldc)
{
    DAG_CORE_GEMM;
    QUARK_Insert_Task(quark, CORE_dgemm_quark, task_flags,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &transB,    VALUE,
        sizeof(int),                        &m,         VALUE,
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
void QUARK_CORE_dgemm2( Quark *quark, Quark_Task_Flags *task_flags,
                        int transA, int transB,
                        int m, int n, int k, int nb,
                        double alpha, double *A, int lda,
                        double *B, int ldb,
                        double beta, double *C, int ldc)
{
    DAG_CORE_GEMM;
    QUARK_Insert_Task(quark, CORE_dgemm_quark, task_flags,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &transB,    VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(double),         &alpha,     VALUE,
        sizeof(double)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(double)*nb*nb,    B,                 INPUT,
        sizeof(int),                        &ldb,       VALUE,
        sizeof(double),         &beta,      VALUE,
        sizeof(double)*nb*nb,    C,                 INOUT | LOCALITY | GATHERV,
        sizeof(int),                        &ldc,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dgemm_quark = PCORE_dgemm_quark
#define CORE_dgemm_quark PCORE_dgemm_quark
#endif
void CORE_dgemm_quark(Quark *quark)
{
    int transA;
    int transB;
    int m;
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

    quark_unpack_args_13(quark, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    cblas_dgemm(
        CblasColMajor,
        (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
        m, n, k,
        (alpha), A, lda,
        B, ldb,
        (beta), C, ldc);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dgemm_f2(Quark *quark, Quark_Task_Flags *task_flags,
                         int transA, int transB,
                         int m, int n, int k, int nb,
                         double alpha, double *A, int lda,
                         double *B, int ldb,
                         double beta, double *C, int ldc,
                         double *fake1, int szefake1, int flag1,
                         double *fake2, int szefake2, int flag2)
{
    DAG_CORE_GEMM;
    QUARK_Insert_Task(quark, CORE_dgemm_f2_quark, task_flags,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &transB,    VALUE,
        sizeof(int),                        &m,         VALUE,
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
        sizeof(double)*szefake1, fake1,             flag1,
        sizeof(double)*szefake2, fake2,             flag2,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dgemm_f2_quark = PCORE_dgemm_f2_quark
#define CORE_dgemm_f2_quark PCORE_dgemm_f2_quark
#endif
void CORE_dgemm_f2_quark(Quark* quark)
{
    int transA;
    int transB;
    int M;
    int N;
    int K;
    double alpha;
    double *A;
    int LDA;
    double *B;
    int LDB;
    double beta;
    double *C;
    int LDC;
    void *fake1, *fake2;
    
    quark_unpack_args_15(quark, transA, transB, M, N, K, alpha, 
                         A, LDA, B, LDB, beta, C, LDC, fake1, fake2);
    cblas_dgemm(
        CblasColMajor,
        (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
        M, N, K,
        (alpha), A, LDA,
        B, LDB,
        (beta), C, LDC);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dgemm_p2(Quark *quark, Quark_Task_Flags *task_flags,
                           int transA, int transB,
                           int m, int n, int k, int nb,
                           double alpha, double *A, int lda,
                           double **B, int ldb,
                           double beta, double *C, int ldc)
{
    DAG_CORE_GEMM;
    QUARK_Insert_Task(quark, CORE_dgemm_p2_quark, task_flags,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &transB,    VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(double),         &alpha,     VALUE,
        sizeof(double)*lda*nb,   A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(double*),         B,                 INPUT,
        sizeof(int),                        &ldb,       VALUE,
        sizeof(double),         &beta,      VALUE,
        sizeof(double)*ldc*nb,    C,                 INOUT,
        sizeof(int),                        &ldc,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dgemm_p2_quark = PCORE_dgemm_p2_quark
#define CORE_dgemm_p2_quark PCORE_dgemm_p2_quark
#endif
void CORE_dgemm_p2_quark(Quark* quark)
{
    int transA;
    int transB;
    int M;
    int N;
    int K;
    double alpha;
    double *A;
    int LDA;
    double **B;
    int LDB;
    double beta;
    double *C;
    int LDC;
    
    quark_unpack_args_13(quark, transA, transB, M, N, K, alpha, 
                         A, LDA, B, LDB, beta, C, LDC);
    cblas_dgemm(
        CblasColMajor,
        (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
        M, N, K,
        (alpha), A, LDA,
        *B, LDB,
        (beta), C, LDC);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dgemm_p3(Quark *quark, Quark_Task_Flags *task_flags,
                           int transA, int transB,
                           int m, int n, int k, int nb,
                           double alpha, double *A, int lda,
                           double *B, int ldb,
                           double beta, double **C, int ldc)
{
    DAG_CORE_GEMM;
    QUARK_Insert_Task(quark, CORE_dgemm_p3_quark, task_flags,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &transB,    VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(double),         &alpha,     VALUE,
        sizeof(double)*lda*nb,   A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(double)*ldb*nb,   B,                 INPUT,
        sizeof(int),                        &ldb,       VALUE,
        sizeof(double),         &beta,      VALUE,
        sizeof(double*),         C,                 INOUT,
        sizeof(int),                        &ldc,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dgemm_p3_quark = PCORE_dgemm_p3_quark
#define CORE_dgemm_p3_quark PCORE_dgemm_p3_quark
#endif
void CORE_dgemm_p3_quark(Quark* quark)
{
    int transA;
    int transB;
    int M;
    int N;
    int K;
    double alpha;
    double *A;
    int LDA;
    double *B;
    int LDB;
    double beta;
    double **C;
    int LDC;
    
    quark_unpack_args_13(quark, transA, transB, M, N, K, alpha, 
                         A, LDA, B, LDB, beta, C, LDC);
    cblas_dgemm(
        CblasColMajor,
        (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
        M, N, K,
        (alpha), A, LDA,
        B, LDB,
        (beta), *C, LDC);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dgemm_p2f1(Quark *quark, Quark_Task_Flags *task_flags,
                           int transA, int transB,
                           int m, int n, int k, int nb,
                           double alpha, double *A, int lda,
                           double **B, int ldb,
                           double beta, double *C, int ldc,
                           double *fake1, int szefake1, int flag1)
{
    DAG_CORE_GEMM;
    QUARK_Insert_Task(quark, CORE_dgemm_p2f1_quark, task_flags,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &transB,    VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(double),         &alpha,     VALUE,
        sizeof(double)*lda*nb,   A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(double*),         B,                 INPUT,
        sizeof(int),                        &ldb,       VALUE,
        sizeof(double),         &beta,      VALUE,
        sizeof(double)*ldc*nb,    C,                 INOUT,
        sizeof(int),                        &ldc,       VALUE,
        sizeof(double)*szefake1, fake1,             flag1,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dgemm_p2f1_quark = PCORE_dgemm_p2f1_quark
#define CORE_dgemm_p2f1_quark PCORE_dgemm_p2f1_quark
#endif
void CORE_dgemm_p2f1_quark(Quark* quark)
{
    int transA;
    int transB;
    int M;
    int N;
    int K;
    double alpha;
    double *A;
    int LDA;
    double **B;
    int LDB;
    double beta;
    double *C;
    int LDC;
    void *fake1;
    
    quark_unpack_args_14(quark, transA, transB, M, N, K, alpha, 
                         A, LDA, B, LDB, beta, C, LDC, fake1);
    cblas_dgemm(
        CblasColMajor,
        (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
        M, N, K,
        (alpha), A, LDA,
        *B, LDB,
        (beta), C, LDC);
}
