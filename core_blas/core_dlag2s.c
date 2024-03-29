/**
 *
 * @file core_dlag2s.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated ds Thu Sep 15 12:09:00 2011
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlag2s = PCORE_dlag2s
#define CORE_dlag2s PCORE_dlag2s
#endif
void CORE_dlag2s(int m, int n,
                 double *A, int lda,
                 float *B, int ldb, int *info)
{
    *info = LAPACKE_dlag2s_work(LAPACK_COL_MAJOR, m, n, A, lda, B, ldb);
}
/***************************************************************************//**
 *
 **/
void QUARK_CORE_dlag2s(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int nb,
                       double *A, int lda,
                       float *B, int ldb,
                       PLASMA_sequence *sequence, PLASMA_request *request)
{
    DAG_CORE_LAG2C;
    QUARK_Insert_Task(quark, CORE_dlag2s_quark, task_flags,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(double)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(float)*nb*nb,    B,                 OUTPUT,
        sizeof(int),                        &ldb,       VALUE,
        sizeof(PLASMA_sequence*),           &sequence,  VALUE,
        sizeof(PLASMA_request*),            &request,   VALUE,
        0);
}


/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlag2s_quark = PCORE_dlag2s_quark
#define CORE_dlag2s_quark PCORE_dlag2s_quark
#endif
void CORE_dlag2s_quark(Quark *quark)
{
    int m;
    int n;
    double *A;
    int lda;
    float *B;
    int ldb;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    int info;

    quark_unpack_args_8(quark, m, n, A, lda, B, ldb, sequence, request);
    info = LAPACKE_dlag2s_work(LAPACK_COL_MAJOR, m, n, A, lda, B, ldb);
    if (sequence->status == PLASMA_SUCCESS && info != 0)
        plasma_sequence_flush(quark, sequence, request, info);
}

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slag2d = PCORE_slag2d
#define CORE_slag2d PCORE_slag2d
#endif
void CORE_slag2d(int m, int n,
                 float *A, int lda,
                 double *B, int ldb)
{
    int info;
    info = LAPACKE_slag2d_work(LAPACK_COL_MAJOR, m, n, A, lda, B, ldb);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_slag2d(Quark *quark, Quark_Task_Flags *task_flags,
                      int m, int n, int nb,
                      float *A, int lda,
                      double *B, int ldb)
{
    QUARK_Insert_Task(quark, CORE_slag2d_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(float)*nb*nb,    A,             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(double)*nb*nb,    B,             INOUT,
        sizeof(int),                        &ldb,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slag2d_quark = PCORE_slag2d_quark
#define CORE_slag2d_quark PCORE_slag2d_quark
#endif
void CORE_slag2d_quark(Quark *quark)
{
    int m;
    int n;
    float *A;
    int lda;
    double *B;
    int ldb;
    int info;

    quark_unpack_args_6(quark, m, n, A, lda, B, ldb);
    info = LAPACKE_slag2d_work(LAPACK_COL_MAJOR, m, n, A, lda, B, ldb);
}

