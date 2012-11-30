/**
 *
 * @file core_dlansy.c
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
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlansy = PCORE_dlansy
#define CORE_dlansy PCORE_dlansy
#endif
void CORE_dlansy(int norm, int uplo, int N,
                 double *A, int LDA,
                 double *work, double *normA)
{
    *normA = LAPACKE_dlansy_work(
        LAPACK_COL_MAJOR,
        lapack_const(norm), lapack_const(uplo),
        N, A, LDA, work);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dlansy(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, int uplo, int N,
                       double *A, int LDA, int szeA,
                       int szeW, double *result)
{
    szeW = max(1, szeW);
    DAG_CORE_LANSY;
    QUARK_Insert_Task(quark, CORE_dlansy_quark, task_flags,
        sizeof(PLASMA_enum),                &norm,  VALUE,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(int),                        &N,     VALUE,
        sizeof(double)*szeA,     A,             INPUT,
        sizeof(int),                        &LDA,   VALUE,
        sizeof(double)*szeW,                 NULL,          SCRATCH,
        sizeof(double),                      result,        OUTPUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlansy_quark = PCORE_dlansy_quark
#define CORE_dlansy_quark PCORE_dlansy_quark
#endif
void CORE_dlansy_quark(Quark *quark)
{
    double *normA;
    int norm;
    int uplo;
    int N;
    double *A;
    int LDA;
    double *work;

    quark_unpack_args_7(quark, normA, norm, uplo, N, A, LDA, work);
    *normA = LAPACKE_dlansy_work(
        LAPACK_COL_MAJOR,
        lapack_const(norm), lapack_const(uplo),
        N, A, LDA, work);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dlansy_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum norm, PLASMA_enum uplo, int N,
                          double *A, int LDA, int szeA,
                          int szeW, double *result,
                          double *fake, int szeF)
{
    szeW = max(1, szeW);
    DAG_CORE_LANSY;
    QUARK_Insert_Task(quark, CORE_dlansy_f1_quark, task_flags,
        sizeof(PLASMA_enum),                &norm,  VALUE,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(int),                        &N,     VALUE,
        sizeof(double)*szeA,     A,             INPUT,
        sizeof(int),                        &LDA,   VALUE,
        sizeof(double)*szeW,                 NULL,          SCRATCH,
        sizeof(double),                      result,        OUTPUT,
        sizeof(double)*szeF,                 fake,          OUTPUT | GATHERV,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlansy_f1_quark = PCORE_dlansy_f1_quark
#define CORE_dlansy_f1_quark PCORE_dlansy_f1_quark
#endif
void CORE_dlansy_f1_quark(Quark *quark)
{
    double *normA;
    int norm;
    int uplo;
    int N;
    double *A;
    int LDA;
    double *work;
    double *fake;

    quark_unpack_args_8(quark, norm, uplo, N, A, LDA, work, normA, fake);
    *normA = LAPACKE_dlansy_work(
        LAPACK_COL_MAJOR,
        lapack_const(norm),
        lapack_const(uplo),
        N, A, LDA, work);
}
