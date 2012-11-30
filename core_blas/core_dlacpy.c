/**
 *
 * @file core_dlacpy.c
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
#pragma weak CORE_dlacpy = PCORE_dlacpy
#define CORE_dlacpy PCORE_dlacpy
#endif
void CORE_dlacpy(PLASMA_enum uplo, int M, int N,
                 double *A, int LDA,
                 double *B, int LDB)
{
    LAPACKE_dlacpy_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        M, N, A, LDA, B, LDB);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dlacpy(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int m, int n, int nb,
                       double *A, int lda,
                       double *B, int ldb)
{
    DAG_CORE_LACPY;
    QUARK_Insert_Task(quark, CORE_dlacpy_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(double)*nb*nb,    A,             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(double)*nb*nb,    B,             OUTPUT,
        sizeof(int),                        &ldb,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlacpy_quark = PCORE_dlacpy_quark
#define CORE_dlacpy_quark PCORE_dlacpy_quark
#endif
void CORE_dlacpy_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int M;
    int N;
    double *A;
    int LDA;
    double *B;
    int LDB;

    quark_unpack_args_7(quark, uplo, M, N, A, LDA, B, LDB);
    LAPACKE_dlacpy_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        M, N, A, LDA, B, LDB);
}

