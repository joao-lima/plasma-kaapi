/**
 *
 * @file core_slacpy.c
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
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_float
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slacpy = PCORE_slacpy
#define CORE_slacpy PCORE_slacpy
#endif
void CORE_slacpy(PLASMA_enum uplo, int M, int N,
                 float *A, int LDA,
                 float *B, int LDB)
{
    LAPACKE_slacpy_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        M, N, A, LDA, B, LDB);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_slacpy(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int m, int n, int nb,
                       float *A, int lda,
                       float *B, int ldb)
{
    DAG_CORE_LACPY;
    QUARK_Insert_Task(quark, CORE_slacpy_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(float)*nb*nb,    A,             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(float)*nb*nb,    B,             OUTPUT,
        sizeof(int),                        &ldb,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slacpy_quark = PCORE_slacpy_quark
#define CORE_slacpy_quark PCORE_slacpy_quark
#endif
void CORE_slacpy_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int M;
    int N;
    float *A;
    int LDA;
    float *B;
    int LDB;

    quark_unpack_args_7(quark, uplo, M, N, A, LDA, B, LDB);
    LAPACKE_slacpy_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        M, N, A, LDA, B, LDB);
}

