/**
 *
 * @file core_clacpy.c
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
 * @generated c Thu Sep 15 12:08:59 2011
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex32_t
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_clacpy = PCORE_clacpy
#define CORE_clacpy PCORE_clacpy
#endif
void CORE_clacpy(PLASMA_enum uplo, int M, int N,
                 PLASMA_Complex32_t *A, int LDA,
                 PLASMA_Complex32_t *B, int LDB)
{
    LAPACKE_clacpy_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        M, N, A, LDA, B, LDB);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_clacpy(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int m, int n, int nb,
                       PLASMA_Complex32_t *A, int lda,
                       PLASMA_Complex32_t *B, int ldb)
{
    DAG_CORE_LACPY;
    QUARK_Insert_Task(quark, CORE_clacpy_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A,             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    B,             OUTPUT,
        sizeof(int),                        &ldb,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_clacpy_quark = PCORE_clacpy_quark
#define CORE_clacpy_quark PCORE_clacpy_quark
#endif
void CORE_clacpy_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int M;
    int N;
    PLASMA_Complex32_t *A;
    int LDA;
    PLASMA_Complex32_t *B;
    int LDB;

    quark_unpack_args_7(quark, uplo, M, N, A, LDA, B, LDB);
    LAPACKE_clacpy_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        M, N, A, LDA, B, LDB);
}

