/**
 *
 * @file core_zlacpy.c
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
 * @precisions normal z -> c d s
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlacpy = PCORE_zlacpy
#define CORE_zlacpy PCORE_zlacpy
#endif
void CORE_zlacpy(PLASMA_enum uplo, int M, int N,
                 PLASMA_Complex64_t *A, int LDA,
                 PLASMA_Complex64_t *B, int LDB)
{
    LAPACKE_zlacpy_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        M, N, A, LDA, B, LDB);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zlacpy(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int m, int n, int nb,
                       PLASMA_Complex64_t *A, int lda,
                       PLASMA_Complex64_t *B, int ldb)
{
    DAG_CORE_LACPY;
    QUARK_Insert_Task(quark, CORE_zlacpy_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A,             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    B,             OUTPUT,
        sizeof(int),                        &ldb,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlacpy_quark = PCORE_zlacpy_quark
#define CORE_zlacpy_quark PCORE_zlacpy_quark
#endif
void CORE_zlacpy_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int M;
    int N;
    PLASMA_Complex64_t *A;
    int LDA;
    PLASMA_Complex64_t *B;
    int LDB;

    quark_unpack_args_7(quark, uplo, M, N, A, LDA, B, LDB);
    LAPACKE_zlacpy_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        M, N, A, LDA, B, LDB);
}

