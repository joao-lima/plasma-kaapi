/**
 *
 * @file core_clauum.c
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
 * @generated c Thu Sep 15 12:08:58 2011
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
#pragma weak CORE_clauum = PCORE_clauum
#define CORE_clauum PCORE_clauum
#endif
void CORE_clauum(int uplo, int N, PLASMA_Complex32_t *A, int LDA)
{
    int info;
    info = LAPACKE_clauum_work(LAPACK_COL_MAJOR, lapack_const(uplo), N, A, LDA );
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_clauum(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int n, int nb,
                       PLASMA_Complex32_t *A, int lda)
{
    DAG_CORE_LAUUM;
    QUARK_Insert_Task(quark, CORE_clauum_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A,             INOUT,
        sizeof(int),                        &lda,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_clauum_quark = PCORE_clauum_quark
#define CORE_clauum_quark PCORE_clauum_quark
#endif
void CORE_clauum_quark(Quark *quark)
{
    int uplo;
    int N;
    PLASMA_Complex32_t *A;
    int LDA;
    int info;

    quark_unpack_args_4(quark, uplo, N, A, LDA);
    info = LAPACKE_clauum_work(LAPACK_COL_MAJOR, lapack_const(uplo), N, A, LDA);
}
