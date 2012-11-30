/**
 *
 * @file core_dlauum.c
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
 * @generated d Thu Sep 15 12:08:58 2011
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
#pragma weak CORE_dlauum = PCORE_dlauum
#define CORE_dlauum PCORE_dlauum
#endif
void CORE_dlauum(int uplo, int N, double *A, int LDA)
{
    int info;
    info = LAPACKE_dlauum_work(LAPACK_COL_MAJOR, lapack_const(uplo), N, A, LDA );
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dlauum(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int n, int nb,
                       double *A, int lda)
{
    DAG_CORE_LAUUM;
    QUARK_Insert_Task(quark, CORE_dlauum_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(double)*nb*nb,    A,             INOUT,
        sizeof(int),                        &lda,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlauum_quark = PCORE_dlauum_quark
#define CORE_dlauum_quark PCORE_dlauum_quark
#endif
void CORE_dlauum_quark(Quark *quark)
{
    int uplo;
    int N;
    double *A;
    int LDA;
    int info;

    quark_unpack_args_4(quark, uplo, N, A, LDA);
    info = LAPACKE_dlauum_work(LAPACK_COL_MAJOR, lapack_const(uplo), N, A, LDA);
}
