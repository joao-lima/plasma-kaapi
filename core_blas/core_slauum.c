/**
 *
 * @file core_slauum.c
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
 * @generated s Thu Sep 15 12:08:58 2011
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
#pragma weak CORE_slauum = PCORE_slauum
#define CORE_slauum PCORE_slauum
#endif
void CORE_slauum(int uplo, int N, float *A, int LDA)
{
    int info;
    info = LAPACKE_slauum_work(LAPACK_COL_MAJOR, lapack_const(uplo), N, A, LDA );
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_slauum(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int n, int nb,
                       float *A, int lda)
{
    DAG_CORE_LAUUM;
    QUARK_Insert_Task(quark, CORE_slauum_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(float)*nb*nb,    A,             INOUT,
        sizeof(int),                        &lda,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slauum_quark = PCORE_slauum_quark
#define CORE_slauum_quark PCORE_slauum_quark
#endif
void CORE_slauum_quark(Quark *quark)
{
    int uplo;
    int N;
    float *A;
    int LDA;
    int info;

    quark_unpack_args_4(quark, uplo, N, A, LDA);
    info = LAPACKE_slauum_work(LAPACK_COL_MAJOR, lapack_const(uplo), N, A, LDA);
}
