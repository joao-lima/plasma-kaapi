/**
 *
 * @file core_zlauum.c
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
#pragma weak CORE_zlauum = PCORE_zlauum
#define CORE_zlauum PCORE_zlauum
#endif
void CORE_zlauum(int uplo, int N, PLASMA_Complex64_t *A, int LDA)
{
    int info;
    info = LAPACKE_zlauum_work(LAPACK_COL_MAJOR, lapack_const(uplo), N, A, LDA );
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zlauum(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int n, int nb,
                       PLASMA_Complex64_t *A, int lda)
{
    DAG_CORE_LAUUM;
    QUARK_Insert_Task(quark, CORE_zlauum_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A,             INOUT,
        sizeof(int),                        &lda,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlauum_quark = PCORE_zlauum_quark
#define CORE_zlauum_quark PCORE_zlauum_quark
#endif
void CORE_zlauum_quark(Quark *quark)
{
    int uplo;
    int N;
    PLASMA_Complex64_t *A;
    int LDA;
    int info;

    quark_unpack_args_4(quark, uplo, N, A, LDA);
    info = LAPACKE_zlauum_work(LAPACK_COL_MAJOR, lapack_const(uplo), N, A, LDA);
}
