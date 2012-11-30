/**
 *
 * @file core_cpotrf.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
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
#pragma weak CORE_cpotrf = PCORE_cpotrf
#define CORE_cpotrf PCORE_cpotrf
#endif
void CORE_cpotrf(int uplo, int N, PLASMA_Complex32_t *A, int LDA, int *INFO)
{
    *INFO = LAPACKE_cpotrf_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        N, A, LDA );
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_cpotrf(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int n, int nb,
                       PLASMA_Complex32_t *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo)
{
    DAG_CORE_POTRF;
    QUARK_Insert_Task(quark, CORE_cpotrf_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A,                 INOUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(PLASMA_sequence*),           &sequence,  VALUE,
        sizeof(PLASMA_request*),            &request,   VALUE,
        sizeof(int),                        &iinfo,     VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cpotrf_quark = PCORE_cpotrf_quark
#define CORE_cpotrf_quark PCORE_cpotrf_quark
#endif
void CORE_cpotrf_quark(Quark *quark)
{
    int uplo;
    int n;
    PLASMA_Complex32_t *A;
    int lda;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    int iinfo;

    int info;

    quark_unpack_args_7(quark, uplo, n, A, lda, sequence, request, iinfo);
    info = LAPACKE_cpotrf_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        n, A, lda);
    if (sequence->status == PLASMA_SUCCESS && info != 0)
      plasma_sequence_flush(quark, sequence, request, iinfo+info);
}
