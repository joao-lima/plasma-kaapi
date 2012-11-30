/**
 *
 * @file core_ctrtri.c
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
#pragma weak CORE_ctrtri = PCORE_ctrtri
#define CORE_ctrtri PCORE_ctrtri
#endif
void CORE_ctrtri(int uplo, int diag, int N, PLASMA_Complex32_t *A, int LDA, int *info)
{
    *info = LAPACKE_ctrtri_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo), lapack_const(diag),
        N, A, LDA);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_ctrtri(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int diag,
                       int n, int nb,
                       PLASMA_Complex32_t *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo)
{
    QUARK_Insert_Task(
        quark, CORE_ctrtri_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &diag,      VALUE,
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
#pragma weak CORE_ctrtri_quark = PCORE_ctrtri_quark
#define CORE_ctrtri_quark PCORE_ctrtri_quark
#endif
void CORE_ctrtri_quark(Quark *quark)
{
    int uplo;
    int diag;
    int N;
    PLASMA_Complex32_t *A;
    int LDA;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    int iinfo;

    int info;

    quark_unpack_args_8(quark, uplo, diag, N, A, LDA, sequence, request, iinfo);
    info = LAPACKE_ctrtri_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo), lapack_const(diag),
        N, A, LDA);
    if ((sequence->status == PLASMA_SUCCESS) && (info != 0))
        plasma_sequence_flush(quark, sequence, request, iinfo + info);
}
