/**
 *
 * @file core_strtri.c
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
#pragma weak CORE_strtri = PCORE_strtri
#define CORE_strtri PCORE_strtri
#endif
void CORE_strtri(int uplo, int diag, int N, float *A, int LDA, int *info)
{
    *info = LAPACKE_strtri_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo), lapack_const(diag),
        N, A, LDA);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_strtri(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int diag,
                       int n, int nb,
                       float *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo)
{
    QUARK_Insert_Task(
        quark, CORE_strtri_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &diag,      VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(float)*nb*nb,    A,                 INOUT,
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
#pragma weak CORE_strtri_quark = PCORE_strtri_quark
#define CORE_strtri_quark PCORE_strtri_quark
#endif
void CORE_strtri_quark(Quark *quark)
{
    int uplo;
    int diag;
    int N;
    float *A;
    int LDA;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    int iinfo;

    int info;

    quark_unpack_args_8(quark, uplo, diag, N, A, LDA, sequence, request, iinfo);
    info = LAPACKE_strtri_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo), lapack_const(diag),
        N, A, LDA);
    if ((sequence->status == PLASMA_SUCCESS) && (info != 0))
        plasma_sequence_flush(quark, sequence, request, iinfo + info);
}
