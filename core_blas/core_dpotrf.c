/**
 *
 * @file core_dpotrf.c
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
#pragma weak CORE_dpotrf = PCORE_dpotrf
#define CORE_dpotrf PCORE_dpotrf
#endif
void CORE_dpotrf(int uplo, int N, double *A, int LDA, int *INFO)
{
#if 1
    *INFO = LAPACKE_dpotrf_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        N, A, LDA );
#else
    *INFO = clapack_dpotrf(CblasColMajor, CblasLower, N, A, LDA);
#endif
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dpotrf(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int n, int nb,
                       double *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo)
{
    DAG_CORE_POTRF;
    QUARK_Insert_Task(quark, CORE_dpotrf_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(double)*nb*nb,    A,                 INOUT,
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
#pragma weak CORE_dpotrf_quark = PCORE_dpotrf_quark
#define CORE_dpotrf_quark PCORE_dpotrf_quark
#endif
void CORE_dpotrf_quark(Quark *quark)
{
    int uplo;
    int n;
    double *A;
    int lda;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    int iinfo;

    int info;

#if 0
  fprintf(stdout, "%s\n", __FUNCTION__);
  fflush(stdout);
#endif
    quark_unpack_args_7(quark, uplo, n, A, lda, sequence, request, iinfo);
#if 0
    info = LAPACKE_dpotrf_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        n, A, lda);
#else
    info = clapack_dpotrf(CblasColMajor, CblasLower, n, A, lda);
#endif
    if (sequence->status == PLASMA_SUCCESS && info != 0)
      plasma_sequence_flush(quark, sequence, request, iinfo+info);
}
