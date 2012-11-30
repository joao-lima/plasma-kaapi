/**
 *
 * @file core_ssygst.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @generated s Thu Sep 15 12:09:01 2011
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
#pragma weak CORE_ssygst = PCORE_ssygst
#define CORE_ssygst PCORE_ssygst
#endif
void CORE_ssygst(int itype, PLASMA_enum uplo, int N, float *A, int LDA, float *B, int LDB, int *INFO)
{
    *INFO = LAPACKE_ssygst_work(
        LAPACK_COL_MAJOR,
        itype,
        lapack_const(uplo),
        N, A, LDA, B, LDB );
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_ssygst(Quark *quark, Quark_Task_Flags *task_flags,
                       int itype, PLASMA_enum uplo, int n,
                       float *A, int lda,
                       float *B, int ldb,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo)
{
    QUARK_Insert_Task(quark, CORE_ssygst_quark, task_flags,
        sizeof(int),                        &itype,      VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(float)*n*n,    A,                 INOUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(float)*n*n,    B,                 INOUT,
        sizeof(int),                        &ldb,       VALUE,
        sizeof(PLASMA_sequence*),           &sequence,  VALUE,
        sizeof(PLASMA_request*),            &request,   VALUE,
        sizeof(int),                        &iinfo,     VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ssygst_quark = PCORE_ssygst_quark
#define CORE_ssygst_quark PCORE_ssygst_quark
#endif
void CORE_ssygst_quark(Quark *quark)
{
    int itype;
    PLASMA_enum uplo;
    int n;
    float *A;
    int lda;
    float *B;
    int ldb;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    int iinfo;

    int info;

    quark_unpack_args_10(quark, itype, uplo, n, A, lda, B, ldb, sequence, request, iinfo);
    info = LAPACKE_ssygst_work(
        LAPACK_COL_MAJOR,
        itype,
        lapack_const(uplo),
        n, A, lda, B, ldb);
    if (sequence->status == PLASMA_SUCCESS && info != 0)
      plasma_sequence_flush(quark, sequence, request, iinfo+info);
}
