/**
 *
 * @file core_dgetrf.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Thu Sep 15 12:09:01 2011
 *
 **/
#include <lapacke.h>
#include "common.h"

#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dgetrf = PCORE_dgetrf
#define CORE_dgetrf PCORE_dgetrf
#endif
int CORE_dgetrf(int m, int n,
                 double *A, int lda,
                 int *IPIV, int *info)
{
    *info = LAPACKE_dgetrf_work(LAPACK_COL_MAJOR, m, n, A, lda, IPIV );
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dgetrf(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int nb,
                       double *A, int lda,
                       int *IPIV,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       PLASMA_bool check_info, int iinfo)
{
    DAG_CORE_GETRF;
    QUARK_Insert_Task(quark, CORE_dgetrf_quark, task_flags,
        sizeof(int),                        &m,             VALUE,
        sizeof(int),                        &n,             VALUE,
        sizeof(double)*nb*nb,    A,                     INOUT | LOCALITY,
        sizeof(int),                        &lda,           VALUE,
        sizeof(int)*nb,                      IPIV,                  OUTPUT,
        sizeof(PLASMA_sequence*),           &sequence,      VALUE,
        sizeof(PLASMA_request*),            &request,       VALUE,
        sizeof(PLASMA_bool),                &check_info,    VALUE,
        sizeof(int),                        &iinfo,         VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dgetrf_quark = PCORE_dgetrf_quark
#define CORE_dgetrf_quark PCORE_dgetrf_quark
#endif
void CORE_dgetrf_quark(Quark *quark)
{
    int m;
    int n;
    double *A;
    int lda;
    int *IPIV;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    PLASMA_bool check_info;
    int iinfo;
    int info;

    quark_unpack_args_9(quark, m, n, A, lda, IPIV, sequence, request, check_info, iinfo);
    info = LAPACKE_dgetrf_work(LAPACK_COL_MAJOR, m, n, A, lda, IPIV );
    if (info != PLASMA_SUCCESS && check_info)
        plasma_sequence_flush(quark, sequence, request, iinfo+info);
}
