/**
 *
 * @file core_cgetrf.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated c Thu Sep 15 12:09:01 2011
 *
 **/
#include <lapacke.h>
#include "common.h"

#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cgetrf = PCORE_cgetrf
#define CORE_cgetrf PCORE_cgetrf
#endif
int CORE_cgetrf(int m, int n,
                 PLASMA_Complex32_t *A, int lda,
                 int *IPIV, int *info)
{
    *info = LAPACKE_cgetrf_work(LAPACK_COL_MAJOR, m, n, A, lda, IPIV );
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_cgetrf(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int nb,
                       PLASMA_Complex32_t *A, int lda,
                       int *IPIV,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       PLASMA_bool check_info, int iinfo)
{
    DAG_CORE_GETRF;
    QUARK_Insert_Task(quark, CORE_cgetrf_quark, task_flags,
        sizeof(int),                        &m,             VALUE,
        sizeof(int),                        &n,             VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A,                     INOUT | LOCALITY,
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
#pragma weak CORE_cgetrf_quark = PCORE_cgetrf_quark
#define CORE_cgetrf_quark PCORE_cgetrf_quark
#endif
void CORE_cgetrf_quark(Quark *quark)
{
    int m;
    int n;
    PLASMA_Complex32_t *A;
    int lda;
    int *IPIV;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    PLASMA_bool check_info;
    int iinfo;
    int info;

    quark_unpack_args_9(quark, m, n, A, lda, IPIV, sequence, request, check_info, iinfo);
    info = LAPACKE_cgetrf_work(LAPACK_COL_MAJOR, m, n, A, lda, IPIV );
    if (info != PLASMA_SUCCESS && check_info)
        plasma_sequence_flush(quark, sequence, request, iinfo+info);
}
