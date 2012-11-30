/**
 *
 * @file core_dabsum.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Thu Sep 15 12:08:59 2011
 *
 **/
#include <cblas.h>
#include <math.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dasum = PCORE_dasum
#define CORE_dasum PCORE_dasum
#endif
void CORE_dasum(int storev, int uplo, int M, int N,
                 double *A, int lda, double *work)
{
    double *tmpA;
    double *tmpW, sum, abs;
    int i,j;

    switch (uplo) {
    case PlasmaUpper:
        for (j = 0; j < N; j++) {
            tmpA = A+(j*lda);
            sum = 0.0;
            for (i = 0; i < j; i++) {
                abs      = fabs(*tmpA);
                sum     += abs;
                work[i] += abs;
                tmpA++;
            }
            work[j] += sum + fabs(*tmpA);
        }
        break;
    case PlasmaLower:
        for (j = 0; j < N; j++) {
            tmpA = A+(j*lda)+j;

            sum = 0.0;
            work[j] += fabs(*tmpA);

            tmpA++;
            for (i = j+1; i < M; i++) {
                abs      = fabs(*tmpA);
                sum     += abs;
                work[i] += abs;
                tmpA++;
            }
            work[j] += sum;
        }
        break;
    case PlasmaUpperLower:
    default:
        if (storev == PlasmaColumnwise) {
            for (j = 0; j < N; j++) {
                /* work[j] += cblas_dasum(M, &(A[j*lda]), 1); */
                tmpA = A+(j*lda);
                for (i = 0; i < M; i++) {
                    work[j] +=  fabs(*tmpA);
                    tmpA++;
                }
            }
        }
        else {
            for (j = 0; j < N; j++) {
                tmpA = A+(j*lda);
                tmpW = work;
                for (i = 0; i < M; i++) {
                    /* work[i] += fabs( A[j*lda+i] );*/
                    *tmpW += fabs( *tmpA );
                    tmpA++; tmpW++;
                }
            }
        }
    }
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dasum(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum storev, PLASMA_enum uplo, int M, int N,
                       double *A, int lda, int szeA,
                       double *work, int szeW)
{
    QUARK_Insert_Task(
        quark, CORE_dasum_quark, task_flags,
        sizeof(PLASMA_enum),                &storev,    VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(int),                        &M,         VALUE,
        sizeof(int),                        &N,         VALUE,
        sizeof(double)*szeA,     A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(double)*szeW,                 work,              INOUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dasum_quark = PCORE_dasum_quark
#define CORE_dasum_quark PCORE_dasum_quark
#endif
void CORE_dasum_quark(Quark *quark)
{
    int storev;
    int uplo;
    int M;
    int N;
    double *A;
    int lda;
    double *work;

    quark_unpack_args_7(quark, storev, uplo, M, N, A, lda, work);
    CORE_dasum(storev, uplo, M, N, A, lda, work);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dasum_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum storev, PLASMA_enum uplo, int M, int N,
                          double *A, int lda, int szeA,
                          double *work, int szeW, double *fake, int szeF)
{
    DAG_CORE_ASUM;
    QUARK_Insert_Task(
        quark, CORE_dasum_f1_quark, task_flags,
        sizeof(PLASMA_enum),                &storev,    VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(int),                        &M,         VALUE,
        sizeof(int),                        &N,         VALUE,
        sizeof(double)*szeA,     A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(double)*szeW,                 work,              INOUT,
        sizeof(double)*szeF,                 fake,              OUTPUT | GATHERV,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dasum_f1_quark = PCORE_dasum_f1_quark
#define CORE_dasum_f1_quark PCORE_dasum_f1_quark
#endif
void CORE_dasum_f1_quark(Quark *quark)
{
    int storev;
    int uplo;
    int M;
    int N;
    double *A;
    int lda;
    double *work;
    double *fake;

    quark_unpack_args_8(quark, storev, uplo, M, N, A, lda, work, fake);
    CORE_dasum(storev, uplo, M, N, A, lda, work);
}
