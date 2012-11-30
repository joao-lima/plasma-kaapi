/**
 *
 * @file core_sabsum.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Thu Sep 15 12:08:59 2011
 *
 **/
#include <cblas.h>
#include <math.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_float
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sasum = PCORE_sasum
#define CORE_sasum PCORE_sasum
#endif
void CORE_sasum(int storev, int uplo, int M, int N,
                 float *A, int lda, float *work)
{
    float *tmpA;
    float *tmpW, sum, abs;
    int i,j;

    switch (uplo) {
    case PlasmaUpper:
        for (j = 0; j < N; j++) {
            tmpA = A+(j*lda);
            sum = 0.0;
            for (i = 0; i < j; i++) {
                abs      = fabsf(*tmpA);
                sum     += abs;
                work[i] += abs;
                tmpA++;
            }
            work[j] += sum + fabsf(*tmpA);
        }
        break;
    case PlasmaLower:
        for (j = 0; j < N; j++) {
            tmpA = A+(j*lda)+j;

            sum = 0.0;
            work[j] += fabsf(*tmpA);

            tmpA++;
            for (i = j+1; i < M; i++) {
                abs      = fabsf(*tmpA);
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
                /* work[j] += cblas_sasum(M, &(A[j*lda]), 1); */
                tmpA = A+(j*lda);
                for (i = 0; i < M; i++) {
                    work[j] +=  fabsf(*tmpA);
                    tmpA++;
                }
            }
        }
        else {
            for (j = 0; j < N; j++) {
                tmpA = A+(j*lda);
                tmpW = work;
                for (i = 0; i < M; i++) {
                    /* work[i] += fabsf( A[j*lda+i] );*/
                    *tmpW += fabsf( *tmpA );
                    tmpA++; tmpW++;
                }
            }
        }
    }
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sasum(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum storev, PLASMA_enum uplo, int M, int N,
                       float *A, int lda, int szeA,
                       float *work, int szeW)
{
    QUARK_Insert_Task(
        quark, CORE_sasum_quark, task_flags,
        sizeof(PLASMA_enum),                &storev,    VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(int),                        &M,         VALUE,
        sizeof(int),                        &N,         VALUE,
        sizeof(float)*szeA,     A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(float)*szeW,                 work,              INOUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sasum_quark = PCORE_sasum_quark
#define CORE_sasum_quark PCORE_sasum_quark
#endif
void CORE_sasum_quark(Quark *quark)
{
    int storev;
    int uplo;
    int M;
    int N;
    float *A;
    int lda;
    float *work;

    quark_unpack_args_7(quark, storev, uplo, M, N, A, lda, work);
    CORE_sasum(storev, uplo, M, N, A, lda, work);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sasum_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum storev, PLASMA_enum uplo, int M, int N,
                          float *A, int lda, int szeA,
                          float *work, int szeW, float *fake, int szeF)
{
    DAG_CORE_ASUM;
    QUARK_Insert_Task(
        quark, CORE_sasum_f1_quark, task_flags,
        sizeof(PLASMA_enum),                &storev,    VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(int),                        &M,         VALUE,
        sizeof(int),                        &N,         VALUE,
        sizeof(float)*szeA,     A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(float)*szeW,                 work,              INOUT,
        sizeof(float)*szeF,                 fake,              OUTPUT | GATHERV,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sasum_f1_quark = PCORE_sasum_f1_quark
#define CORE_sasum_f1_quark PCORE_sasum_f1_quark
#endif
void CORE_sasum_f1_quark(Quark *quark)
{
    int storev;
    int uplo;
    int M;
    int N;
    float *A;
    int lda;
    float *work;
    float *fake;

    quark_unpack_args_8(quark, storev, uplo, M, N, A, lda, work, fake);
    CORE_sasum(storev, uplo, M, N, A, lda, work);
}
