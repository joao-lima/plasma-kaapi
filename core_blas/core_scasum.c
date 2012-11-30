/**
 *
 * @file core_cabsum.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated c Thu Sep 15 12:08:59 2011
 *
 **/
#include <cblas.h>
#include <math.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex32_t
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_scasum = PCORE_scasum
#define CORE_scasum PCORE_scasum
#endif
void CORE_scasum(int storev, int uplo, int M, int N,
                 PLASMA_Complex32_t *A, int lda, float *work)
{
    PLASMA_Complex32_t *tmpA;
    float *tmpW, sum, abs;
    int i,j;

    switch (uplo) {
    case PlasmaUpper:
        for (j = 0; j < N; j++) {
            tmpA = A+(j*lda);
            sum = 0.0;
            for (i = 0; i < j; i++) {
                abs      = cabsf(*tmpA);
                sum     += abs;
                work[i] += abs;
                tmpA++;
            }
            work[j] += sum + cabsf(*tmpA);
        }
        break;
    case PlasmaLower:
        for (j = 0; j < N; j++) {
            tmpA = A+(j*lda)+j;

            sum = 0.0;
            work[j] += cabsf(*tmpA);

            tmpA++;
            for (i = j+1; i < M; i++) {
                abs      = cabsf(*tmpA);
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
                /* work[j] += cblas_scasum(M, &(A[j*lda]), 1); */
                tmpA = A+(j*lda);
                for (i = 0; i < M; i++) {
                    work[j] +=  cabsf(*tmpA);
                    tmpA++;
                }
            }
        }
        else {
            for (j = 0; j < N; j++) {
                tmpA = A+(j*lda);
                tmpW = work;
                for (i = 0; i < M; i++) {
                    /* work[i] += cabsf( A[j*lda+i] );*/
                    *tmpW += cabsf( *tmpA );
                    tmpA++; tmpW++;
                }
            }
        }
    }
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_scasum(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum storev, PLASMA_enum uplo, int M, int N,
                       PLASMA_Complex32_t *A, int lda, int szeA,
                       float *work, int szeW)
{
    QUARK_Insert_Task(
        quark, CORE_scasum_quark, task_flags,
        sizeof(PLASMA_enum),                &storev,    VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(int),                        &M,         VALUE,
        sizeof(int),                        &N,         VALUE,
        sizeof(PLASMA_Complex32_t)*szeA,     A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(float)*szeW,                 work,              INOUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_scasum_quark = PCORE_scasum_quark
#define CORE_scasum_quark PCORE_scasum_quark
#endif
void CORE_scasum_quark(Quark *quark)
{
    int storev;
    int uplo;
    int M;
    int N;
    PLASMA_Complex32_t *A;
    int lda;
    float *work;

    quark_unpack_args_7(quark, storev, uplo, M, N, A, lda, work);
    CORE_scasum(storev, uplo, M, N, A, lda, work);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_scasum_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum storev, PLASMA_enum uplo, int M, int N,
                          PLASMA_Complex32_t *A, int lda, int szeA,
                          float *work, int szeW, float *fake, int szeF)
{
    DAG_CORE_ASUM;
    QUARK_Insert_Task(
        quark, CORE_scasum_f1_quark, task_flags,
        sizeof(PLASMA_enum),                &storev,    VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(int),                        &M,         VALUE,
        sizeof(int),                        &N,         VALUE,
        sizeof(PLASMA_Complex32_t)*szeA,     A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(float)*szeW,                 work,              INOUT,
        sizeof(float)*szeF,                 fake,              OUTPUT | GATHERV,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_scasum_f1_quark = PCORE_scasum_f1_quark
#define CORE_scasum_f1_quark PCORE_scasum_f1_quark
#endif
void CORE_scasum_f1_quark(Quark *quark)
{
    int storev;
    int uplo;
    int M;
    int N;
    PLASMA_Complex32_t *A;
    int lda;
    float *work;
    float *fake;

    quark_unpack_args_8(quark, storev, uplo, M, N, A, lda, work, fake);
    CORE_scasum(storev, uplo, M, N, A, lda, work);
}
