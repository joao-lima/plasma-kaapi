/**
 *
 * @file core_clanhe.c
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
 * @generated c Thu Sep 15 12:09:00 2011
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
#pragma weak CORE_clanhe = PCORE_clanhe
#define CORE_clanhe PCORE_clanhe
#endif
void CORE_clanhe(int norm, int uplo, int N,
                 PLASMA_Complex32_t *A, int LDA,
                 float *work, float *normA)
{
    *normA = LAPACKE_clanhe_work(
        LAPACK_COL_MAJOR,
        lapack_const(norm), lapack_const(uplo),
        N, A, LDA, work);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_clanhe(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, int uplo, int N,
                       PLASMA_Complex32_t *A, int LDA, int szeA,
                       int szeW, float *result)
{
    szeW = max(1, szeW);
    DAG_CORE_LANHE;
    QUARK_Insert_Task(quark, CORE_clanhe_quark, task_flags,
        sizeof(PLASMA_enum),                &norm,  VALUE,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(int),                        &N,     VALUE,
        sizeof(PLASMA_Complex32_t)*szeA,     A,             INPUT,
        sizeof(int),                        &LDA,   VALUE,
        sizeof(float)*szeW,                 NULL,          SCRATCH,
        sizeof(float),                     result,         OUTPUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_clanhe_quark = PCORE_clanhe_quark
#define CORE_clanhe_quark PCORE_clanhe_quark
#endif
void CORE_clanhe_quark(Quark *quark)
{
    float *normA;
    int norm;
    int uplo;
    int N;
    PLASMA_Complex32_t *A;
    int LDA;
    float *work;

    quark_unpack_args_7(quark, normA, norm, uplo, N, A, LDA, work);
    *normA = LAPACKE_clanhe_work(
        LAPACK_COL_MAJOR,
        lapack_const(norm), lapack_const(uplo),
        N, A, LDA, work);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_clanhe_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum norm, PLASMA_enum uplo, int N,
                          PLASMA_Complex32_t *A, int LDA, int szeA,
                          int szeW, float *result,
                          float *fake, int szeF)
{
    szeW = max(1, szeW);
    DAG_CORE_LANHE;
    QUARK_Insert_Task(quark, CORE_clanhe_f1_quark, task_flags,
        sizeof(PLASMA_enum),                &norm, VALUE,
        sizeof(PLASMA_enum),                &uplo, VALUE,
        sizeof(int),                        &N,    VALUE,
        sizeof(PLASMA_Complex32_t)*szeA,     A,             INPUT,
        sizeof(int),                        &LDA,  VALUE,
        sizeof(float)*szeW,                 NULL,          SCRATCH,
        sizeof(float),                      result,        OUTPUT,
        sizeof(float)*szeF,                 fake,          OUTPUT | GATHERV,
        0);
}
/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_clanhe_f1_quark = PCORE_clanhe_f1_quark
#define CORE_clanhe_f1_quark PCORE_clanhe_f1_quark
#endif
void CORE_clanhe_f1_quark(Quark *quark)
{
    float *normA;
    int norm;
    int uplo;
    int N;
    PLASMA_Complex32_t *A;
    int LDA;
    float *work;
    float *fake;

    quark_unpack_args_8(quark, norm, uplo, N, A, LDA, work, normA, fake);
    *normA = LAPACKE_clanhe_work(
        LAPACK_COL_MAJOR,
        lapack_const(norm), lapack_const(uplo),
        N, A, LDA, work);
}

