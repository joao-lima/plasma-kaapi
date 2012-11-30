/**
 *
 * @file core_slaset2.c
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

/***************************************************************************/
/**
 *
 * @ingroup CORE_float
 *
 *  CORE_slaset2 - Sets the elements of the matrix A to alpha.
 *  Not LAPACK compliant! Read below.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies which elements of the matrix are to be set 
 *          = PlasmaUpper: STRICT Upper part of A is set to alpha;
 *          = PlasmaLower: STRICT Lower part of A is set to alpha;
 *          = PlasmaUpperLower: ALL elements of A are set to alpha.
 *          Not LAPACK Compliant.
 *
 * @param[in] M
 *          The number of rows of the matrix A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the matrix A.  N >= 0.
 *
 * @param[in] alpha
 *         The constant to which the elements are to be set.
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile A.
 *         On exit, A has been set to alpha accordingly.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slaset2 = PCORE_slaset2
#define CORE_slaset2 PCORE_slaset2
#endif
void CORE_slaset2(PLASMA_enum uplo, int M, int N,
                  float alpha, float *A, int LDA)
{
    if (uplo == PlasmaUpper) {
        LAPACKE_slaset_work(
            LAPACK_COL_MAJOR,
            lapack_const(uplo),
            M, N-1, alpha, alpha, A+LDA, LDA);
    }
    else if (uplo == PlasmaLower) {
        LAPACKE_slaset_work(
            LAPACK_COL_MAJOR,
            lapack_const(uplo),
            M-1, N, alpha, alpha, A+1, LDA);
    }
    else {
        LAPACKE_slaset_work(
            LAPACK_COL_MAJOR,
            lapack_const(uplo),
            M, N, alpha, alpha, A, LDA);
    }
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_slaset2(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int M, int N,
                       float alpha, float *A, int LDA)
{
    QUARK_Insert_Task(quark, CORE_slaset2_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(int),                        &M,     VALUE,
        sizeof(int),                        &N,     VALUE,
        sizeof(float),         &alpha, VALUE,
        sizeof(float)*M*N,     A,      OUTPUT,
        sizeof(int),                        &LDA,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slaset2_quark = PCORE_slaset2_quark
#define CORE_slaset2_quark PCORE_slaset2_quark
#endif
void CORE_slaset2_quark(Quark *quark)
{
    int uplo;
    int M;
    int N;
    float alpha;
    float *A;
    int LDA;

    quark_unpack_args_6(quark, uplo, M, N, alpha, A, LDA);
    CORE_slaset2(uplo, M, N, alpha, A, LDA);
}
