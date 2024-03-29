/**
 *
 * @file core_claset2.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @generated c Thu Sep 15 12:09:01 2011
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************/
/**
 *
 * @ingroup CORE_PLASMA_Complex32_t
 *
 *  CORE_claset2 - Sets the elements of the matrix A to alpha.
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
#pragma weak CORE_claset2 = PCORE_claset2
#define CORE_claset2 PCORE_claset2
#endif
void CORE_claset2(PLASMA_enum uplo, int M, int N,
                  PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int LDA)
{
    if (uplo == PlasmaUpper) {
        LAPACKE_claset_work(
            LAPACK_COL_MAJOR,
            lapack_const(uplo),
            M, N-1, alpha, alpha, A+LDA, LDA);
    }
    else if (uplo == PlasmaLower) {
        LAPACKE_claset_work(
            LAPACK_COL_MAJOR,
            lapack_const(uplo),
            M-1, N, alpha, alpha, A+1, LDA);
    }
    else {
        LAPACKE_claset_work(
            LAPACK_COL_MAJOR,
            lapack_const(uplo),
            M, N, alpha, alpha, A, LDA);
    }
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_claset2(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int M, int N,
                       PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int LDA)
{
    QUARK_Insert_Task(quark, CORE_claset2_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(int),                        &M,     VALUE,
        sizeof(int),                        &N,     VALUE,
        sizeof(PLASMA_Complex32_t),         &alpha, VALUE,
        sizeof(PLASMA_Complex32_t)*M*N,     A,      OUTPUT,
        sizeof(int),                        &LDA,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_claset2_quark = PCORE_claset2_quark
#define CORE_claset2_quark PCORE_claset2_quark
#endif
void CORE_claset2_quark(Quark *quark)
{
    int uplo;
    int M;
    int N;
    PLASMA_Complex32_t alpha;
    PLASMA_Complex32_t *A;
    int LDA;

    quark_unpack_args_6(quark, uplo, M, N, alpha, A, LDA);
    CORE_claset2(uplo, M, N, alpha, A, LDA);
}
