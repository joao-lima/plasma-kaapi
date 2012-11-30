/**
 *
 * @file core_claset.c
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
 *  CORE_claset - Sets the elements of the matrix A on the diagonal
 *  to beta and on the off-diagonals to alpha
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies which elements of the matrix are to be set 
 *          = PlasmaUpper: Upper part of A is set;
 *          = PlasmaLower: Lower part of A is set;
 *          = PlasmaUpperLower: ALL elements of A are set.
 *
 * @param[in] M
 *          The number of rows of the matrix A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the matrix A.  N >= 0.
 *
 * @param[in] alpha
 *         The constant to which the off-diagonal elements are to be set.
 *
 * @param[in] beta
 *         The constant to which the diagonal elements are to be set.
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile A.
 *         On exit, A has been set accordingly.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_claset = PCORE_claset
#define CORE_claset PCORE_claset
#endif
void CORE_claset(PLASMA_enum uplo, int M, int N,
                 PLASMA_Complex32_t alpha, PLASMA_Complex32_t beta, 
                 PLASMA_Complex32_t *A, int LDA)
{
    LAPACKE_claset_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        M, N, alpha, beta, A, LDA);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_claset(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int M, int N,
                       PLASMA_Complex32_t alpha, PLASMA_Complex32_t beta, 
                       PLASMA_Complex32_t *A, int LDA)
{
    QUARK_Insert_Task(quark, CORE_claset_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(int),                        &M,     VALUE,
        sizeof(int),                        &N,     VALUE,
        sizeof(PLASMA_Complex32_t),         &alpha, VALUE,
        sizeof(PLASMA_Complex32_t),         &beta,  VALUE,
        sizeof(PLASMA_Complex32_t)*M*N,     A,      OUTPUT,
        sizeof(int),                        &LDA,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_claset_quark = PCORE_claset_quark
#define CORE_claset_quark PCORE_claset_quark
#endif
void CORE_claset_quark(Quark *quark)
{
    int uplo;
    int M;
    int N;
    PLASMA_Complex32_t alpha;
    PLASMA_Complex32_t beta;
    PLASMA_Complex32_t *A;
    int LDA;

    quark_unpack_args_7(quark, uplo, M, N, alpha, beta, A, LDA);
    LAPACKE_claset_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        M, N, alpha, beta, A, LDA);
}
