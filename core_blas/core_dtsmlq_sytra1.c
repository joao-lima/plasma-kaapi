/**
 *
 * @file core_dtsmlq_sytra1.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Azzam Haidar
 * @date 2010-11-15
 * @generated d Thu Sep 15 12:09:03 2011
 *
 **/
#include <lapacke.h>
#include "common.h"
#undef COMPLEX
#define REAL

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 * CORE_dtsmlq_sytra1: see CORE_dtsmlq
 *
 * This kernel applies a Right transformation on | A1' A2 |
 * and does not handle the transpose of A1.
 * Needs therefore to make the explicit transpose of A1 before
 * and after the application of the block of reflectors
 * Can be further optimized by changing accordingly the underneath
 * kernel ztsrfb!
 *
 *******************************************************************************
 *
 * @param[in] side
 *         @arg PlasmaLeft  : apply Q or Q**T from the Left;
 *         @arg PlasmaRight : apply Q or Q**T from the Right.
 *
 * @param[in] trans
 *         @arg PlasmaNoTrans   :  No transpose, apply Q;
 *         @arg PlasmaTrans :  ConjTranspose, apply Q**T.
 *
 * @param[in] M1
 *         The number of rows of the tile A1. M1 >= 0.
 *
 * @param[in] N1
 *         The number of columns of the tile A1. N1 >= 0.
 *
 * @param[in] M2
 *         The number of rows of the tile A2. M2 >= 0.
 *         M2 = M1 if side == PlasmaRight.
 *
 * @param[in] N2
 *         The number of columns of the tile A2. N2 >= 0.
 *         N2 = N1 if side == PlasmaLeft.
 *
 * @param[in] K
 *         The number of elementary reflectors whose product defines
 *         the matrix Q.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A1
 *         On entry, the M1-by-N1 tile A1.
 *         On exit, A1 is overwritten by the application of Q.
 *
 * @param[in] LDA1
 *         The leading dimension of the array A1. LDA1 >= max(1,M1).
 *
 * @param[in,out] A2
 *         On entry, the M2-by-N2 tile A2.
 *         On exit, A2 is overwritten by the application of Q.
 *
 * @param[in] LDA2
 *         The leading dimension of the tile A2. LDA2 >= max(1,M2).
 *
 * @param[in] V
 *         The i-th row must contain the vector which defines the
 *         elementary reflector H(i), for i = 1,2,...,k, as returned by
 *         CORE_DTSLQT in the first k rows of its array argument V.
 *
 * @param[in] LDV
 *         The leading dimension of the array V. LDV >= max(1,K).
 *
 * @param[out] T
 *         The IB-by-N1 triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[out] WORK
 *         Workspace array of size 
 *             LDWORK-by-M1 if side == PlasmaLeft
 *             LDWORK-by-IB if side == PlasmaRight
 *
 * @param[in] LDWORK
 *         The leading dimension of the array WORK. 
 *             LDWORK >= max(1,IB) if side == PlasmaLeft
 *             LDWORK >= max(1,N1) if side == PlasmaRight
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dtsmlq_sytra1 = PCORE_dtsmlq_sytra1
#define CORE_dtsmlq_sytra1 PCORE_dtsmlq_sytra1
#define CORE_dtsmlq PCORE_dtsmlq
int  CORE_dtsmlq(int side, int trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 double *A1, int LDA1,
                 double *A2, int LDA2,
                 double *V, int LDV,
                 double *T, int LDT,
                 double *WORK, int LDWORK);
#endif
int CORE_dtsmlq_sytra1( int side, int trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        double *A1, int lda1,
                        double *A2, int lda2,
                        double *V, int ldv,
                        double *T, int ldt,
                        double *WORK, int ldwork)
{
    int i, j;

    if ( (m1 != n1) ) {
        coreblas_error(3, "Illegal value of M1, N1");
        return -3;
    }

    /* in-place transposition of A1 */
    for (j = 0; j < n1; j++){
        A1[j + j*lda1] = (A1[j + j*lda1]);

        for (i = j+1; i < m1; i++){
            *WORK = *(A1 + i + j*lda1);
            *(A1 + i + j*lda1) = (*(A1 + j + i*lda1));
            *(A1 + j + i*lda1) = (*WORK);
        }
    }

    CORE_dtsmlq(side, trans, m1, n1, m2, n2, k, ib, 
                A1, lda1, A2, lda2, 
                V,  ldv,  T,  ldt, 
                WORK, ldwork);

    /* in-place transposition of A1 */
    for (j = 0; j < n1; j++){
        A1[j + j*lda1] = (A1[j + j*lda1]);

        for (i = j+1; i < m1; i++){
            *WORK = *(A1 + i + j*lda1);
            *(A1 + i + j*lda1) = (*(A1 + j + i*lda1));
            *(A1 + j + i*lda1) = (*WORK);
        }
    }
    
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dtsmlq_sytra1(Quark *quark, Quark_Task_Flags *task_flags,
                              int side, int trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              double *A1, int lda1,
                              double *A2, int lda2,
                              double *V, int ldv,
                              double *T, int ldt)
{
    int ldwork = side == PlasmaLeft ? ib : nb;

    QUARK_Insert_Task(quark, CORE_dtsmlq_sytra1_quark, task_flags,
        sizeof(PLASMA_enum),                &side,  VALUE,
        sizeof(PLASMA_enum),                &trans, VALUE,
        sizeof(int),                        &m1,    VALUE,
        sizeof(int),                        &n1,    VALUE,
        sizeof(int),                        &m2,    VALUE,
        sizeof(int),                        &n2,    VALUE,
        sizeof(int),                        &k,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(double)*nb*nb,    A1,            INOUT|QUARK_REGION_U|QUARK_REGION_D,
        sizeof(int),                        &lda1,  VALUE,
        sizeof(double)*nb*nb,    A2,            INOUT,
        sizeof(int),                        &lda2,  VALUE,
        sizeof(double)*nb*nb,    V,             INPUT,
        sizeof(int),                        &ldv,   VALUE,
        sizeof(double)*ib*nb,    T,             INPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(double)*ib*nb,    NULL,          SCRATCH,
        sizeof(int),                        &ldwork, VALUE,
        0);
}

/***************************************************************************//**
 * This kernel applies a Right transformation on | A1' A2 |
 * and does not handle the transpose of A1.
 * Needs therefore to make the explicit transpose of A1 before
 * and after the application of the block of reflectors
 * Can be further optimized by changing accordingly the underneath
 * kernel ztsrfb!
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dtsmlq_sytra1_quark = PCORE_dtsmlq_sytra1_quark
#define CORE_dtsmlq_sytra1_quark PCORE_dtsmlq_sytra1_quark
#endif
void CORE_dtsmlq_sytra1_quark(Quark *quark)
{
    int side;
    int trans;
    int m1;
    int n1;
    int m2;
    int n2;
    int k;
    int ib;
    double *A1;
    int lda1;
    double *A2;
    int lda2;
    double *V;
    int ldv;
    double *T;
    int ldt;
    double *WORK;
    int ldwork;

    quark_unpack_args_18(quark, side, trans, m1, n1, m2, n2, k, ib, 
                         A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
    CORE_dtsmlq_sytra1(side, trans, m1, n1, m2, n2, k, ib, 
                       A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
}

