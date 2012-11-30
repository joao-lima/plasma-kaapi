/**
 *
 * @file core_dssssm.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @date 2010-11-15
 * @generated d Thu Sep 15 12:08:58 2011
 *
 **/
#include <cblas.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 *  CORE_dtstrf computes an LU factorization of a complex matrix formed
 *  by an upper triangular M1-by-N1 tile U on top of a M2-by-N2 tile A
 *  (N1 == N2) using partial pivoting with row interchanges.
 *
 *  This is the right-looking Level 2.5 BLAS version of the algorithm.
 *
 *******************************************************************************
 *
 * @param[in] M1
 *         The number of rows of the tile A1.  M1 >= 0.
 *
 * @param[in] N1
 *         The number of columns of the tile A1.  N1 >= 0.
 *
 * @param[in] M2
 *         The number of rows of the tile A2.  M2 >= 0.
 *
 * @param[in] N2
 *         The number of columns of the tile A2.  N2 >= 0.
 *
 * @param[in] K
 *         The number of columns of the tiles L1 and L2.  K >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A1
 *         On entry, the M1-by-N1 tile A1.
 *         On exit, A1 is overwritten by the application of L.
 *
 * @param[in] LDA1
 *         The leading dimension of the array A1.  LDA1 >= max(1,M1).
 *
 * @param[in,out] A2
 *         On entry, the M2-by-N2 tile A2.
 *         On exit, A2 is overwritten by the application of L.
 *
 * @param[in] LDA2
 *         The leading dimension of the array A2.  LDA2 >= max(1,M2).
 *
 * @param[in] L1
 *         The IB-by-K lower triangular tile as returned by CORE_dtstrf.
 *
 * @param[in] LDL1
 *         The leading dimension of the array L1.  LDL1 >= max(1,IB).
 *
 * @param[in] L2
 *         The M2-by-N2 tile as returned by CORE_dtstrf.
 *
 * @param[in] LDL2
 *         The leading dimension of the array L2.  LDL2 >= max(1,M2).
 *
 * @param[in] IPIV
 *         as returned by CORE_dtstrf.
 *
 *******************************************************************************
 *
 * @return
 *         \retval PLASMA_SUCCESS successful exit
 *         \retval <0 if INFO = -k, the k-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dssssm = PCORE_dssssm
#define CORE_dssssm PCORE_dssssm
#endif
int CORE_dssssm(int M1, int N1, int M2, int N2, int K, int IB,
                double *A1, int LDA1,
                double *A2, int LDA2,
                double *L1, int LDL1,
                double *L2, int LDL2,
                int *IPIV)
{
    static double zone  = 1.0;
    static double mzone =-1.0;

    int i, ii, sb;
    int im, ip;

    /* Check input arguments */
    if (M1 < 0) {
        coreblas_error(1, "Illegal value of M1");
        return -1;
    }
    if (N1 < 0) {
        coreblas_error(2, "Illegal value of N1");
        return -2;
    }
    if (M2 < 0) {
        coreblas_error(3, "Illegal value of M2");
        return -3;
    }
    if (N2 < 0) {
        coreblas_error(4, "Illegal value of N2");
        return -4;
    }
    if (K < 0) {
        coreblas_error(5, "Illegal value of K");
        return -5;
    }
    if (IB < 0) {
        coreblas_error(6, "Illegal value of IB");
        return -6;
    }
    if (LDA1 < max(1,M1)) {
        coreblas_error(8, "Illegal value of LDA1");
        return -8;
    }
    if (LDA2 < max(1,M2)) {
        coreblas_error(10, "Illegal value of LDA2");
        return -10;
    }
    if (LDL1 < max(1,IB)) {
        coreblas_error(12, "Illegal value of LDL1");
        return -12;
    }
    if (LDL2 < max(1,M2)) {
        coreblas_error(14, "Illegal value of LDL2");
        return -14;
    }

    /* Quick return */
    if ((M1 == 0) || (N1 == 0) || (M2 == 0) || (N2 == 0) || (K == 0) || (IB == 0))
        return PLASMA_SUCCESS;

    ip = 0;

    for(ii = 0; ii < K; ii += IB) {
        sb = min(K-ii, IB);

        for(i = 0; i < sb; i++) {
            im = IPIV[ip]-1;

            if (im != (ii+i)) {
                im = im - M1;
                cblas_dswap(N1, &A1[ii+i], LDA1, &A2[im], LDA2);
            }
            ip = ip + 1;
        }

        cblas_dtrsm(
            CblasColMajor, CblasLeft, CblasLower,
            CblasNoTrans, CblasUnit,
            sb, N1, (zone),
            &L1[LDL1*ii], LDL1,
            &A1[ii], LDA1);

        cblas_dgemm(
            CblasColMajor, CblasNoTrans, CblasNoTrans,
            M2, N2, sb,
            (mzone), &L2[LDL2*ii], LDL2,
            &A1[ii], LDA1,
            (zone), A2, LDA2);
    }
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dssssm(Quark *quark, Quark_Task_Flags *task_flags,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       double *L1, int ldl1,
                       double *L2, int ldl2,
                       int *IPIV)
{
    DAG_CORE_SSSSM;
    QUARK_Insert_Task(quark, CORE_dssssm_quark, task_flags,
        sizeof(int),                        &m1,    VALUE,
        sizeof(int),                        &n1,    VALUE,
        sizeof(int),                        &m2,    VALUE,
        sizeof(int),                        &n2,    VALUE,
        sizeof(int),                        &k,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(double)*nb*nb,    A1,            INOUT,
        sizeof(int),                        &lda1,  VALUE,
        sizeof(double)*nb*nb,    A2,            INOUT | LOCALITY,
        sizeof(int),                        &lda2,  VALUE,
        sizeof(double)*ib*nb,    L1,            INPUT,
        sizeof(int),                        &ldl1,  VALUE,
        sizeof(double)*ib*nb,    L2,            INPUT,
        sizeof(int),                        &ldl2,  VALUE,
        sizeof(int)*nb,                      IPIV,          INPUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dssssm_quark = PCORE_dssssm_quark
#define CORE_dssssm_quark PCORE_dssssm_quark
#endif
void CORE_dssssm_quark(Quark *quark)
{
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
    double *L1;
    int ldl1;
    double *L2;
    int ldl2;
    int *IPIV;

    quark_unpack_args_15(quark, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, L1, ldl1, L2, ldl2, IPIV);
    CORE_dssssm(m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, L1, ldl1, L2, ldl2, IPIV);
}
