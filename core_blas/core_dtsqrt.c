/**
 *
 * @file core_dtsqrt.c
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
 * @generated d Thu Sep 15 12:08:56 2011
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
 * CORE_dtsqrt computes a QR factorization of a rectangular matrix
 * formed by coupling a complex N-by-N upper triangular tile A1
 * on top of a complex M-by-N tile A2:
 *
 *    | A1 | = Q * R
 *    | A2 |
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of columns of the tile A2. M >= 0.
 *
 * @param[in] N
 *         The number of rows of the tile A1.
 *         The number of columns of the tiles A1 and A2. N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A1
 *         On entry, the N-by-N tile A1.
 *         On exit, the elements on and above the diagonal of the array
 *         contain the N-by-N upper trapezoidal tile R;
 *         the elements below the diagonal are not referenced.
 *
 * @param[in] LDA1
 *         The leading dimension of the array A1. LDA1 >= max(1,N).
 *
 * @param[in,out] A2
 *         On entry, the M-by-N tile A2.
 *         On exit, all the elements with the array TAU, represent
 *         the unitary tile Q as a product of elementary reflectors
 *         (see Further Details).
 *
 * @param[in] LDA2
 *         The leading dimension of the tile A2. LDA2 >= max(1,M).
 *
 * @param[out] T
 *         The IB-by-N triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[out] TAU
 *         The scalar factors of the elementary reflectors (see Further
 *         Details).
 *
 * @param[out] WORK
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dtsqrt = PCORE_dtsqrt
#define CORE_dtsqrt PCORE_dtsqrt
#define CORE_dtsmqr PCORE_dtsmqr
int  CORE_dtsmqr(int side, int trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 double *A1, int LDA1,
                 double *A2, int LDA2,
                 double *V, int LDV,
                 double *T, int LDT,
                 double *WORK, int LDWORK);
#endif

int CORE_dtsqrt(int M, int N, int IB,
                double *A1, int LDA1,
                double *A2, int LDA2,
                double *T, int LDT,
                double *TAU, double *WORK)
{
    static double zone  = 1.0;
    static double zzero = 0.0;

    double alpha;
    int i, ii, sb;

    /* Check input arguments */
    if (M < 0) {
        coreblas_error(1, "Illegal value of M");
        return -1;
    }
    if (N < 0) {
        coreblas_error(2, "Illegal value of N");
        return -2;
    }
    if (IB < 0) {
        coreblas_error(3, "Illegal value of IB");
        return -3;
    }
    if ((LDA2 < max(1,M)) && (M > 0)) {
        coreblas_error(8, "Illegal value of LDA2");
        return -8;
    }

    /* Quick return */
    if ((M == 0) || (N == 0) || (IB == 0))
        return PLASMA_SUCCESS;

    for(ii = 0; ii < N; ii += IB) {
        sb = min(N-ii, IB);
        for(i = 0; i < sb; i++) {
            /*
             * Generate elementary reflector H( II*IB+I ) to annihilate
             * A( II*IB+I:M, II*IB+I )
             */
            LAPACKE_dlarfg_work(M+1, &A1[LDA1*(ii+i)+ii+i], &A2[LDA2*(ii+i)], 1, &TAU[ii+i]);

            if (ii+i+1 < N) {
                /*
                 * Apply H( II*IB+I ) to A( II*IB+I:M, II*IB+I+1:II*IB+IB ) from the left
                 */
                alpha = -(TAU[ii+i]);
                cblas_dcopy(
                    sb-i-1,
                    &A1[LDA1*(ii+i+1)+(ii+i)], LDA1,
                    WORK, 1);
#ifdef COMPLEX
                LAPACKE_dlacgv_work(sb-i-1, WORK, 1);
#endif
                cblas_dgemv(
                    CblasColMajor, (CBLAS_TRANSPOSE)PlasmaTrans,
                    M, sb-i-1,
                    (zone), &A2[LDA2*(ii+i+1)], LDA2,
                    &A2[LDA2*(ii+i)], 1,
                    (zone), WORK, 1);
#ifdef COMPLEX
                LAPACKE_dlacgv_work(sb-i-1, WORK, 1 );
#endif
                cblas_daxpy(
                    sb-i-1, (alpha),
                    WORK, 1,
                    &A1[LDA1*(ii+i+1)+ii+i], LDA1);
#ifdef COMPLEX
                LAPACKE_dlacgv_work(sb-i-1, WORK, 1 );
#endif
                cblas_dger(
                    CblasColMajor, M, sb-i-1, (alpha),
                    &A2[LDA2*(ii+i)], 1,
                    WORK, 1,
                    &A2[LDA2*(ii+i+1)], LDA2);
            }
            /*
             * Calculate T
             */
            alpha = -TAU[ii+i];
            cblas_dgemv(
                CblasColMajor, (CBLAS_TRANSPOSE)PlasmaTrans, M, i,
                (alpha), &A2[LDA2*ii], LDA2,
                &A2[LDA2*(ii+i)], 1,
                (zzero), &T[LDT*(ii+i)], 1);

            cblas_dtrmv(
                CblasColMajor, (CBLAS_UPLO)PlasmaUpper,
                (CBLAS_TRANSPOSE)PlasmaNoTrans, (CBLAS_DIAG)PlasmaNonUnit, i,
                &T[LDT*ii], LDT,
                &T[LDT*(ii+i)], 1);

            T[LDT*(ii+i)+i] = TAU[ii+i];
        }
        if (N > ii+sb) {
            CORE_dtsmqr(
                PlasmaLeft, PlasmaTrans,
                sb, N-(ii+sb), M, N-(ii+sb), IB, IB,
                &A1[LDA1*(ii+sb)+ii], LDA1,
                &A2[LDA2*(ii+sb)], LDA2,
                &A2[LDA2*ii], LDA2,
                &T[LDT*ii], LDT,
                WORK, sb);
        }
    }
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dtsqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       double *T, int ldt)
{
    DAG_CORE_TSQRT;
    QUARK_Insert_Task(quark, CORE_dtsqrt_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(double)*nb*nb,    A1,            INOUT | QUARK_REGION_D | QUARK_REGION_U,
        sizeof(int),                        &lda1,  VALUE,
        sizeof(double)*nb*nb,    A2,            INOUT | LOCALITY,
        sizeof(int),                        &lda2,  VALUE,
        sizeof(double)*ib*nb,    T,             OUTPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(double)*nb,       NULL,          SCRATCH,
        sizeof(double)*ib*nb,    NULL,          SCRATCH,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dtsqrt_quark = PCORE_dtsqrt_quark
#define CORE_dtsqrt_quark PCORE_dtsqrt_quark
#endif
void CORE_dtsqrt_quark(Quark *quark)
{
    int m;
    int n;
    int ib;
    double *A1;
    int lda1;
    double *A2;
    int lda2;
    double *T;
    int ldt;
    double *TAU;
    double *WORK;

    quark_unpack_args_11(quark, m, n, ib, A1, lda1, A2, lda2, T, ldt, TAU, WORK);
    CORE_dtsqrt(m, n, ib, A1, lda1, A2, lda2, T, ldt, TAU, WORK);
}
