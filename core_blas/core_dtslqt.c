/**
 *
 * @file core_dtslqt.c
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
 * @generated d Thu Sep 15 12:08:57 2011
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
 *  CORE_dtslqt computes a LQ factorization of a rectangular matrix
 *  formed by coupling side-by-side a complex M-by-M
 *  lower triangular tile A1 and a complex M-by-N tile A2:
 *
 *    | A1 A2 | = L * Q
 *
 *  The tile Q is represented as a product of elementary reflectors
 *
 *    Q = H(k)' . . . H(2)' H(1)', where k = min(M,N).
 *
 *  Each H(i) has the form
 *
 *    H(i) = I - tau * v * v'
 *
 *  where tau is a complex scalar, and v is a complex vector with
 *  v(1:i-1) = 0 and v(i) = 1; g(v(i+1:n)) is stored on exit in
 *  A2(i,1:n), and tau in TAU(i).
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of rows of the tile A1 and A2. M >= 0.
 *         The number of columns of the tile A1.
 *
 * @param[in] N
 *         The number of columns of the tile A2. N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A1
 *         On entry, the M-by-M tile A1.
 *         On exit, the elements on and below the diagonal of the array
 *         contain the M-by-M lower trapezoidal tile L;
 *         the elements above the diagonal are not referenced.
 *
 * @param[in] LDA1
 *         The leading dimension of the array A1. LDA1 >= max(1,M).
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
#pragma weak CORE_dtslqt = PCORE_dtslqt
#define CORE_dtslqt PCORE_dtslqt
#define CORE_dtsmlq PCORE_dtsmlq
int  CORE_dtsmlq(int side, int trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 double *A1, int LDA1,
                 double *A2, int LDA2,
                 double *V, int LDV,
                 double *T, int LDT,
                 double *WORK, int LDWORK);
#endif
int CORE_dtslqt(int M, int N, int IB,
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

    for(ii = 0; ii < M; ii += IB) {
        sb = min(M-ii, IB);
        for(i = 0; i < sb; i++) {
            /*
             * Generate elementary reflector H( II*IB+I ) to annihilate A( II*IB+I, II*IB+I:N ).
             */
#ifdef COMPLEX
            LAPACKE_dlacgv_work(N, &A2[ii+i], LDA2);
            LAPACKE_dlacgv_work(1, &A1[LDA1*(ii+i)+ii+i], LDA1);
#endif
            LAPACKE_dlarfg_work(N+1, &A1[LDA1*(ii+i)+ii+i], &A2[ii+i], LDA2, &TAU[ii+i]);

            alpha = -(TAU[ii+i]);
            if (ii+i+1 < M) {
                /*
                 * Apply H( II+I-1 ) to A( II+I:II+IB-1, II+I-1:N  ) from the right.
                 */
                cblas_dcopy(
                    sb-i-1,
                    &A1[LDA1*(ii+i)+(ii+i+1)], 1,
                    WORK, 1);

                cblas_dgemv(
                    CblasColMajor, (CBLAS_TRANSPOSE)PlasmaNoTrans,
                    sb-i-1, N,
                    (zone), &A2[ii+i+1], LDA2,
                    &A2[ii+i], LDA2,
                    (zone), WORK, 1);

                cblas_daxpy(
                    sb-i-1, (alpha),
                    WORK, 1,
                    &A1[LDA1*(ii+i)+ii+i+1], 1);

                cblas_dger(
                    CblasColMajor, sb-i-1, N,
                    (alpha), WORK, 1,
                    &A2[ii+i], LDA2,
                    &A2[ii+i+1], LDA2);
            }
            /*
             * Calculate T.
             */
            cblas_dgemv(
                CblasColMajor, (CBLAS_TRANSPOSE)PlasmaNoTrans, i, N,
                (alpha), &A2[ii], LDA2,
                &A2[ii+i], LDA2,
                (zzero), &T[LDT*(ii+i)], 1);
#ifdef COMPLEX
            LAPACKE_dlacgv_work(N, &A2[ii+i], LDA2 );
            LAPACKE_dlacgv_work(1, &A1[LDA1*(ii+i)+ii+i], LDA1 );
#endif
            cblas_dtrmv(
                CblasColMajor, (CBLAS_UPLO)PlasmaUpper,
                (CBLAS_TRANSPOSE)PlasmaNoTrans, (CBLAS_DIAG)PlasmaNonUnit, i,
                &T[LDT*ii], LDT,
                &T[LDT*(ii+i)], 1);

            T[LDT*(ii+i)+i] = TAU[ii+i];
        }
        if (M > ii+sb) {
            CORE_dtsmlq(
                PlasmaRight, PlasmaTrans,
                M-(ii+sb), sb, M-(ii+sb), N, IB, IB,
                &A1[LDA1*ii+ii+sb], LDA1,
                &A2[ii+sb], LDA2,
                &A2[ii], LDA2,
                &T[LDT*ii], LDT,
                WORK, LDA1);
        }
    }
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dtslqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       double *T, int ldt)
{
    DAG_CORE_TSLQT;
    QUARK_Insert_Task(quark, CORE_dtslqt_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(double)*nb*nb,    A1,            INOUT | QUARK_REGION_D | QUARK_REGION_L,
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
#pragma weak CORE_dtslqt_quark = PCORE_dtslqt_quark
#define CORE_dtslqt_quark PCORE_dtslqt_quark
#endif
void CORE_dtslqt_quark(Quark *quark)
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
    CORE_dtslqt(m, n, ib, A1, lda1, A2, lda2, T, ldt, TAU, WORK);
}
