/**
 *
 * @file core_sttlqt.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Thu Sep 15 12:08:57 2011
 *
 **/
#include <lapacke.h>
#include "common.h"
#undef COMPLEX
#define REAL

/***************************************************************************//**
 *
 * @ingroup CORE_float
 *
 *  CORE_sttlqt computes a LQ factorization of a rectangular matrix
 *  formed by coupling side-by-side a complex M-by-M lower triangular tile A1
 *  and a complex M-by-N lower triangular tile A2:
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
 *         The leading dimension of the array A1.  LDA1 >= max(1,N).
 *
 * @param[in,out] A2
 *         On entry, the M-by-N lower triangular tile A2.
 *         On exit, the elements on and below the diagonal of the array
 *         with the array TAU, represent
 *         the unitary tile Q as a product of elementary reflectors
 *         (see Further Details).
 *
 * @param[in] LDA2
 *         The leading dimension of the array A2.  LDA2 >= max(1,M).
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
 * @param[in,out] WORK
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sttlqt = PCORE_sttlqt
#define CORE_sttlqt PCORE_sttlqt
#define CORE_sttrfb PCORE_sttrfb
int  CORE_sttrfb(int side, int trans, int direct, int storev,
                 int M1, int N1, int M2, int N2, int IB,
                 float *A1, int LDA1,
                 float *A2, int LDA2,
                 float *V, int LDV,
                 float *T, int LDT,
                 float *WORK, int LDWORK);
#endif
int CORE_sttlqt(int M, int N, int IB,
                float *A1, int LDA1,
                float *A2, int LDA2,
                float *T, int LDT,
                float *TAU, float *WORK)
{
    static float zone  = 1.0;
    static float zzero = 0.0;

    float alpha;
    int i, j, ii, sb, mi, ni;

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
        coreblas_error(7, "Illegal value of LDA2");
        return -7;
    }

    /* Quick return */
    if ((M == 0) || (N == 0) || (IB == 0))
        return PLASMA_SUCCESS;

    for(ii = 0; ii < M; ii += IB) {
        sb = min(M-ii, IB);
        for(i = 0; i < sb; i++) {
            /*
             * Generate elementary reflector H( II*IB+I ) to annihilate A( II*IB+I, II*IB+I:M ).
             */
            ni = ii+i+1;
#ifdef COMPLEX
            LAPACKE_slacgv_work(ni, &A2[ii+i], LDA2);
            LAPACKE_slacgv_work(1,  &A1[LDA1*(ii+i)+ii+i], LDA1);
#endif
            //tmp = ii+i+2;
            LAPACKE_slarfg_work(ni+1, &A1[LDA1*(ii+i)+ii+i], &A2[ii+i], LDA2, &TAU[ii+i]);

            if (sb-i-1 > 0) { //*DB
            //if (ii+i+1 < M) {
                /*
                 * Apply H( II+I-1 ) to A( II+I:II+IB-1, II+I-1:M  ) from the right.
                 */
                mi = sb-i-1;
                cblas_scopy(
                    mi,
                    &A1[LDA1*(ii+i)+(ii+i+1)], 1,
                    WORK, 1);

                cblas_sgemv(
                    CblasColMajor, (CBLAS_TRANSPOSE)PlasmaNoTrans,
                    mi, ni,
                    (zone), &A2[ii+i+1], LDA2,
                    &A2[ii+i], LDA2,
                    (zone), WORK, 1);

                alpha = -(TAU[ii+i]);
                cblas_saxpy(
                    mi, (alpha),
                    WORK, 1,
                    &A1[LDA1*(ii+i)+ii+i+1], 1);

                cblas_sger(
                    CblasColMajor, mi, ni,
                    (alpha), WORK, 1,
                    &A2[ii+i], LDA2,
                    &A2[ii+i+1], LDA2);
            }
            /*
             * Calculate T.
             */

            if (i > 0 ) {

                cblas_scopy(i, &A2[LDA2*ii+ii+i], LDA2, &WORK[ii], 1);

                cblas_strmv(
                    CblasColMajor, (CBLAS_UPLO)PlasmaLower,
                    (CBLAS_TRANSPOSE)PlasmaNoTrans, (CBLAS_DIAG)PlasmaNonUnit,
                    i, &A2[LDA2*ii+ii], LDA2,
                    &WORK[ii], 1);

                alpha = -(TAU[ii+i]);

                for(j = 0; j < i; j++) {
                    WORK[ii+j] = alpha * WORK[ii+j];
                }

                if (ii > 0) {
                    cblas_sgemv(
                        CblasColMajor, (CBLAS_TRANSPOSE)PlasmaNoTrans, i, ii,
                        (alpha), &A2[ii], LDA2,
                        &A2[ii+i], LDA2,
                        (zzero), WORK, 1);

                    cblas_saxpy(i, (zone), &WORK[ii], 1, WORK, 1);
                }

                cblas_scopy(i, WORK, 1, &T[LDT*(ii+i)], 1);

                cblas_strmv(
                    CblasColMajor, (CBLAS_UPLO)PlasmaUpper,
                    (CBLAS_TRANSPOSE)PlasmaNoTrans, (CBLAS_DIAG)PlasmaNonUnit,
                    i, &T[LDT*ii], LDT,
                    &T[LDT*(ii+i)], 1);

            }

#ifdef COMPLEX
            LAPACKE_slacgv_work(ni, &A2[ii+i], LDA2 );
            LAPACKE_slacgv_work(1,  &A1[LDA1*(ii+i)+ii+i], LDA1 );
#endif

            T[LDT*(ii+i)+i] = TAU[ii+i];
        }

        /* Apply Q' to the rest of the matrix to the right */
        if (M > ii+sb) {
            CORE_sttrfb(
                PlasmaRight, PlasmaNoTrans,
                PlasmaForward, PlasmaRowwise,
                M-(ii+sb), sb, M-(ii+sb), ii+sb, sb,
                &A1[LDA1*ii+ii+sb], LDA1,
                &A2[ii+sb], LDA2,
                &A2[ii], LDA2,
                &T[LDT*ii], LDT,
                WORK, N);

        }
    }

    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sttlqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       float *A1, int lda1,
                       float *A2, int lda2,
                       float *T, int ldt)
{
    DAG_CORE_TTLQT;
    QUARK_Insert_Task(quark, CORE_sttlqt_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(float)*nb*nb,    A1,            INOUT|QUARK_REGION_D|QUARK_REGION_L,
        sizeof(int),                        &lda1,  VALUE,
        sizeof(float)*nb*nb,    A2,            INOUT|QUARK_REGION_D|QUARK_REGION_L,
        sizeof(int),                        &lda2,  VALUE,
        sizeof(float)*ib*nb,    T,             OUTPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(float)*nb,       NULL,          SCRATCH,
        sizeof(float)*ib*nb,    NULL,          SCRATCH,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sttlqt_quark = PCORE_sttlqt_quark
#define CORE_sttlqt_quark PCORE_sttlqt_quark
#endif
void CORE_sttlqt_quark(Quark *quark)
{
    int m;
    int n;
    int ib;
    float *A1;
    int lda1;
    float *A2;
    int lda2;
    float *T;
    int ldt;
    float *TAU;
    float *WORK;

    quark_unpack_args_11(quark, m, n, ib, A1, lda1, A2, lda2, T, ldt, TAU, WORK);
    CORE_sttlqt(m, n, ib, A1, lda1, A2, lda2, T, ldt, TAU, WORK);
}
