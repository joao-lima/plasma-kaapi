/**
 *
 * @file core_zttlqt.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include <lapacke.h>
#include "common.h"
#undef REAL
#define COMPLEX

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 *  CORE_zttlqt computes a LQ factorization of a rectangular matrix
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
 *  v(1:i-1) = 0 and v(i) = 1; conjg(v(i+1:n)) is stored on exit in
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
#pragma weak CORE_zttlqt = PCORE_zttlqt
#define CORE_zttlqt PCORE_zttlqt
#define CORE_zttrfb PCORE_zttrfb
int  CORE_zttrfb(int side, int trans, int direct, int storev,
                 int M1, int N1, int M2, int N2, int IB,
                 PLASMA_Complex64_t *A1, int LDA1,
                 PLASMA_Complex64_t *A2, int LDA2,
                 PLASMA_Complex64_t *V, int LDV,
                 PLASMA_Complex64_t *T, int LDT,
                 PLASMA_Complex64_t *WORK, int LDWORK);
#endif
int CORE_zttlqt(int M, int N, int IB,
                PLASMA_Complex64_t *A1, int LDA1,
                PLASMA_Complex64_t *A2, int LDA2,
                PLASMA_Complex64_t *T, int LDT,
                PLASMA_Complex64_t *TAU, PLASMA_Complex64_t *WORK)
{
    static PLASMA_Complex64_t zone  = 1.0;
    static PLASMA_Complex64_t zzero = 0.0;

    PLASMA_Complex64_t alpha;
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
            LAPACKE_zlacgv_work(ni, &A2[ii+i], LDA2);
            LAPACKE_zlacgv_work(1,  &A1[LDA1*(ii+i)+ii+i], LDA1);
#endif
            //tmp = ii+i+2;
            LAPACKE_zlarfg_work(ni+1, &A1[LDA1*(ii+i)+ii+i], &A2[ii+i], LDA2, &TAU[ii+i]);

            if (sb-i-1 > 0) { //*DB
            //if (ii+i+1 < M) {
                /*
                 * Apply H( II+I-1 ) to A( II+I:II+IB-1, II+I-1:M  ) from the right.
                 */
                mi = sb-i-1;
                cblas_zcopy(
                    mi,
                    &A1[LDA1*(ii+i)+(ii+i+1)], 1,
                    WORK, 1);

                cblas_zgemv(
                    CblasColMajor, (CBLAS_TRANSPOSE)PlasmaNoTrans,
                    mi, ni,
                    CBLAS_SADDR(zone), &A2[ii+i+1], LDA2,
                    &A2[ii+i], LDA2,
                    CBLAS_SADDR(zone), WORK, 1);

                alpha = -(TAU[ii+i]);
                cblas_zaxpy(
                    mi, CBLAS_SADDR(alpha),
                    WORK, 1,
                    &A1[LDA1*(ii+i)+ii+i+1], 1);

                cblas_zgerc(
                    CblasColMajor, mi, ni,
                    CBLAS_SADDR(alpha), WORK, 1,
                    &A2[ii+i], LDA2,
                    &A2[ii+i+1], LDA2);
            }
            /*
             * Calculate T.
             */

            if (i > 0 ) {

                cblas_zcopy(i, &A2[LDA2*ii+ii+i], LDA2, &WORK[ii], 1);

                cblas_ztrmv(
                    CblasColMajor, (CBLAS_UPLO)PlasmaLower,
                    (CBLAS_TRANSPOSE)PlasmaNoTrans, (CBLAS_DIAG)PlasmaNonUnit,
                    i, &A2[LDA2*ii+ii], LDA2,
                    &WORK[ii], 1);

                alpha = -(TAU[ii+i]);

                for(j = 0; j < i; j++) {
                    WORK[ii+j] = alpha * WORK[ii+j];
                }

                if (ii > 0) {
                    cblas_zgemv(
                        CblasColMajor, (CBLAS_TRANSPOSE)PlasmaNoTrans, i, ii,
                        CBLAS_SADDR(alpha), &A2[ii], LDA2,
                        &A2[ii+i], LDA2,
                        CBLAS_SADDR(zzero), WORK, 1);

                    cblas_zaxpy(i, CBLAS_SADDR(zone), &WORK[ii], 1, WORK, 1);
                }

                cblas_zcopy(i, WORK, 1, &T[LDT*(ii+i)], 1);

                cblas_ztrmv(
                    CblasColMajor, (CBLAS_UPLO)PlasmaUpper,
                    (CBLAS_TRANSPOSE)PlasmaNoTrans, (CBLAS_DIAG)PlasmaNonUnit,
                    i, &T[LDT*ii], LDT,
                    &T[LDT*(ii+i)], 1);

            }

#ifdef COMPLEX
            LAPACKE_zlacgv_work(ni, &A2[ii+i], LDA2 );
            LAPACKE_zlacgv_work(1,  &A1[LDA1*(ii+i)+ii+i], LDA1 );
#endif

            T[LDT*(ii+i)+i] = TAU[ii+i];
        }

        /* Apply Q' to the rest of the matrix to the right */
        if (M > ii+sb) {
            CORE_zttrfb(
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
void QUARK_CORE_zttlqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       PLASMA_Complex64_t *A1, int lda1,
                       PLASMA_Complex64_t *A2, int lda2,
                       PLASMA_Complex64_t *T, int ldt)
{
    DAG_CORE_TTLQT;
    QUARK_Insert_Task(quark, CORE_zttlqt_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A1,            INOUT|QUARK_REGION_D|QUARK_REGION_L,
        sizeof(int),                        &lda1,  VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A2,            INOUT|QUARK_REGION_D|QUARK_REGION_L,
        sizeof(int),                        &lda2,  VALUE,
        sizeof(PLASMA_Complex64_t)*ib*nb,    T,             OUTPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(PLASMA_Complex64_t)*nb,       NULL,          SCRATCH,
        sizeof(PLASMA_Complex64_t)*ib*nb,    NULL,          SCRATCH,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zttlqt_quark = PCORE_zttlqt_quark
#define CORE_zttlqt_quark PCORE_zttlqt_quark
#endif
void CORE_zttlqt_quark(Quark *quark)
{
    int m;
    int n;
    int ib;
    PLASMA_Complex64_t *A1;
    int lda1;
    PLASMA_Complex64_t *A2;
    int lda2;
    PLASMA_Complex64_t *T;
    int ldt;
    PLASMA_Complex64_t *TAU;
    PLASMA_Complex64_t *WORK;

    quark_unpack_args_11(quark, m, n, ib, A1, lda1, A2, lda2, T, ldt, TAU, WORK);
    CORE_zttlqt(m, n, ib, A1, lda1, A2, lda2, T, ldt, TAU, WORK);
}
