/**
 *
 * @file core_cttqrt.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Dulceneia Becker
 * @date 2010-11-15
 * @generated c Thu Sep 15 12:08:57 2011
 *
 **/
#include <lapacke.h>
#include "common.h"
#undef REAL
#define COMPLEX

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex32_t
 *
 *  CORE_cttqrt computes a QR factorization of a rectangular matrix
 *  formed by coupling a complex N-by-N upper triangular tile A1
 *  on top of a complex M-by-N upper triangular tile A2:
 *
 *    | A1 | = Q * R
 *    | A2 |
 *
 *  The tile Q is represented as a product of elementary reflectors
 *
 *    Q = H(1) H(2) . . . H(k), where k = min(M,N).
 *
 *  Each H(i) has the form
 *
 *    H(i) = I - tau * v * v'
 *
 *  where tau is a complex scalar, and v is a complex vector with
 *  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A2(1:m,i),
 *  and tau in TAU(i).
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of rows of the tile A2.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A1 and A2.  N >= 0.
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
 *         The leading dimension of the array A1.  LDA1 >= max(1,N).
 *
 * @param[in,out] A2
 *         On entry, the M-by-N upper triangular tile A2.
 *         On exit, the elements on and above the diagonal of the array
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
#pragma weak CORE_cttqrt = PCORE_cttqrt
#define CORE_cttqrt PCORE_cttqrt
#define CORE_cttrfb PCORE_cttrfb
int  CORE_cttrfb(int side, int trans, int direct, int storev,
                 int M1, int N1, int M2, int N2, int IB,
                 PLASMA_Complex32_t *A1, int LDA1,
                 PLASMA_Complex32_t *A2, int LDA2,
                 PLASMA_Complex32_t *V, int LDV,
                 PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *WORK, int LDWORK);
#endif
int CORE_cttqrt(int M, int N, int IB,
                PLASMA_Complex32_t *A1, int LDA1,
                PLASMA_Complex32_t *A2, int LDA2,
                PLASMA_Complex32_t *T, int LDT,
                PLASMA_Complex32_t *TAU, PLASMA_Complex32_t *WORK)
{
    static PLASMA_Complex32_t zone  = 1.0;
    static PLASMA_Complex32_t zzero = 0.0;
    static int                ione  = 1;

    PLASMA_Complex32_t alpha;
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

    for(ii = 0; ii < N; ii += IB) {
        sb = min(N-ii, IB);
        for(i = 0; i < sb; i++) {
            /*
             * Generate elementary reflector H( II*IB+I ) to annihilate
             * A( II*IB+I:mi, II*IB+I ).
             */
            mi = ii + i + 1;
            LAPACKE_clarfg_work(mi+1, &A1[LDA1*(ii+i)+ii+i], &A2[LDA2*(ii+i)], ione, &TAU[ii+i]);

            if (sb-i-1>0) {
                /*
                 * Apply H( II*IB+I ) to A( II*IB+I:M, II*IB+I+1:II*IB+IB ) from the left.
                 */
                ni = sb-i-1;
                cblas_ccopy(
                    ni,
                    &A1[LDA1*(ii+i+1)+(ii+i)], LDA1,
                    WORK, 1);

#ifdef COMPLEX
                LAPACKE_clacgv_work(ni, WORK, ione);
#endif
                cblas_cgemv(
                    CblasColMajor, (CBLAS_TRANSPOSE)PlasmaConjTrans,
                    mi, ni,
                    CBLAS_SADDR(zone), &A2[LDA2*(ii+i+1)], LDA2,
                    &A2[LDA2*(ii+i)], 1,
                    CBLAS_SADDR(zone), WORK, 1);
#ifdef COMPLEX
                LAPACKE_clacgv_work(ni, WORK, ione);
#endif
                alpha = -conjf(TAU[ii+i]);
                cblas_caxpy(
                    ni, CBLAS_SADDR(alpha),
                    WORK, 1,
                    &A1[LDA1*(ii+i+1)+ii+i], LDA1);
#ifdef COMPLEX
                LAPACKE_clacgv_work(ni, WORK, ione);
#endif
                cblas_cgerc(
                    CblasColMajor, mi, ni,
                    CBLAS_SADDR(alpha), &A2[LDA2*(ii+i)], 1,
                    WORK, 1,
                    &A2[LDA2*(ii+i+1)], LDA2);
            }

            /*
             * Calculate T.
             */

            if (i > 0 ) {

                cblas_ccopy(i, &A2[LDA2*(ii+i)+ii], 1, &WORK[ii], 1);
    
                cblas_ctrmv(
                    CblasColMajor, (CBLAS_UPLO)PlasmaUpper,
                    (CBLAS_TRANSPOSE)PlasmaConjTrans, (CBLAS_DIAG)PlasmaNonUnit,
                    i, &A2[LDA2*ii+ii], LDA2,
                    &WORK[ii], 1);
    
                alpha = -(TAU[ii+i]);
                for(j = 0; j < i; j++) {
                    WORK[ii+j] = alpha * WORK[ii+j];
                }
    
                if (ii > 0) {
                    cblas_cgemv(
                        CblasColMajor, (CBLAS_TRANSPOSE)PlasmaConjTrans, ii, i,
                        CBLAS_SADDR(alpha), &A2[LDA2*ii], LDA2,
                        &A2[LDA2*(ii+i)], 1,
                        CBLAS_SADDR(zzero), WORK, 1);
    
                    cblas_caxpy(i, CBLAS_SADDR(zone), &WORK[ii], 1, WORK, 1);
                }
    
                cblas_ccopy(i, WORK, 1, &T[LDT*(ii+i)], 1);
    
                cblas_ctrmv(
                    CblasColMajor, (CBLAS_UPLO)PlasmaUpper,
                    (CBLAS_TRANSPOSE)PlasmaNoTrans, (CBLAS_DIAG)PlasmaNonUnit,
                    i, &T[LDT*ii], LDT,
                    &T[LDT*(ii+i)], 1);

            }

            T[LDT*(ii+i)+i] = TAU[ii+i];
        }

        /* Apply Q' to the rest of the matrix to the left  */
        if (N > ii+sb) {
            CORE_cttrfb(
                PlasmaLeft, PlasmaConjTrans,
                PlasmaForward, PlasmaColumnwise,
                sb, N-(ii+sb), ii+sb, N-(ii+sb), sb,
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
void QUARK_CORE_cttqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       PLASMA_Complex32_t *A1, int lda1,
                       PLASMA_Complex32_t *A2, int lda2,
                       PLASMA_Complex32_t *T, int ldt)
{
    DAG_CORE_TTQRT;
    QUARK_Insert_Task(quark, CORE_cttqrt_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A1,            INOUT|QUARK_REGION_D|QUARK_REGION_U,
        sizeof(int),                        &lda1,  VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A2,            INOUT|QUARK_REGION_D|QUARK_REGION_U | LOCALITY,
        sizeof(int),                        &lda2,  VALUE,
        sizeof(PLASMA_Complex32_t)*ib*nb,    T,             OUTPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(PLASMA_Complex32_t)*nb,       NULL,          SCRATCH,
        sizeof(PLASMA_Complex32_t)*ib*nb,    NULL,          SCRATCH,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cttqrt_quark = PCORE_cttqrt_quark
#define CORE_cttqrt_quark PCORE_cttqrt_quark
#endif
void CORE_cttqrt_quark(Quark *quark)
{
    int m;
    int n;
    int ib;
    PLASMA_Complex32_t *A1;
    int lda1;
    PLASMA_Complex32_t *A2;
    int lda2;
    PLASMA_Complex32_t *T;
    int ldt;
    PLASMA_Complex32_t *TAU;
    PLASMA_Complex32_t *WORK;

    quark_unpack_args_11(quark, m, n, ib, A1, lda1, A2, lda2, T, ldt, TAU, WORK);
    CORE_cttqrt(m, n, ib, A1, lda1, A2, lda2, T, ldt, TAU, WORK);
}
