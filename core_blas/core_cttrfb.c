/**
 *
 * @file core_cttrfb.c
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
#include <cblas.h>
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex32_t
 *
 *  CORE_cttrfb applies a complex upper triangular block reflector H
 *  or its transpose H' to a
 *  complex rectangular matrix formed by coupling two tiles A1 and A2.
 *
 *******************************************************************************
 *
 * @param[in] side
 *         @arg PlasmaLeft  : apply Q or Q**H from the Left;
 *         @arg PlasmaRight : apply Q or Q**H from the Right.
 *
 * @param[in] trans
 *         @arg PlasmaNoTrans   : No transpose, apply Q;
 *         @arg PlasmaConjTrans : ConjTranspose, apply Q**H.
 *
 * @param[in] direct
 *         Indicates how H is formed from a product of elementary
 *         reflectors
 *         @arg PlasmaForward  : H = H(1) H(2) . . . H(k) (Forward)
 *         @arg PlasmaBackward : H = H(k) . . . H(2) H(1) (Backward)
 *
 * @param[in] storev
 *         Indicates how the vectors which define the elementary
 *         reflectors are stored:
 *         @arg PlasmaColumnwise
 *         @arg PlasmaRowwise
 *
 * @param[in] M1
 *         The number of columns of the tile A1. M1 >= 0.
 *
 * @param[in] N1
 *         The number of rows of the tile A1. N1 >= 0.
 *
 * @param[in] M2
 *         The number of columns of the tile A2. M2 >= 0.
 *
 * @param[in] N2
 *         The number of rows of the tile A2. N2 >= 0.
 *
 * @param[in] K
 *          The order of the matrix T (= the number of elementary
 *          reflectors whose product defines the block reflector).
 *
 * @param[in,out] A1
 *         On entry, the M1-by-N1 tile A1.
 *         On exit, A1 is overwritten by the application of Q.
 *
 * @param[in] LDA1
 *         The leading dimension of the array A1. LDA1 >= max(1,N1).
 *
 * @param[in,out] A2
 *         On entry, the M2-by-N2 tile A2.
 *         On exit, A2 is overwritten by the application of Q.
 *
 * @param[in] LDA2
 *         The leading dimension of the tile A2. LDA2 >= max(1,N2).
 *
 * @param[in] V
 *         (LDV,K) if STOREV = 'C'
 *         (LDV,M2) if STOREV = 'R' and SIDE = 'L'
 *         (LDV,N2) if STOREV = 'R' and SIDE = 'R'
 *         Matrix V.
 *
 * @param[in] LDV
 *         The leading dimension of the array V.
 *         If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M2);
 *         if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N2);
 *         if STOREV = 'R', LDV >= K.
 *
 * @param[out] T
 *         The triangular K-by-K matrix T in the representation of the
 *         block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= K.
 *
 * @param[in,out] WORK
 *
 * @param[in] LDWORK
 *         The dimension of the array WORK.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cttrfb = PCORE_cttrfb
#define CORE_cttrfb PCORE_cttrfb
#endif
int CORE_cttrfb(int side, int trans, int direct, int storev,
                int M1, int N1, int M2, int N2, int K,
                PLASMA_Complex32_t *A1, int LDA1,
                PLASMA_Complex32_t *A2, int LDA2,
                PLASMA_Complex32_t *V, int LDV,
                PLASMA_Complex32_t *T, int LDT,
                PLASMA_Complex32_t *WORK, int LDWORK)
{
    static PLASMA_Complex32_t zone  =  1.0;
    static PLASMA_Complex32_t mzone = -1.0;

    int j, vi;

    /* Check input arguments */
    if (M1 < 0) {
        coreblas_error(5, "Illegal value of M1");
        return -5;
    }
    if (N1 < 0) {
        coreblas_error(6, "Illegal value of N1");
        return -6;
    }
    if ((M2 < 0) ||
        ( (side == PlasmaRight) && (M1 != M2) ) ) {
        coreblas_error(7, "Illegal value of M2");
        return -7;
    }
    if ((N2 < 0) ||
        ( (side == PlasmaLeft) && (N1 != N2) ) ) {
        coreblas_error(8, "Illegal value of N2");
        return -8;
    }
    if (K < 0) {
        coreblas_error(9, "Illegal value of K");
        return -9;
    }

    /* Quick return */
    if ((M1 == 0) || (N1 == 0) || (M2 == 0) || (N2 == 0) || (K == 0))
        return PLASMA_SUCCESS;

    if (storev == PlasmaColumnwise) {

        if (direct == PlasmaForward) {
            /*
             * Let  V =  ( V1 )    (first K rows)
             *           ( V2 )
             * where V2 is non-unit upper triangular
             */
            if (side == PlasmaLeft) {

                /*
                 * Colwise / Forward / Left
                 * -------------------------
                 *
                 * Form  H * A  or  H' * A  where  A = ( A1 )
                 *                                     ( A2 )
                 * where A2 = ( A2_1 )
                 *            ( A2_2 )
                 */

                /*
                 * W = A1 + V' * A2
                 */

                /*
                 * W = A2_2
                 */
                LAPACKE_clacpy_work(LAPACK_COL_MAJOR,
                    lapack_const(PlasmaUpperLower),
                    K, N2,
                    &A2[M2-K], LDA2, WORK, LDWORK);

                /*
                 * W = V2' * A2_2
                 */
                cblas_ctrmm(
                    CblasColMajor, CblasLeft, CblasUpper,
                    CblasConjTrans, CblasNonUnit, K, N2,
                    CBLAS_SADDR(zone), &V[M2-K], LDV,
                    WORK, LDWORK);

                if (M2 > K) {
                    /*
                     * W = W + V1' * A2_1
                     */
                    cblas_cgemm(
                        CblasColMajor, CblasConjTrans, CblasNoTrans,
                        K, N2, M2-K,
                        CBLAS_SADDR(zone), V, LDV,
                        A2, LDA2,
                        CBLAS_SADDR(zone), WORK, LDWORK);
                }

                /*
                 * W = A1 + W
                 */
                for(j = 0; j < N1; j++) {
                    cblas_caxpy(K, CBLAS_SADDR(zone),
                            &A1[LDA1*j], 1,
                            &WORK[LDWORK*j], 1);
                }

                /*
                 * A2 = A2 - V * T * W -> W = T * W, A2 = A2 - V * W
                 */

                /*
                 * W = T * W
                 */
                cblas_ctrmm(
                    CblasColMajor, CblasLeft, CblasUpper,
                    (CBLAS_TRANSPOSE)trans, CblasNonUnit, K, N2,
                    CBLAS_SADDR(zone), T, LDT, WORK, LDWORK);
                /*
                 * A1 = A1 - W
                 */
                for(j = 0; j < N1; j++) {
                    cblas_caxpy(K, CBLAS_SADDR(mzone),
                            &WORK[LDWORK*j], 1,
                            &A1[LDA1*j], 1);
                }

                /*
                 * A2_1 = A2_1 - V1 * W
                 */
                if (M2 > K) {
                    cblas_cgemm(
                        CblasColMajor, CblasNoTrans, CblasNoTrans,
                        M2-K, N2, K,
                        CBLAS_SADDR(mzone), V, LDV,
                        WORK, LDWORK,
                        CBLAS_SADDR(zone), A2, LDA2);
                }

                /*
                 * W = - V2 * W
                 */
                cblas_ctrmm(
                    CblasColMajor, CblasLeft, CblasUpper,
                    CblasNoTrans, CblasNonUnit, K, N2,
                    CBLAS_SADDR(mzone), &V[M2-K], LDV,
                    WORK, LDWORK);

                /*
                 * A2_2 = A2_2 + W
                 */
                for(j = 0; j < N2; j++) {
                    cblas_caxpy(
                        K, CBLAS_SADDR(zone),
                        &WORK[LDWORK*j], 1,
                        &A2[LDA2*j+(M2-K)], 1);
                }
            }
            else {
                /* 
                 * Colwise / Forward / Right
                 * -------------------------
                 *
                 * Form  H * A  or  H' * A  where:
                 * 
                 *   A  = ( A1 A2 )
                 *
                 *   A2 = ( A2_1 : A2_2 )
                 *
                 *         A2_1 is M2 x (M2-K)
                 *         A2_2 is M2 x K
                 *
                 *   V = ( V_1 ) 
                 *       ( V_2 )
                 *
                 *         V_1 is full and (N2-K) x K
                 *         V_2 is upper triangular and K x K
                 */

                /*
                 * W = ( A1 + A2_1*V_1 + A2_2*V_2 ) * op(T)
                 *
                 *    W  is M x K
                 *    A1 is M x K
                 *    A2 is M x N2 split as (A2_1 A2_2) such as
                 *       A2_1 is (N2-K) x K
                 *       A2_2 is M x K
                 */

                /* W = A2_2 */
                LAPACKE_clacpy_work(LAPACK_COL_MAJOR,
                    lapack_const(PlasmaUpperLower),
                    M2, K,
                    &A2[LDA2*(N2-K)], LDA2, WORK, LDWORK);

                /* W = W * V_2 --> W = A2_2 * V_2 */
                cblas_ctrmm(
                    CblasColMajor, CblasRight, CblasUpper,
                    CblasNoTrans, CblasNonUnit, M2, K,
                    CBLAS_SADDR(zone), &V[N2-K], LDV,
                    WORK, LDWORK);

                /* W = W + A2_1 * V_1 */
                if (N2 > K) {
                    cblas_cgemm(
                        CblasColMajor, CblasNoTrans, CblasNoTrans,
                        M2, K, N2-K,
                        CBLAS_SADDR(zone), A2, LDA2,
                        V, LDV,
                        CBLAS_SADDR(zone), WORK, LDWORK);
                }

                /* W = A1 + W */
                for (j = 0; j < K; j++) {
                    cblas_caxpy(M1, CBLAS_SADDR(zone), &A1[LDA1*j], 1,
                                &WORK[LDWORK*j], 1);
                }

                /* W = W * T --> ( A1 + A2_1*V_1 + A2_2*V_2 ) * op(T) */
                cblas_ctrmm(
                    CblasColMajor, CblasRight, CblasUpper,
                    (CBLAS_TRANSPOSE)trans, CblasNonUnit, M2, K,
                    CBLAS_SADDR(zone), T, LDT, WORK, LDWORK);

                /*
                 * A1 = A1 - W
                 */

                for(j = 0; j < K; j++) {
                    cblas_caxpy(M1, CBLAS_SADDR(mzone), &WORK[LDWORK*j], 1,
                                &A1[LDA1*j], 1);
                }

                /*
                 * A2 = A2 - W * V' --> A2 - W*V_1' - W*V_2'
                 */

                /* A2 = A2 - W * V_1' */
                if (N2 > K) {
                    cblas_cgemm(
                        CblasColMajor, CblasNoTrans, CblasConjTrans,
                        M2, N2-K, K,
                        CBLAS_SADDR(mzone), WORK, LDWORK,
                        V, LDV,
                        CBLAS_SADDR(zone), A2, LDA2);
                }

                /* A2 =  A2 -  W * V_2' */
                cblas_ctrmm(
                    CblasColMajor, CblasRight, CblasUpper,
                    CblasConjTrans, CblasNonUnit, M2, K,
                    CBLAS_SADDR(mzone), &V[N2-K], LDV,
                    WORK, LDWORK);

                for(j = 0; j < K; j++) {
                    cblas_caxpy(
                        M2, CBLAS_SADDR(zone),
                        &WORK[LDWORK*j], 1,
                        &A2[LDA2*(j+N2-K)], 1);
                }
            }
        }
        else {
            coreblas_error(3, "Not implemented (ColWise / Backward / Left or Right)");
            return PLASMA_ERR_NOT_SUPPORTED;
        }
    }
    else {
        /*
         * Rowwise
         */
        if (direct == PlasmaForward) {
            /*
             * Let  V =  ( V1 V2 )    (V1: first K cols)
             * 
             * where V2 is non-unit lower triangular
             */

            if (side == PlasmaLeft) {

                /*
                 * Rowwise / Forward / Left
                 * -------------------------
                 *
                 * Form  H * A  or  H' * A  where  A = ( A1 )
                 *                                     ( A2 )
                 * where A2 = ( A2_1 )
                 *            ( A2_2 )
                 */

                /* V_2 first element */
                vi = LDV*(M2-K);

                /*
                 * W = A1 + V * A2
                 */

                /*
                 * W = A2_2
                 */
                LAPACKE_clacpy_work(LAPACK_COL_MAJOR,
                    lapack_const(PlasmaUpperLower),
                    K, N2,
                    &A2[M2-K], LDA2, WORK, LDWORK);

                /*
                 * W = V2 * A2_2
                 */
                //**DB CblasColMajor, CblasLeft, CblasUpper,
                cblas_ctrmm(
                    CblasColMajor, CblasLeft, CblasLower,
                    CblasNoTrans, CblasNonUnit, K, N2,
                    CBLAS_SADDR(zone), &V[vi], LDV,
                    WORK, LDWORK);

                if (M2 > K) {
                    /*
                     * W = W + V1 * A2_1
                     */
                    cblas_cgemm(
                        CblasColMajor, CblasNoTrans, CblasNoTrans,
                        K, N2, M2-K,
                        CBLAS_SADDR(zone), V, LDV,
                        A2, LDA2,
                        CBLAS_SADDR(zone), WORK, LDWORK);
                }

                /*
                 * W = A1 + W
                 */
                for(j = 0; j < N1; j++) {
                    cblas_caxpy(
                            K, CBLAS_SADDR(zone),
                            &A1[LDA1*j], 1,
                            &WORK[LDWORK*j], 1);
                }

                /*
                 * W = T * W
                 */
                cblas_ctrmm(
                    CblasColMajor, CblasLeft, CblasUpper,
                    (CBLAS_TRANSPOSE)trans, CblasNonUnit, K, N2,
                    CBLAS_SADDR(zone), T, LDT, WORK, LDWORK);

                /*
                 * A1 = A1 - W
                 */
                for(j = 0; j < N1; j++) {
                    cblas_caxpy(
                            K, CBLAS_SADDR(mzone),
                            &WORK[LDWORK*j], 1,
                            &A1[LDA1*j], 1);
                }

                /*
                 * A2 = A2 - V' * T * W -> A2 = A2 - V' * W
                 */

                /*
                 * A2_1 = A2_1 - V1' * W
                 */
                if (M2 > K) {
                    cblas_cgemm(
                        CblasColMajor, CblasConjTrans, CblasNoTrans,
                        M2-K, N2, K,
                        CBLAS_SADDR(mzone), V, LDV,
                        WORK, LDWORK,
                        CBLAS_SADDR(zone), A2, LDA2);
                }

                /*
                 * W = - V2' * W
                 */
                cblas_ctrmm(
                    CblasColMajor, CblasLeft, CblasLower,
                    CblasConjTrans, CblasNonUnit, K, N2,
                    CBLAS_SADDR(mzone), &V[vi], LDV,
                    WORK, LDWORK);

                /*
                 * A2_2 = A2_2 + W
                 */
                for(j = 0; j < N2; j++) {
                    cblas_caxpy(
                        K, CBLAS_SADDR(zone),
                        &WORK[LDWORK*j], 1,
                        &A2[LDA2*j+(M2-K)], 1);
                }
            }
            else {
                /* 
                 * Rowwise / Forward / Right
                 * -------------------------
                 *
                 * Form  H * A  or  H' * A  where:
                 * 
                 *   A  = ( A1 A2 )
                 *
                 *   A2 = ( A2_1 : A2_2 )
                 *
                 *         A2_1 is M2 x (M2-K)
                 *         A2_2 is M2 x K
                 *
                 *   V = ( V_1 ) 
                 *       ( V_2 )
                 *
                 *         V_1 is full and (N2-K) x K
                 *         V_2 is lower triangular and K x K
                 */

                /*
                 * W = ( A1 + A2_1*V_1 + A2_2*V_2 ) * op(T)
                 *
                 *    W  is M x K
                 *    A1 is M x K
                 *    A2 is M x N2 split as (A2_1 A2_2) such as
                 *       A2_1 is (N2-K) x K
                 *       A2_2 is M x K
                 */

                /* V_2 and A2_2 first element */
                vi = LDV*(N2-K);

                /* W = A2_2 */
                LAPACKE_clacpy_work(LAPACK_COL_MAJOR,
                    lapack_const(PlasmaUpperLower),
                    M2, K,
                    &A2[LDA2*(N2-K)], LDA2, WORK, LDWORK);

                /* W = W * V_2' --> W = A2_2 * V_2' */
                cblas_ctrmm(
                    CblasColMajor, CblasRight, CblasLower,
                    CblasConjTrans, CblasNonUnit, M2, K,
                    CBLAS_SADDR(zone), &V[vi], LDV,
                    WORK, LDWORK);

                /* W = W + A2_1 * V_1' */
                if (N2 > K) {
                    cblas_cgemm(
                        CblasColMajor, CblasNoTrans, CblasConjTrans,
                        M2, K, N2-K,
                        CBLAS_SADDR(zone), A2, LDA2,
                        V, LDV,
                        CBLAS_SADDR(zone), WORK, LDWORK);
                }

                /* W = A1 + W */
                for (j = 0; j < K; j++) {
                    cblas_caxpy(M1, CBLAS_SADDR(zone), &A1[LDA1*j], 1,
                                &WORK[LDWORK*j], 1);
                }

                /* W = W * op(T) --> ( A1 + A2_1*V_1 + A2_2*V_2 ) * op(T) */
                cblas_ctrmm(
                    CblasColMajor, CblasRight, CblasUpper,
                    (CBLAS_TRANSPOSE)trans, CblasNonUnit, M2, K,
                    CBLAS_SADDR(zone), T, LDT, WORK, LDWORK);

                /*
                 * A1 = A1 - W
                 */

                for(j = 0; j < K; j++) {
                    cblas_caxpy(M1, CBLAS_SADDR(mzone), &WORK[LDWORK*j], 1,
                                &A1[LDA1*j], 1);
                }

                /*
                 * A2 = A2 - W * V --> A2 - W*V_1 - W*V_2
                 */

                /* A2 = A2 - W * V_1 */
                if (N2 > K) {
                    cblas_cgemm(
                        CblasColMajor, CblasNoTrans, CblasNoTrans,
                        M2, N2-K, K,
                        CBLAS_SADDR(mzone), WORK, LDWORK,
                        V, LDV,
                        CBLAS_SADDR(zone), A2, LDA2);
                }

                /* A2 =  A2 -  W * V_2 */
                cblas_ctrmm(
                    CblasColMajor, CblasRight, CblasLower,
                    CblasNoTrans, CblasNonUnit, M2, K,
                    CBLAS_SADDR(mzone), &V[vi], LDV,
                    WORK, LDWORK);

                for(j = 0; j < K; j++) {
                    cblas_caxpy(
                        M2, CBLAS_SADDR(zone),
                        &WORK[LDWORK*j], 1,
                        &A2[LDA2*(N2-K+j)], 1);
                }

            }
 
        }
        else {
            coreblas_error(3, "Not implemented (Rowwise / Backward / Left or Right)");
            return PLASMA_ERR_NOT_SUPPORTED;
        }
    }
    return PLASMA_SUCCESS;
}
