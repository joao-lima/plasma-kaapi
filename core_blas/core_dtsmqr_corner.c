/**
 *
 * @file core_dtsmqr_corner.c
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
 *  CORE_dtsmqr_corner: see CORE_dtsmqr
 *
 * This kernel applies left and right transformations as depicted below:
 * |I -VT'V'| * | A1 A2'| * |I - VTV'|
 *              | A2 A3 |
 * where A1 and A3 are symmetric matrices.
 * Only the lower part is referenced.
 * This is an adhoc implementation, can be further optimized...
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
 *         CORE_DTSQRT in the first k columns of its array argument V.
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
 *             LDWORK-by-N1 if side == PlasmaLeft
 *             LDWORK-by-IB if side == PlasmaRight
 *
 * @param[in] LDWORK
 *         The leading dimension of the array WORK. 
 *             LDWORK >= max(1,IB) if side == PlasmaLeft
 *             LDWORK >= max(1,M1) if side == PlasmaRight
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dtsmqr_corner = PCORE_dtsmqr_corner
#define CORE_dtsmqr_corner PCORE_dtsmqr_corner
#define CORE_dtsmqr PCORE_dtsmqr
int  CORE_dtsmqr(int side, int trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 double *A1, int LDA1,
                 double *A2, int LDA2,
                 double *V, int LDV,
                 double *T, int LDT,
                 double *WORK, int LDWORK);
#endif
int CORE_dtsmqr_corner( int m1, int n1, int m2, int n2, int m3, int n3,
                        int k, int ib, int nb,
                        double *A1, int lda1,
                        double *A2, int lda2,
                        double *A3, int lda3,
                        double *V, int ldv,
                        double *T, int ldt,
                        double *WORK, int ldwork)
{
    int i, j;
    PLASMA_enum side, trans;

    if ( m1 != n1 ) {
        coreblas_error(1, "Illegal value of M1, N1");
        return -1;
    }

    /*  Rebuild the symmetric block: WORK <- A1 */
    for (j = 0; j < n1; j++)
        for (i = j; i < m1; i++){
            *(WORK + i + j*ldwork) = *(A1 + i + j*lda1);
            if (i > j){
                *(WORK + j + i*ldwork) =  ( *(WORK + i + j*ldwork) );
            }
        }
    
    /*  Copy the transpose of A2: WORK+nb*ldwork <- A2' */
    for (j = 0; j < n2; j++)
        for (i = 0; i < m2; i++){
            *(WORK + j + (i + nb) * ldwork) = ( *(A2 + i + j*lda2) );
        }

    side  = PlasmaLeft;
    trans = PlasmaTrans;

    /*  Left application on |A1| */
    /*                      |A2| */
    CORE_dtsmqr(side, trans, m1, n1, m2, n2, k, ib, 
                WORK, ldwork, A2, lda2, 
                V, ldv, T, ldt, 
                WORK + 3*nb*ldwork, ldwork);

    /*  Rebuild the symmetric block: WORK+2*nb*ldwork <- A3 */
    for (j = 0; j < n3; j++)
        for (i = j; i < m3; i++){
            *(WORK + i + (j + 2*nb) * ldwork) = *(A3 + i + j*lda3);
            if (i != j){
                *(WORK + j + (i + 2*nb) * ldwork) =  ( *(WORK + i + (j + 2*nb) * ldwork) );
            }
        }
    /*  Left application on | A2'| */
    /*                      | A3 | */
    CORE_dtsmqr(side, trans, n2, m2, m3, n3, k, ib, 
                WORK+nb*ldwork, ldwork, WORK+2*nb*ldwork, ldwork, 
                V, ldv, T, ldt, 
                WORK + 3*nb*ldwork, ldwork);

    side  = PlasmaRight;
    trans = PlasmaNoTrans;

    /*  Right application on | A1 A2' | */
    CORE_dtsmqr(side, trans, m1, n1, n2, m2, k, ib, 
                WORK, ldwork, WORK+nb*ldwork, ldwork, 
                V, ldv, T, ldt, 
                WORK + 3*nb*ldwork, ldwork);

    /*  Copy back the final result to the lower part of A1 */
    /*  A1 = WORK */
    for (j = 0; j < n1; j++)
        for (i = j; i < m1; i++)
            *(A1 + i + j*lda1) = *(WORK + i + j*ldwork);

    /*  Right application on | A2 A3 | */
    CORE_dtsmqr(side, trans, m2, n2, m3, n3, k, ib, 
                A2, lda2, WORK+2*nb*ldwork, ldwork, 
                V,  ldv,  T, ldt, 
                WORK + 3*nb*ldwork, ldwork);

    /*  Copy back the final result to the lower part of A3 */
    /*  A3 = WORK+2*nb*ldwork */
    for (j = 0; j < n3; j++)
        for (i = j; i < m3; i++)
            *(A3 + i + j*lda3) = *(WORK + i + (j+ 2*nb) * ldwork);

    return PLASMA_SUCCESS;
}
/***************************************************************************//**
 *
 **/
void QUARK_CORE_dtsmqr_corner(Quark *quark, Quark_Task_Flags *task_flags,
                         int m1, int n1, int m2, int n2, int m3, int n3, int k, int ib, int nb,
                         double *A1, int lda1,
                         double *A2, int lda2,
                         double *A3, int lda3,
                         double *V, int ldv,
                         double *T, int ldt)
{
    int ldwork = nb;

    QUARK_Insert_Task(quark, CORE_dtsmqr_corner_quark, task_flags,
        sizeof(int),                        &m1,    VALUE,
        sizeof(int),                        &n1,    VALUE,
        sizeof(int),                        &m2,    VALUE,
        sizeof(int),                        &n2,    VALUE,
        sizeof(int),                        &m3,    VALUE,
        sizeof(int),                        &n3,    VALUE,
        sizeof(int),                        &k,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(int),                        &nb,    VALUE,
        sizeof(double)*nb*nb,    A1,            INOUT|QUARK_REGION_D|QUARK_REGION_L,
        sizeof(int),                        &lda1,  VALUE,
        sizeof(double)*nb*nb,    A2,            INOUT,
        sizeof(int),                        &lda2,  VALUE,
        sizeof(double)*nb*nb,    A3,            INOUT|QUARK_REGION_D|QUARK_REGION_L,
        sizeof(int),                        &lda3,  VALUE,
        sizeof(double)*nb*nb,    V,             INPUT,
        sizeof(int),                        &ldv,   VALUE,
        sizeof(double)*ib*nb,    T,             INPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(double)*4*nb*nb,    NULL,          SCRATCH,
        sizeof(int),                        &ldwork, VALUE,
        0);
}


#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dtsmqr_corner_quark = PCORE_dtsmqr_corner_quark
#define CORE_dtsmqr_corner_quark PCORE_dtsmqr_corner_quark
#endif
void CORE_dtsmqr_corner_quark(Quark *quark)
{
    int m1;
    int n1;
    int m2;
    int n2;
    int m3;
    int n3;
    int k;
    int ib;
    int nb;
    double *A1;
    int lda1;
    double *A2;
    int lda2;
    double *A3;
    int lda3;
    double *V;
    int ldv;
    double *T;
    int ldt;
    double *WORK;
    int ldwork;

    quark_unpack_args_21(quark, m1, n1, m2, n2, m3, n3, k, ib, nb, 
                         A1, lda1, A2, lda2, A3, lda3, V, ldv, T, ldt, WORK, ldwork);
    CORE_dtsmqr_corner(m1, n1, m2, n2, m3, n3, k, ib, nb, 
                       A1, lda1, A2, lda2, A3, lda3, V, ldv, T, ldt, WORK, ldwork);
}
