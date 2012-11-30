/**
 *
 * @file core_zherfb.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Hatem Ltaief
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
 *  CORE_zherfb overwrites the symmetric complex N-by-N tile C with
 *
 *    Q**T*C*Q
 *
 *  where Q is a complex unitary matrix defined as the product of k
 *  elementary reflectors
 *
 *    Q = H(1) H(2) . . . H(k)
 *
 *  as returned by CORE_zgeqrt. Only PlasmaLower supported!
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *         @arg PlasmaLower : the upper part of the symmetric matrix C 
 *                            is not referenced.
 *         @arg PlasmaUpper : the lower part of the symmetric matrix C 
 *                            is not referenced (not supported).
 *
 * @param[in] n
 *         The number of rows/columns of the tile C.  N >= 0.
 *
 * @param[in] k
 *         The number of elementary reflectors whose product defines
 *         the matrix Q. K >= 0.
 *
 * @param[in] ib
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in] nb
 *         The blocking size.  NB >= 0.
 *
 * @param[in] A
 *         The i-th column must contain the vector which defines the
 *         elementary reflector H(i), for i = 1,2,...,k, as returned by
 *         CORE_zgeqrt in the first k columns of its array argument A.
 *
 * @param[in] lda
 *         The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[out] T
 *         The IB-by-K triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] ldt
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[in,out] C
 *         On entry, the symmetric N-by-N tile C.
 *         On exit, C is overwritten by Q**T*C*Q.
 *
 * @param[in] ldc
 *         The leading dimension of the array C. LDC >= max(1,M).
 *
 * @param[in,out] WORK
 *         On exit, if INFO = 0, WORK(1) returns the optimal LDWORK.
 *
 * @param[in] ldwork
 *         The dimension of the array WORK. LDWORK >= max(1,N);
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zherfb = PCORE_zherfb
#define CORE_zherfb PCORE_zherfb
#define CORE_zunmlq PCORE_zunmlq
#define CORE_zunmqr PCORE_zunmqr
int  CORE_zunmlq(int side, int trans,
                 int M, int N, int IB, int K,
                 PLASMA_Complex64_t *V, int LDV,
                 PLASMA_Complex64_t *T, int LDT,
                 PLASMA_Complex64_t *C, int LDC,
                 PLASMA_Complex64_t *WORK, int LDWORK);
int  CORE_zunmqr(int side, int trans,
                 int M, int N, int K, int IB,
                 PLASMA_Complex64_t *V, int LDV,
                 PLASMA_Complex64_t *T, int LDT,
                 PLASMA_Complex64_t *C, int LDC,
                 PLASMA_Complex64_t *WORK, int LDWORK);
#endif
int CORE_zherfb( PLASMA_enum uplo, int n,
                 int k, int ib, int nb,
                 PLASMA_Complex64_t *A, int lda,
                 PLASMA_Complex64_t *T, int ldt,
                 PLASMA_Complex64_t *C, int ldc,
                 PLASMA_Complex64_t *WORK, int ldwork )
{
    int i, j;

    if (uplo == PlasmaLower) {
        /* Rebuild the symmetric block: WORK <- C */
        for (j = 0; j < n; j++)
            for (i = j; i < n; i++){
                *(WORK + i + j * ldwork) = *(C + i + j*ldc);
                if (i > j){
                    *(WORK + j + i * ldwork) =  *(WORK + i + j * ldwork);
#ifdef COMPLEX
                    LAPACKE_zlacgv_work(1, WORK + j + i * ldwork, ldwork);
#endif
                }
            }
        
        /* Left */
        CORE_zunmqr(PlasmaLeft, PlasmaConjTrans, n, n, k, ib, 
                    A, lda, T, ldt, WORK, ldwork, WORK+nb*ldwork, ldwork);
        /* Right */
        CORE_zunmqr(PlasmaRight, PlasmaNoTrans, n, n, k, ib, 
                    A, lda, T, ldt, WORK, ldwork, WORK+nb*ldwork, ldwork);
        
        /* 
         * Copy back the final result to the lower part of C 
         */
        /* C = WORK */
        for (j = 0; j < n; j++)
            for (i = j; i < n; i++)
                *(C + i + j*ldc) = *(WORK + i + j * ldwork);
    }
    else {
        /* Rebuild the symmetric block: WORK <- C */
        for (i = 0; i < n; i++)
            for (j = i; j < n; j++){
                *(WORK + i + j * ldwork) = *(C + i + j*ldc);
                if (j > i){
                    *(WORK + j + i * ldwork) =  *(WORK + i + j * ldwork);
#ifdef COMPLEX
                    LAPACKE_zlacgv_work(1, WORK + j + i * ldwork, ldwork);
#endif
                }
            }
        
        /* Right */
        CORE_zunmlq(PlasmaRight, PlasmaConjTrans, n, n, k, ib, 
                    A, lda, T, ldt, WORK, ldwork, WORK+nb*ldwork, ldwork);
        /* Left */
        CORE_zunmlq(PlasmaLeft, PlasmaNoTrans, n, n, k, ib, 
                    A, lda, T, ldt, WORK, ldwork, WORK+nb*ldwork, ldwork);
        
        /* 
         * Copy back the final result to the upper part of C 
         */
        /* C = WORK */
        for (i = 0; i < n; i++)
            for (j = i; j < n; j++)
                *(C + i + j*ldc) = *(WORK + i + j * ldwork);
    }
    return 0;
}


/***************************************************************************//**
 * This kernel is just a workaround for now... will be deleted eventually
 * and replaced by the one above (Piotr's Task)
 **/
void QUARK_CORE_zherfb(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo,
                       int n, int k, int ib, int nb,
                       PLASMA_Complex64_t *A, int lda,
                       PLASMA_Complex64_t *T, int ldt,
                       PLASMA_Complex64_t *C, int ldc)
{
    QUARK_Insert_Task(
        quark, CORE_zherfb_quark, task_flags,
        sizeof(PLASMA_enum),                     &uplo,  VALUE,
        sizeof(int),                             &n,     VALUE,
        sizeof(int),                             &k,     VALUE,
        sizeof(int),                             &ib,    VALUE,
        sizeof(int),                             &nb,    VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,        A,          uplo == PlasmaUpper ? INOUT|QUARK_REGION_U : INOUT|QUARK_REGION_L,
        sizeof(int),                             &lda,   VALUE,
        sizeof(PLASMA_Complex64_t)*ib*nb,        T,          INPUT,
        sizeof(int),                             &ldt,   VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,        C,          uplo == PlasmaUpper ? INOUT|QUARK_REGION_D|QUARK_REGION_U : INOUT|QUARK_REGION_D|QUARK_REGION_L,
        sizeof(int),                             &ldc,   VALUE,
        sizeof(PLASMA_Complex64_t)*2*nb*nb,    NULL,         SCRATCH,
        sizeof(int),                             &nb,    VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zherfb_quark = PCORE_zherfb_quark
#define CORE_zherfb_quark PCORE_zherfb_quark
#endif
void CORE_zherfb_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int n;
    int k;
    int ib;
    int nb;
    PLASMA_Complex64_t *A;
    int lda;
    PLASMA_Complex64_t *T;
    int ldt;
    PLASMA_Complex64_t *C;
    int ldc;
    PLASMA_Complex64_t *WORK;
    int ldwork;

    quark_unpack_args_13(quark, uplo, n, k, ib, nb, A, lda, T, ldt, C, ldc, WORK, ldwork);
    CORE_zherfb(uplo, n, k, ib, nb, A, lda, T, ldt, C, ldc, WORK, ldwork);
}
