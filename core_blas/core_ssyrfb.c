/**
 *
 * @file core_ssyrfb.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @generated s Thu Sep 15 12:09:01 2011
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
 *  CORE_ssyrfb overwrites the symmetric complex N-by-N tile C with
 *
 *    Q**T*C*Q
 *
 *  where Q is a complex unitary matrix defined as the product of k
 *  elementary reflectors
 *
 *    Q = H(1) H(2) . . . H(k)
 *
 *  as returned by CORE_sgeqrt. Only PlasmaLower supported!
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
 *         CORE_sgeqrt in the first k columns of its array argument A.
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
#pragma weak CORE_ssyrfb = PCORE_ssyrfb
#define CORE_ssyrfb PCORE_ssyrfb
#define CORE_sormlq PCORE_sormlq
#define CORE_sormqr PCORE_sormqr
int  CORE_sormlq(int side, int trans,
                 int M, int N, int IB, int K,
                 float *V, int LDV,
                 float *T, int LDT,
                 float *C, int LDC,
                 float *WORK, int LDWORK);
int  CORE_sormqr(int side, int trans,
                 int M, int N, int K, int IB,
                 float *V, int LDV,
                 float *T, int LDT,
                 float *C, int LDC,
                 float *WORK, int LDWORK);
#endif
int CORE_ssyrfb( PLASMA_enum uplo, int n,
                 int k, int ib, int nb,
                 float *A, int lda,
                 float *T, int ldt,
                 float *C, int ldc,
                 float *WORK, int ldwork )
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
                    LAPACKE_slacgv_work(1, WORK + j + i * ldwork, ldwork);
#endif
                }
            }
        
        /* Left */
        CORE_sormqr(PlasmaLeft, PlasmaTrans, n, n, k, ib, 
                    A, lda, T, ldt, WORK, ldwork, WORK+nb*ldwork, ldwork);
        /* Right */
        CORE_sormqr(PlasmaRight, PlasmaNoTrans, n, n, k, ib, 
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
                    LAPACKE_slacgv_work(1, WORK + j + i * ldwork, ldwork);
#endif
                }
            }
        
        /* Right */
        CORE_sormlq(PlasmaRight, PlasmaTrans, n, n, k, ib, 
                    A, lda, T, ldt, WORK, ldwork, WORK+nb*ldwork, ldwork);
        /* Left */
        CORE_sormlq(PlasmaLeft, PlasmaNoTrans, n, n, k, ib, 
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
void QUARK_CORE_ssyrfb(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo,
                       int n, int k, int ib, int nb,
                       float *A, int lda,
                       float *T, int ldt,
                       float *C, int ldc)
{
    QUARK_Insert_Task(
        quark, CORE_ssyrfb_quark, task_flags,
        sizeof(PLASMA_enum),                     &uplo,  VALUE,
        sizeof(int),                             &n,     VALUE,
        sizeof(int),                             &k,     VALUE,
        sizeof(int),                             &ib,    VALUE,
        sizeof(int),                             &nb,    VALUE,
        sizeof(float)*nb*nb,        A,          uplo == PlasmaUpper ? INOUT|QUARK_REGION_U : INOUT|QUARK_REGION_L,
        sizeof(int),                             &lda,   VALUE,
        sizeof(float)*ib*nb,        T,          INPUT,
        sizeof(int),                             &ldt,   VALUE,
        sizeof(float)*nb*nb,        C,          uplo == PlasmaUpper ? INOUT|QUARK_REGION_D|QUARK_REGION_U : INOUT|QUARK_REGION_D|QUARK_REGION_L,
        sizeof(int),                             &ldc,   VALUE,
        sizeof(float)*2*nb*nb,    NULL,         SCRATCH,
        sizeof(int),                             &nb,    VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ssyrfb_quark = PCORE_ssyrfb_quark
#define CORE_ssyrfb_quark PCORE_ssyrfb_quark
#endif
void CORE_ssyrfb_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int n;
    int k;
    int ib;
    int nb;
    float *A;
    int lda;
    float *T;
    int ldt;
    float *C;
    int ldc;
    float *WORK;
    int ldwork;

    quark_unpack_args_13(quark, uplo, n, k, ib, nb, A, lda, T, ldt, C, ldc, WORK, ldwork);
    CORE_ssyrfb(uplo, n, k, ib, nb, A, lda, T, ldt, C, ldc, WORK, ldwork);
}
