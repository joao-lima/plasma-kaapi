/**
 *
 * @file core_claswp.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated c Thu Sep 15 12:09:00 2011
 *
 **/
#include <lapacke.h>
#include "common.h"

#define A(m, n) BLKADDR(descA, PLASMA_Complex32_t, m, n)

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex32_t
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_claswp = PCORE_claswp
#define CORE_claswp PCORE_claswp
#endif
void CORE_claswp(int N, PLASMA_Complex32_t *A, int LDA, int I1,  int I2, int *IPIV, int INC)
{
    LAPACKE_claswp_work( LAPACK_COL_MAJOR, N, A, LDA, I1, I2, IPIV, INC );
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_claswp(Quark *quark, Quark_Task_Flags *task_flags,
                       int n, PLASMA_Complex32_t *A, int lda, 
                       int i1,  int i2, int *ipiv, int inc)
{
    DAG_CORE_LASWP;
    QUARK_Insert_Task(
        quark, CORE_claswp_quark, task_flags,
        sizeof(int),                      &n,    VALUE,
        sizeof(PLASMA_Complex32_t)*lda*n,  A,        INOUT | LOCALITY,
        sizeof(int),                      &lda,  VALUE,
        sizeof(int),                      &i1,   VALUE,
        sizeof(int),                      &i2,   VALUE,
        sizeof(int)*n,                     ipiv,     INPUT,
        sizeof(int),                      &inc,  VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_claswp_quark = PCORE_claswp_quark
#define CORE_claswp_quark PCORE_claswp_quark
#endif
void CORE_claswp_quark(Quark *quark)
{
    int n, lda, i1, i2, inc;
    int *ipiv;
    PLASMA_Complex32_t *A;
    
    quark_unpack_args_7(quark, n, A, lda, i1, i2, ipiv, inc);
    LAPACKE_claswp_work(LAPACK_COL_MAJOR, n, A, lda, i1, i2, ipiv, inc );
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_claswp_f2(Quark *quark, Quark_Task_Flags *task_flags,
                          int n, PLASMA_Complex32_t *A, int lda, 
                          int i1,  int i2, int *ipiv, int inc,
                          PLASMA_Complex32_t *fake1, int szefake1, int flag1,
                          PLASMA_Complex32_t *fake2, int szefake2, int flag2)
{
    DAG_CORE_LASWP;
    QUARK_Insert_Task(
        quark, CORE_claswp_f2_quark, task_flags,
        sizeof(int),                        &n,     VALUE,
        sizeof(PLASMA_Complex32_t)*lda*n,    A,         INOUT | LOCALITY,
        sizeof(int),                        &lda,   VALUE,
        sizeof(int),                        &i1,    VALUE,
        sizeof(int),                        &i2,    VALUE,
        sizeof(int)*n,                       ipiv,      INPUT,
        sizeof(int),                        &inc,   VALUE,
        sizeof(PLASMA_Complex32_t)*szefake1, fake1,     flag1,
        sizeof(PLASMA_Complex32_t)*szefake2, fake2,     flag2,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_claswp_f2_quark = PCORE_claswp_f2_quark
#define CORE_claswp_f2_quark PCORE_claswp_f2_quark
#endif
void CORE_claswp_f2_quark(Quark* quark)
{
    int n, lda, i1, i2, inc;
    int *ipiv;
    PLASMA_Complex32_t *A;
    void *fake1, *fake2;
    
    quark_unpack_args_9(quark, n, A, lda, i1, i2, ipiv, inc, fake1, fake2);
    LAPACKE_claswp_work(LAPACK_COL_MAJOR, n, A, lda, i1, i2, ipiv, inc );
}

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex32_t
 *
 *  CORE_claswp_ontile apply the claswp function on a matrix stored in
 *  tile layout
 *
 *******************************************************************************
 *
 *  @param[in,out] A
 *          The descriptor of the matrix A to permute.
 *
 *  @param[in] i1
 *          The first element of IPIV for which a row interchange will
 *          be done.
 *
 *  @param[in] i2
 *          The last element of IPIV for which a row interchange will
 *          be done.
 *
 *  @param[in] ipiv
 *          The pivot indices; Only the element in position i1 to i2
 *          are accessed. The pivot are offset by A.i.
 *
 *  @param[in] inc
 *          The increment between successive values of IPIV.  If IPIV
 *          is negative, the pivots are applied in reverse order.
 *
 *******************************************************************************
 */
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_claswp_ontile = PCORE_claswp_ontile
#define CORE_claswp_ontile PCORE_claswp_ontile
#endif
int CORE_claswp_ontile(PLASMA_desc descA, int i1, int i2, int *ipiv, int inc)
{
    int i, j, ip, it;
    PLASMA_Complex32_t *A1;
    int lda1, lda2;

    /* Change i1 to C notation */
    i1--;
    if ( descA.nt > 1 ) {
        coreblas_error(1, "Illegal value of descA.nt");
        return -1;
    }
    if ( i1 < 0 ) {
        coreblas_error(2, "Illegal value of i1");
        return -2;
    }
    if ( (i2 < i1) || (i2 > descA.m) ) {
        coreblas_error(3, "Illegal value of i2");
        return -3;
    }
    if ( ! ( (i2 - i1 - i1%descA.mb -1) < descA.mb ) ) {
        coreblas_error(2, "Illegal value of i1,i2. They have to be part of the same block.");
        return -3;
    }

    it = i1 / descA.mb;
    if (inc > 0) {
        A1 = A(it, 0);
        lda1 = BLKLDD(descA, 0);

        for (j = i1; j < i2; ++j, ipiv+=inc) {
            ip = (*ipiv) - descA.i - 1;
            if ( ip != j )
            {
                it = ip / descA.mb;
                i  = ip % descA.mb;
                lda2 = BLKLDD(descA, it);
                cblas_cswap(descA.n, A1       + j, lda1,
                                     A(it, 0) + i, lda2 );
            }
        }   
    }
    else
    {
        A1 = A(it, 0);
        lda1 = BLKLDD(descA, descA.mt-1);
        
        i1--;
        ipiv = &ipiv[(1-i2)*inc];
        for (j = i2-1; j > i1; --j, ipiv+=inc) {
            ip = (*ipiv) - descA.i - 1;
            if ( ip != j )
            {
                it = ip / descA.mb;
                i  = ip % descA.mb;
                lda2 = BLKLDD(descA, it);
                cblas_cswap(descA.n, A1       + j, lda1,
                                     A(it, 0) + i, lda2 );
            }
        }
    }

    return PLASMA_SUCCESS;
}
/***************************************************************************//**
 *
 **/
void QUARK_CORE_claswp_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, PLASMA_Complex32_t *Aij, 
                              int i1,  int i2, int *ipiv, int inc, PLASMA_Complex32_t *fakepanel)
{
    DAG_CORE_LASWP;
    QUARK_Insert_Task(
        quark, CORE_claswp_ontile_quark, task_flags,
        sizeof(PLASMA_desc),              &descA,     VALUE,
        sizeof(PLASMA_Complex32_t)*1,      Aij,           INOUT | LOCALITY,
        sizeof(int),                      &i1,        VALUE,
        sizeof(int),                      &i2,        VALUE,
        sizeof(int)*(i2-i1+1)*abs(inc),   ipiv,           INPUT,
        sizeof(int),                      &inc,       VALUE,
        sizeof(PLASMA_Complex32_t)*1,      fakepanel,     INOUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_claswp_ontile_quark = PCORE_claswp_ontile_quark
#define CORE_claswp_ontile_quark PCORE_claswp_ontile_quark
#endif
void CORE_claswp_ontile_quark(Quark *quark)
{
    int i1, i2, inc;
    int *ipiv;
    PLASMA_Complex32_t *A, *fake;
    PLASMA_desc descA;

    quark_unpack_args_7(quark, descA, A, i1, i2, ipiv, inc, fake);
    CORE_claswp_ontile(descA, i1, i2, ipiv, inc);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_claswp_ontile_f2(Quark *quark, Quark_Task_Flags *task_flags,
                                 PLASMA_desc descA, PLASMA_Complex32_t *Aij, 
                                 int i1,  int i2, int *ipiv, int inc, 
                                 PLASMA_Complex32_t *fake1, int szefake1, int flag1,
                                 PLASMA_Complex32_t *fake2, int szefake2, int flag2)
{
    DAG_CORE_LASWP;
    QUARK_Insert_Task(
        quark, CORE_claswp_ontile_f2_quark, task_flags,
        sizeof(PLASMA_desc),                &descA, VALUE,
        sizeof(PLASMA_Complex32_t)*1,        Aij,       INOUT | LOCALITY,
        sizeof(int),                        &i1,    VALUE,
        sizeof(int),                        &i2,    VALUE,
        sizeof(int)*(i2-i1+1)*abs(inc),      ipiv,      INPUT,
        sizeof(int),                        &inc,   VALUE,
        sizeof(PLASMA_Complex32_t)*szefake1, fake1, flag1,
        sizeof(PLASMA_Complex32_t)*szefake2, fake2, flag2,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_claswp_ontile_f2_quark = PCORE_claswp_ontile_f2_quark
#define CORE_claswp_ontile_f2_quark PCORE_claswp_ontile_f2_quark
#endif
void CORE_claswp_ontile_f2_quark(Quark *quark)
{
    int i1, i2, inc;
    int *ipiv;
    PLASMA_Complex32_t *A;
    PLASMA_desc descA;
    void *fake1, *fake2;

    quark_unpack_args_8(quark, descA, A, i1, i2, ipiv, inc, fake1, fake2);
    CORE_claswp_ontile(descA, i1, i2, ipiv, inc);
}

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex32_t
 *
 *  CORE_cswptr_ontile apply the claswp function on a matrix stored in
 *  tile layout, followed by a ctrsm on the first tile of the panel.
 *
 *******************************************************************************
 *
 *  @param[in,out] A
 *          The descriptor of the matrix A to permute.
 *
 *  @param[in] i1
 *          The first element of IPIV for which a row interchange will
 *          be done.
 *
 *  @param[in] i2
 *          The last element of IPIV for which a row interchange will
 *          be done.
 *
 *  @param[in] ipiv
 *          The pivot indices; Only the element in position i1 to i2
 *          are accessed. The pivot are offset by A.i.
 *
 *  @param[in] inc
 *          The increment between successive values of IPIV.  If IPIV
 *          is negative, the pivots are applied in reverse order.
 *
 *******************************************************************************
 */
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cswptr_ontile = PCORE_cswptr_ontile
#define CORE_cswptr_ontile PCORE_cswptr_ontile
#endif
int CORE_cswptr_ontile(PLASMA_desc descA, int i1, int i2, int *ipiv, int inc,
                       PLASMA_Complex32_t *Akk, int ldak)
{
    PLASMA_Complex32_t zone  = 1.0;
    int lda;
    int m = descA.mt == 1 ? descA.m : descA.mb;

    if ( descA.nt > 1 ) {
        coreblas_error(1, "Illegal value of descA.nt");
        return -1;
    }
    if ( i1 < 1 ) {
        coreblas_error(2, "Illegal value of i1");
        return -2;
    }
    if ( (i2 < i1) || (i2 > m) ) {
        coreblas_error(3, "Illegal value of i2");
        return -3;
    }

    CORE_claswp_ontile(descA, i1, i2, ipiv, inc);

    lda = BLKLDD(descA, 0);
    cblas_ctrsm( CblasColMajor, CblasLeft, CblasLower, 
                 CblasNoTrans, CblasUnit,
                 m, descA.n, CBLAS_SADDR(zone), 
                 Akk,     ldak, 
                 A(0, 0), lda );

    return PLASMA_SUCCESS;
}
/***************************************************************************//**
 *
 **/
void QUARK_CORE_cswptr_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, PLASMA_Complex32_t *Aij, 
                              int i1,  int i2, int *ipiv, int inc, 
                              PLASMA_Complex32_t *Akk, int ldak)
{
    DAG_CORE_TRSM;
    QUARK_Insert_Task(
        quark, CORE_cswptr_ontile_quark, task_flags,
        sizeof(PLASMA_desc),              &descA, VALUE,
        sizeof(PLASMA_Complex32_t)*1,      Aij,       INOUT | LOCALITY,
        sizeof(int),                      &i1,    VALUE,
        sizeof(int),                      &i2,    VALUE,
        sizeof(int)*(i2-i1+1)*abs(inc),    ipiv,      INPUT,
        sizeof(int),                      &inc,   VALUE,
        sizeof(PLASMA_Complex32_t)*ldak,   Akk,       INPUT,
        sizeof(int),                      &ldak,  VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cswptr_ontile_quark = PCORE_cswptr_ontile_quark
#define CORE_cswptr_ontile_quark PCORE_cswptr_ontile_quark
#endif
void CORE_cswptr_ontile_quark(Quark *quark)
{
    int i1, i2, inc, ldak;
    int *ipiv;
    PLASMA_Complex32_t *A, *Akk;
    PLASMA_desc descA;

    quark_unpack_args_8(quark, descA, A, i1, i2, ipiv, inc, Akk, ldak);
    CORE_cswptr_ontile(descA, i1, i2, ipiv, inc, Akk, ldak);
}

