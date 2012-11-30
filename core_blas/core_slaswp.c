/**
 *
 * @file core_slaswp.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Thu Sep 15 12:09:01 2011
 *
 **/
#include <lapacke.h>
#include "common.h"

#define A(m, n) BLKADDR(descA, float, m, n)

/***************************************************************************//**
 *
 * @ingroup CORE_float
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slaswp = PCORE_slaswp
#define CORE_slaswp PCORE_slaswp
#endif
void CORE_slaswp(int N, float *A, int LDA, int I1,  int I2, int *IPIV, int INC)
{
    LAPACKE_slaswp_work( LAPACK_COL_MAJOR, N, A, LDA, I1, I2, IPIV, INC );
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_slaswp(Quark *quark, Quark_Task_Flags *task_flags,
                       int n, float *A, int lda, 
                       int i1,  int i2, int *ipiv, int inc)
{
    DAG_CORE_LASWP;
    QUARK_Insert_Task(
        quark, CORE_slaswp_quark, task_flags,
        sizeof(int),                      &n,    VALUE,
        sizeof(float)*lda*n,  A,        INOUT | LOCALITY,
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
#pragma weak CORE_slaswp_quark = PCORE_slaswp_quark
#define CORE_slaswp_quark PCORE_slaswp_quark
#endif
void CORE_slaswp_quark(Quark *quark)
{
    int n, lda, i1, i2, inc;
    int *ipiv;
    float *A;
    
    quark_unpack_args_7(quark, n, A, lda, i1, i2, ipiv, inc);
    LAPACKE_slaswp_work(LAPACK_COL_MAJOR, n, A, lda, i1, i2, ipiv, inc );
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_slaswp_f2(Quark *quark, Quark_Task_Flags *task_flags,
                          int n, float *A, int lda, 
                          int i1,  int i2, int *ipiv, int inc,
                          float *fake1, int szefake1, int flag1,
                          float *fake2, int szefake2, int flag2)
{
    DAG_CORE_LASWP;
    QUARK_Insert_Task(
        quark, CORE_slaswp_f2_quark, task_flags,
        sizeof(int),                        &n,     VALUE,
        sizeof(float)*lda*n,    A,         INOUT | LOCALITY,
        sizeof(int),                        &lda,   VALUE,
        sizeof(int),                        &i1,    VALUE,
        sizeof(int),                        &i2,    VALUE,
        sizeof(int)*n,                       ipiv,      INPUT,
        sizeof(int),                        &inc,   VALUE,
        sizeof(float)*szefake1, fake1,     flag1,
        sizeof(float)*szefake2, fake2,     flag2,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slaswp_f2_quark = PCORE_slaswp_f2_quark
#define CORE_slaswp_f2_quark PCORE_slaswp_f2_quark
#endif
void CORE_slaswp_f2_quark(Quark* quark)
{
    int n, lda, i1, i2, inc;
    int *ipiv;
    float *A;
    void *fake1, *fake2;
    
    quark_unpack_args_9(quark, n, A, lda, i1, i2, ipiv, inc, fake1, fake2);
    LAPACKE_slaswp_work(LAPACK_COL_MAJOR, n, A, lda, i1, i2, ipiv, inc );
}

/***************************************************************************//**
 *
 * @ingroup CORE_float
 *
 *  CORE_slaswp_ontile apply the slaswp function on a matrix stored in
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
#pragma weak CORE_slaswp_ontile = PCORE_slaswp_ontile
#define CORE_slaswp_ontile PCORE_slaswp_ontile
#endif
int CORE_slaswp_ontile(PLASMA_desc descA, int i1, int i2, int *ipiv, int inc)
{
    int i, j, ip, it;
    float *A1;
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
                cblas_sswap(descA.n, A1       + j, lda1,
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
                cblas_sswap(descA.n, A1       + j, lda1,
                                     A(it, 0) + i, lda2 );
            }
        }
    }

    return PLASMA_SUCCESS;
}
/***************************************************************************//**
 *
 **/
void QUARK_CORE_slaswp_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, float *Aij, 
                              int i1,  int i2, int *ipiv, int inc, float *fakepanel)
{
    DAG_CORE_LASWP;
    QUARK_Insert_Task(
        quark, CORE_slaswp_ontile_quark, task_flags,
        sizeof(PLASMA_desc),              &descA,     VALUE,
        sizeof(float)*1,      Aij,           INOUT | LOCALITY,
        sizeof(int),                      &i1,        VALUE,
        sizeof(int),                      &i2,        VALUE,
        sizeof(int)*(i2-i1+1)*abs(inc),   ipiv,           INPUT,
        sizeof(int),                      &inc,       VALUE,
        sizeof(float)*1,      fakepanel,     INOUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slaswp_ontile_quark = PCORE_slaswp_ontile_quark
#define CORE_slaswp_ontile_quark PCORE_slaswp_ontile_quark
#endif
void CORE_slaswp_ontile_quark(Quark *quark)
{
    int i1, i2, inc;
    int *ipiv;
    float *A, *fake;
    PLASMA_desc descA;

    quark_unpack_args_7(quark, descA, A, i1, i2, ipiv, inc, fake);
    CORE_slaswp_ontile(descA, i1, i2, ipiv, inc);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_slaswp_ontile_f2(Quark *quark, Quark_Task_Flags *task_flags,
                                 PLASMA_desc descA, float *Aij, 
                                 int i1,  int i2, int *ipiv, int inc, 
                                 float *fake1, int szefake1, int flag1,
                                 float *fake2, int szefake2, int flag2)
{
    DAG_CORE_LASWP;
    QUARK_Insert_Task(
        quark, CORE_slaswp_ontile_f2_quark, task_flags,
        sizeof(PLASMA_desc),                &descA, VALUE,
        sizeof(float)*1,        Aij,       INOUT | LOCALITY,
        sizeof(int),                        &i1,    VALUE,
        sizeof(int),                        &i2,    VALUE,
        sizeof(int)*(i2-i1+1)*abs(inc),      ipiv,      INPUT,
        sizeof(int),                        &inc,   VALUE,
        sizeof(float)*szefake1, fake1, flag1,
        sizeof(float)*szefake2, fake2, flag2,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slaswp_ontile_f2_quark = PCORE_slaswp_ontile_f2_quark
#define CORE_slaswp_ontile_f2_quark PCORE_slaswp_ontile_f2_quark
#endif
void CORE_slaswp_ontile_f2_quark(Quark *quark)
{
    int i1, i2, inc;
    int *ipiv;
    float *A;
    PLASMA_desc descA;
    void *fake1, *fake2;

    quark_unpack_args_8(quark, descA, A, i1, i2, ipiv, inc, fake1, fake2);
    CORE_slaswp_ontile(descA, i1, i2, ipiv, inc);
}

/***************************************************************************//**
 *
 * @ingroup CORE_float
 *
 *  CORE_sswptr_ontile apply the slaswp function on a matrix stored in
 *  tile layout, followed by a strsm on the first tile of the panel.
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
#pragma weak CORE_sswptr_ontile = PCORE_sswptr_ontile
#define CORE_sswptr_ontile PCORE_sswptr_ontile
#endif
int CORE_sswptr_ontile(PLASMA_desc descA, int i1, int i2, int *ipiv, int inc,
                       float *Akk, int ldak)
{
    float zone  = 1.0;
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

    CORE_slaswp_ontile(descA, i1, i2, ipiv, inc);

    lda = BLKLDD(descA, 0);
    cblas_strsm( CblasColMajor, CblasLeft, CblasLower, 
                 CblasNoTrans, CblasUnit,
                 m, descA.n, (zone), 
                 Akk,     ldak, 
                 A(0, 0), lda );

    return PLASMA_SUCCESS;
}
/***************************************************************************//**
 *
 **/
void QUARK_CORE_sswptr_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, float *Aij, 
                              int i1,  int i2, int *ipiv, int inc, 
                              float *Akk, int ldak)
{
    DAG_CORE_TRSM;
    QUARK_Insert_Task(
        quark, CORE_sswptr_ontile_quark, task_flags,
        sizeof(PLASMA_desc),              &descA, VALUE,
        sizeof(float)*1,      Aij,       INOUT | LOCALITY,
        sizeof(int),                      &i1,    VALUE,
        sizeof(int),                      &i2,    VALUE,
        sizeof(int)*(i2-i1+1)*abs(inc),    ipiv,      INPUT,
        sizeof(int),                      &inc,   VALUE,
        sizeof(float)*ldak,   Akk,       INPUT,
        sizeof(int),                      &ldak,  VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sswptr_ontile_quark = PCORE_sswptr_ontile_quark
#define CORE_sswptr_ontile_quark PCORE_sswptr_ontile_quark
#endif
void CORE_sswptr_ontile_quark(Quark *quark)
{
    int i1, i2, inc, ldak;
    int *ipiv;
    float *A, *Akk;
    PLASMA_desc descA;

    quark_unpack_args_8(quark, descA, A, i1, i2, ipiv, inc, Akk, ldak);
    CORE_sswptr_ontile(descA, i1, i2, ipiv, inc, Akk, ldak);
}

