/**
 *
 * @file core_cgetrip.c
 *
 *  PLASMA InPlaceTransformation module
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 *  This work is the implementation of an inplace transformation 
 *  based on the GKK algorithm by Gustavson, Karlsson, Kagstrom 
 *  and its fortran implementation.
 *
 * @version 2.4.2
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 * @generated c Thu Sep 15 12:09:00 2011
 *
 **/

#include <stdlib.h>
#include "common.h"
#include "quark.h"

/** ****************************************************************************
 *
 * @ingroup InPlaceTransformation
 * @ingroup CORE_PLASMA_Complex32_t
 *
 *  CORE_cgetrip transposes a m-by-n matrix in place using an extra
 *      workspace of size m-by-n.  
 *      Note : For square tile, workspace is not used.
 *
 *******************************************************************************
 *
 * @param[in] m
 *         Number of lines of tile A
 *
 * @param[in] n
 *         Number of columns of tile A
 *
 * @param[in,out] A
 *         Tile of size m-by-n
 *         On exit, A = trans(A)
 *
 * @param[out] W
 *         Workspace of size n-by-m if n != m, NULL otherwise.
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cgetrip = PCORE_cgetrip
#define CORE_cgetrip PCORE_cgetrip
#endif
void CORE_cgetrip(int m, int n, PLASMA_Complex32_t *A, PLASMA_Complex32_t *W) {
    PLASMA_Complex32_t t;
    int i, j;
    
    if( m != n ) {
        /* rectangular transposition (use workspace) */
        for (i=0; i<m; i++) {
            for (j=0; j<n; j++) {
                W[j+i*n] = A[i+j*m];
            }
        }
        memcpy(A, W, m*n*sizeof(PLASMA_Complex32_t));
    }
    else {
        /* square transposition (swap pairwise) */
        for (i=0; i<m; i++) {
            for (j=i+1; j<n; j++) {
                t        = A[j+i*n];
                A[j+i*n] = A[i+j*m];
                A[i+j*m] = t;
            }
        }
    }
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_cgetrip(Quark *quark, Quark_Task_Flags *task_flags,
                        int m, int n, PLASMA_Complex32_t *A, int szeA)
{
    DAG_CORE_GETRIP;
    QUARK_Insert_Task(quark, CORE_cgetrip_quark, task_flags,
        sizeof(int),                     &m,   VALUE,
        sizeof(int),                     &n,   VALUE,
        sizeof(PLASMA_Complex32_t)*szeA, A,        INOUT,
        sizeof(PLASMA_Complex32_t)*szeA, NULL,     SCRATCH,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cgetrip_quark = PCORE_cgetrip_quark
#define CORE_cgetrip_quark PCORE_cgetrip_quark
#endif
void CORE_cgetrip_quark(Quark *quark)
{
    int m;
    int n;
    PLASMA_Complex32_t *A;
    PLASMA_Complex32_t *W;

    quark_unpack_args_4(quark, m, n, A, W);
    CORE_cgetrip(m, n, A, W);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_cgetrip_f1(Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n, 
                           PLASMA_Complex32_t *A,    int szeA,
                           PLASMA_Complex32_t *fake, int szeF, int paramF)
{
    DAG_CORE_GETRIP;
    QUARK_Insert_Task(
        quark, CORE_cgetrip_f1_quark, task_flags,
        sizeof(int),                     &m,   VALUE,
        sizeof(int),                     &n,   VALUE,
        sizeof(PLASMA_Complex32_t)*szeA, A,        INOUT,
        sizeof(PLASMA_Complex32_t)*szeA, NULL,     SCRATCH,
        sizeof(PLASMA_Complex32_t)*szeF, fake,     paramF,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cgetrip_f1_quark = PCORE_cgetrip_f1_quark
#define CORE_cgetrip_f1_quark PCORE_cgetrip_f1_quark
#endif
void CORE_cgetrip_f1_quark(Quark *quark)
{
    int m;
    int n;
    PLASMA_Complex32_t *A;
    PLASMA_Complex32_t *W;
    PLASMA_Complex32_t *fake;

    quark_unpack_args_5(quark, m, n, A, W, fake);
    CORE_cgetrip(m, n, A, W);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_cgetrip_f2(Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n, 
                           PLASMA_Complex32_t *A,    int szeA,
                           PLASMA_Complex32_t *fake1, int szeF1, int paramF1,
                           PLASMA_Complex32_t *fake2, int szeF2, int paramF2)
{
    DAG_CORE_GETRIP;
    QUARK_Insert_Task(
        quark, CORE_cgetrip_f2_quark, task_flags,
        sizeof(int),                     &m,   VALUE,
        sizeof(int),                     &n,   VALUE,
        sizeof(PLASMA_Complex32_t)*szeA, A,        INOUT,
        sizeof(PLASMA_Complex32_t)*szeA, NULL,     SCRATCH,
        sizeof(PLASMA_Complex32_t)*szeF1, fake1,     paramF1,
        sizeof(PLASMA_Complex32_t)*szeF2, fake2,     paramF2,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cgetrip_f2_quark = PCORE_cgetrip_f2_quark
#define CORE_cgetrip_f2_quark PCORE_cgetrip_f2_quark
#endif
void CORE_cgetrip_f2_quark(Quark *quark)
{
    int m;
    int n;
    PLASMA_Complex32_t *A;
    PLASMA_Complex32_t *W;
    PLASMA_Complex32_t *fake1;
    PLASMA_Complex32_t *fake2;

    quark_unpack_args_6(quark, m, n, A, W, fake1, fake2);
    CORE_cgetrip(m, n, A, W);
}
