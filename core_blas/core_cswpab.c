/**
 *
 * @file core_cswpab.c
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
 *  CORE_cswpab swaps two adjacent contiguous blocks of data.
 *
 *      n1                     n2
 * +-------------+-------------------------------+
 *
 * become :
 *              n2                        n1
 * +-------------------------------+-------------+
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *         Array of size i+n1+n2.
 *         On entry, a block of size n1 followed by a block of size n2.
 *         On exit, the block of size n1 follows the block of size n2.
 *
 * @param[in] i
 *         First block starts at A[i].
 *
 * @param[in] n1
 *         Size of the first block to swap.
 *
 * @param[in] n2
 *         Size of the second block to swap.
 *
 * @param[out] work
 *         Workspace array of size min(n1, n2).
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cswpab = PCORE_cswpab
#define CORE_cswpab PCORE_cswpab
#endif
void CORE_cswpab(int i, int n1, int n2, PLASMA_Complex32_t *A, PLASMA_Complex32_t *work) {
    PLASMA_Complex32_t *A0 = &(A[i]);
    PLASMA_Complex32_t *A1 = &(A[i+n1]);
    PLASMA_Complex32_t *A2 = &(A[i+n2]);
    int j;
    
    if( n1 < n2 ) {
        memcpy(work,  A0, n1*sizeof(PLASMA_Complex32_t));
        for (j=0; j<n2; j++)
            A0[j] = A1[j];
        memcpy(A2, work,  n1*sizeof(PLASMA_Complex32_t));
    } else {
        memcpy(work,  A1, n2*sizeof(PLASMA_Complex32_t));
        for (j=n1-1; j>-1; j--)
            A2[j] = A0[j];
        memcpy(A0, work,  n2*sizeof(PLASMA_Complex32_t));
    }
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_cswpab(Quark *quark, Quark_Task_Flags *task_flags,
                       int i, int n1, int n2,
                       PLASMA_Complex32_t *A, int szeA)
{
    DAG_CORE_SWPAB;
    QUARK_Insert_Task(
        quark, CORE_cswpab_quark, task_flags,
        sizeof(int),                           &i,   VALUE,
        sizeof(int),                           &n1,  VALUE,
        sizeof(int),                           &n2,  VALUE,
        sizeof(PLASMA_Complex32_t)*szeA,       A,            INOUT,
        sizeof(PLASMA_Complex32_t)*min(n1,n2), NULL,         SCRATCH,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cswpab_quark = PCORE_cswpab_quark
#define CORE_cswpab_quark PCORE_cswpab_quark
#endif
void CORE_cswpab_quark(Quark *quark)
{
    int i;
    int n1;
    int n2;
    PLASMA_Complex32_t *A;
    PLASMA_Complex32_t *work;
    
    quark_unpack_args_5(quark, i, n1, n2, A, work);
    CORE_cswpab( i, n1, n2, A, work);
}
