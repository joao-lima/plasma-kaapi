/**
 *
 * @file core_dshift.c
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
 * @generated d Thu Sep 15 12:09:00 2011
 *
 **/

#include <stdlib.h>
#include "common.h"
#include "quark.h"

/** ****************************************************************************
 *
 * @ingroup InPlaceTransformation
 * @ingroup CORE_double
 *
 *  CORE_dshiftw Shift a linear chain of block using a supplied workspace
 *      by following the cycle defined by:  k_(i+1) = (k_i * m) % q;
 *
 *******************************************************************************
 *
 * @param[in] s
 *         Start value in the cycle
 *
 * @param[in] cl
 *         Cycle length
 *         if cl == 0, all the permutations from the cycle are done
 *         else the cycle is split onto several threads and the number
 *         of permutation to do has to be specified to not get overlap
 *
 * @param[in] m
 *         Number of lines of tile A
 *
 * @param[in] n
 *         Number of columns of tile A
 *
 * @param[in] L
 *         Length of each block of data to move
 *
 * @param[in,out] A
 *         Matrix of size m-by-n with each element of size L.
 *         On exit, A = A', where A' contains the permutations
 *
 * @param[in] W
 *         Array of size L. On entry, must contain:
 *         W(:) = A(s*L:s*L+L-1)
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dshiftw = PCORE_dshiftw
#define CORE_dshiftw PCORE_dshiftw
#endif
void CORE_dshiftw(int s, int cl, int m, int n, int L, double *A, double *W) {
    int64_t k, k1;
    int     i, j, q, kL, k1L;

    q = m * n - 1;
    k = s;

    if( cl != 0 ) {
        for (i=1; i<cl; i++) {
            k1 = (k * m) % (int64_t)q;
            
            /* A(k*L:k*L+L-1) = A(k1*L:k1*L+L-1) */
            kL  = k *L;
            k1L = k1*L;

            for(j=0; j<L; j++) {
                A[kL+j] = A[k1L+j];
            }
            k = k1;
        }
    } 
    else {
        while (1) {
            k1 = (k * m) % (int64_t)q;
            if( k1 == s ) 
                break;
            
            /* A(k*L:k*L+L-1) = A(k1*L:k1*L+L-1) */
            kL  = k *L;
            k1L = k1*L;
            for (j=0; j<L; j++) {
                A[kL+j] = A[k1L+j];
            }
            k = k1;
        }
    }
    memcpy(&(A[k*L]), W, L*sizeof(double));
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dshiftw(Quark *quark, Quark_Task_Flags *task_flags,
                        int s, int cl, int m, int n, int L, double *A, double *W)
{
    DAG_CORE_SHIFTW;
    QUARK_Insert_Task(quark, CORE_dshiftw_quark, task_flags,
        sizeof(int),                      &s,   VALUE,
        sizeof(int),                      &cl,  VALUE,
        sizeof(int),                      &m,   VALUE,
        sizeof(int),                      &n,   VALUE,
        sizeof(int),                      &L,   VALUE,
        sizeof(double)*m*n*L, A,        INOUT,
        sizeof(double)*L,     W,        INPUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dshiftw_quark = PCORE_dshiftw_quark
#define CORE_dshiftw_quark PCORE_dshiftw_quark
#endif
void CORE_dshiftw_quark(Quark *quark)
{
    int s;
    int cl;
    int m;
    int n;
    int L;
    double *A;
    double *W;

    quark_unpack_args_7(quark, s, cl, m, n, L, A, W);
    CORE_dshiftw(s, cl, m, n, L, A, W);
}

/** ****************************************************************************
 *
 * @ingroup InPlaceTransformation
 *
 *  CORE_dshift Shift a cycle of block. Same as core_dshiftw but you
 *     don't need to provide the workspace.  As a matter of fact, the
 *     cycle cannot be split anymore to keep data coherency.
 *
 *******************************************************************************
 *
 * @param[in] s
 *         Start value in the cycle
 *
 * @param[in] m
 *         Number of lines of tile A
 *
 * @param[in] n
 *         Number of columns of tile A
 *
 * @param[in] L
 *         Length of each block of data to move
 *
 * @param[in,out] A
 *         Matrix of size m-by-n with each element of size L.
 *         On exit, A = A', where A' contains the permutations
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dshift = PCORE_dshift
#define CORE_dshift PCORE_dshift
#endif
void CORE_dshift(int s, int m, int n, int L, double *A) {
    double *W;

    W = (double*)malloc(L * sizeof(double));
    memcpy(W, &(A[s*L]), L*sizeof(double));
    CORE_dshiftw(s, 0, m, n, L, A, W);
    free(W);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dshift(Quark *quark, Quark_Task_Flags *task_flags,
                       int s, int m, int n, int L, double *A)
{
    DAG_CORE_SHIFT;
    QUARK_Insert_Task(quark, CORE_dshift_quark, task_flags,
        sizeof(int),                      &s,    VALUE,
        sizeof(int),                      &m,    VALUE,
        sizeof(int),                      &n,    VALUE,
        sizeof(int),                      &L,    VALUE,
        sizeof(double)*m*n*L, A,        INOUT | GATHERV,
        sizeof(double)*L,     NULL,     SCRATCH,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dshift_quark = PCORE_dshift_quark
#define CORE_dshift_quark PCORE_dshift_quark
#endif
void CORE_dshift_quark(Quark *quark)
{
    int s;
    int m;
    int n;
    int L;
    double *A;
    double *W;

    quark_unpack_args_6(quark, s, m, n, L, A, W);
    memcpy(W, &(A[s*L]), L*sizeof(double));
    CORE_dshiftw(s, 0, m, n, L, A, W);
}

