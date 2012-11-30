/**
 *
 * @file core_sbrdalg.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Azzam Haidar
 * @date 2011-05-15
 * @generated s Thu Sep 15 12:09:02 2011
 *
 **/
#include <lapacke.h>
#include "common.h"


/***************************************************************************//**
 *
 * @ingroup CORE_float
 *
 *  CORE_sbrdalg is a part of the bidiagonal reduction algorithm (bulgechasing).
 *  It correspond to a local driver of the kernels that should be executed on a 
 *  single core.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *         @arg PlasmaLower:
 *         @arg PlasmaUpper: 
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in] NB
 *          The size of the Bandwidth of the matrix A,
 *          which correspond to the tile size. NB >= 0.
 *
 * @param[in] pA
 *          A pointer to the descriptor of the matrix A.
 *
 * @param[out] V
 *          float array, dimension (N).
 *          The scalar elementary reflectors are written in this 
 *          array. So it is used as a workspace for V at each step 
 *          of the bulge chasing algorithm.
 *
 * @param[out] TAU
 *          float array, dimension (N).
 *          The scalar factors of the elementary reflectors are written  
 *          in thisarray. So it is used as a workspace for TAU at each step 
 *          of the bulge chasing algorithm.
 *
 * @param[in] i
 *          Integer that refer to the current sweep. (outer loop).
 *
 * @param[in] j
 *          Integer that refer to the sweep to chase.(inner loop).
 *
 * @param[in] m
 *          Integer that refer to a sweep step, to ensure order dependencies. 
 *
 * @param[in] grsiz
 *          Integer that refer to the size of a group.
 *          group mean the number of kernel that should be executed sequentially
 *          on the same core. 
 *          group size is a trade-off between locality (cache reuse) and parallelism.
 *          a small group size increase parallelism while a large group size increase 
 *          cache reuse.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sbrdalg = PCORE_sbrdalg
#define CORE_sbrdalg PCORE_sbrdalg
#endif
void CORE_sbrdalg(PLASMA_enum uplo, int N, int NB, 
                  PLASMA_desc *pA, float *V, float *TAU,
                  int i, int j, int m, int grsiz)
{
    int    k, shift=3;
    int    myid, colpt, stind, edind, blklastind, stepercol;
    size_t eltsize;
    PLASMA_desc A = *pA;

    eltsize = plasma_element_size(A.dtyp);

    k = shift / grsiz;
    stepercol = (k*grsiz == shift) ? k : k+1;
    for (k = 0; k < grsiz; k++){
        myid = (i-j)*(stepercol*grsiz) +(m-1)*grsiz + k+1;
        if(myid%2 ==0) {
            colpt      = (myid/2) * NB + 1 + j - 1;
            stind      = colpt - NB + 1;
            edind      = min(colpt, N);
            blklastind = colpt;
        } else {
            colpt      = ((myid+1)/2)*NB + 1 +j -1 ;
            stind      = colpt-NB+1;
            edind      = min(colpt,N);
            if( (stind>=edind-1) && (edind==N) )
                blklastind = N;
            else
                blklastind = 0;
        }

        if( myid == 1 )
           CORE_sgbelr(uplo, N, &A, V, TAU, stind, edind, eltsize);
        else if(myid%2 == 0)
           CORE_sgbrce(uplo, N, &A, V, TAU, stind, edind, eltsize);
        else /*if(myid%2 == 1)*/
           CORE_sgblrx(uplo, N, &A, V, TAU, stind, edind, eltsize);

        if(blklastind >= (N-1))  break;
    }
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sbrdalg(Quark *quark, Quark_Task_Flags *task_flags,
                        int uplo,
                        int N, int NB,
                        PLASMA_desc *A,
                        float *V,
                        float *TAU,
                        int i, int j, int m, int grsiz, int BAND,
                        int *PCOL, int *ACOL, int *MCOL)
{
    QUARK_Insert_Task(quark, CORE_sbrdalg_quark,   task_flags,
        sizeof(int),               &uplo,               VALUE,
        sizeof(int),                  &N,               VALUE,
        sizeof(int),                 &NB,               VALUE,
        sizeof(PLASMA_desc),           A,               NODEP,
        sizeof(float),    V,               NODEP,
        sizeof(float),    TAU,               NODEP,
        sizeof(int),                  &i,               VALUE,
        sizeof(int),                  &j,               VALUE,
        sizeof(int),                  &m,               VALUE,
        sizeof(int),              &grsiz,               VALUE,
        sizeof(int),                PCOL,               INPUT,
        sizeof(int),                ACOL,               INPUT,
        sizeof(int),                MCOL,              OUTPUT | LOCALITY,
        0);

}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sbrdalg_quark = PCORE_sbrdalg_quark
#define CORE_sbrdalg_quark PCORE_sbrdalg_quark
#endif
void CORE_sbrdalg_quark(Quark *quark)
{
    PLASMA_desc *pA;
    float *V;
    float *TAU;
    int    uplo;
    int    N, NB;
    int    i, j, m, grsiz;

    quark_unpack_args_10(quark, uplo, N, NB, pA, V, TAU, i, j, m, grsiz);
    CORE_sbrdalg(uplo, N, NB, pA, V, TAU, i, j, m, grsiz);
}