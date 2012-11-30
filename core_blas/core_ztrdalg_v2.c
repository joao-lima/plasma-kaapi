/**
 *
 * @file core_ztrdalg_v2.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Azzam Haidar
 * @date 2011-05-15
 * @precisions normal z -> c d s
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 *  CORE_ztrdalg_v2 is a part of the tridiagonal reduction algorithm (bulgechasing)
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
 *          PLASMA_Complex64_t array, dimension (N).
 *          The scalar elementary reflectors are written in this
 *          array. So it is used as a workspace for V at each step
 *          of the bulge chasing algorithm.
 *
 * @param[out] TAU
 *          PLASMA_Complex64_t array, dimension (N).
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
#pragma weak CORE_ztrdalg_v2 = PCORE_ztrdalg_v2
#define CORE_ztrdalg_v2 PCORE_ztrdalg_v2
#endif
void CORE_ztrdalg_v2(PLASMA_enum        uplo,
                     PLASMA_desc        *pA,
             PLASMA_Complex64_t *V,
             PLASMA_Complex64_t *TAU,
                     int grsiz, int lcsweep, int id, int blksweep)
{
  PLASMA_desc A = *pA;
  size_t eltsize = plasma_element_size(A.dtyp);
  int   N, NB;
  int   i, blkid, st, ed, KDM1;
  int   NT=pA->nt;


  N      = A.m;
  NB     = A.mb;
  KDM1   = NB-1;


  /* code for all tiles */
  for (i = 0; i < grsiz ; i++) {
     blkid = id+i;
     st    = min(blkid*NB+lcsweep+1, N-1);
     ed    = min(st+KDM1, N-1);
     /*printf("  COUCOU voici N %5d   NB %5d      st %5d     ed %5d lcsweep %5d     id %5d   blkid %5d\n",N, NB, st, ed, lcsweep, id, blkid);*/

     if(st==ed) /* quick return in case of last tile */
        return;

     st =st +1; /* because kernel are still in fortran way */
     ed =ed +1;
     if(blkid==blksweep){
        CORE_zhbelr(uplo, N, &A, V, TAU, st, ed, eltsize);
        if(id!=(NT-1))CORE_zhbrce(uplo, N, &A, V, TAU, st, ed, eltsize);
     }else{
        CORE_zhblrx(uplo, N, &A, V, TAU, st, ed, eltsize);
        if(id!=(NT-1))CORE_zhbrce(uplo, N, &A, V, TAU, st, ed, eltsize);
     }
  }

}

/***************************************************************************//**
 *
 **/
#define A(m,n) BLKADDR(A, PLASMA_Complex64_t, m, n)
void QUARK_CORE_ztrdalg_v2(Quark *quark, Quark_Task_Flags *task_flags,
                        int uplo,
                        PLASMA_desc *pA,
                        PLASMA_Complex64_t *V,
                        PLASMA_Complex64_t *TAU,
                        int grsiz, int lcsweep, int id, int blksweep)
{
  Quark_Task *MYTASK;
  PLASMA_desc A=*pA;
  int ii, cur_id,  NT=pA->nt;

           //printf("coucou from quark function id %d    lcsweep %d    blksweep %d   grsiz %d  NT %d\n", id, lcsweep, blksweep, grsiz, NT);
           MYTASK = QUARK_Task_Init( quark, CORE_ztrdalg_v2_quark,   task_flags);
           QUARK_Task_Pack_Arg(quark, MYTASK,  sizeof(int),                  &uplo,               VALUE );
           QUARK_Task_Pack_Arg(quark, MYTASK,  sizeof(PLASMA_desc),             pA,               NODEP );
           QUARK_Task_Pack_Arg(quark, MYTASK,  sizeof(PLASMA_Complex64_t),       V,               NODEP );
           QUARK_Task_Pack_Arg(quark, MYTASK,  sizeof(PLASMA_Complex64_t),     TAU,               NODEP );
           QUARK_Task_Pack_Arg(quark, MYTASK,  sizeof(int),                 &grsiz,               VALUE );
           QUARK_Task_Pack_Arg(quark, MYTASK,  sizeof(int),               &lcsweep,               VALUE );
           QUARK_Task_Pack_Arg(quark, MYTASK,  sizeof(int),                    &id,               VALUE );
           QUARK_Task_Pack_Arg(quark, MYTASK,  sizeof(int),               &blksweep,               VALUE );

           QUARK_Task_Pack_Arg(quark, MYTASK, sizeof(PLASMA_Complex64_t),    A(id,   id  ),           INOUT );
           QUARK_Task_Pack_Arg(quark, MYTASK, sizeof(PLASMA_Complex64_t),    A(id+1, id  ),           INOUT );
           if( id<(NT-1) )
              QUARK_Task_Pack_Arg(quark, MYTASK, sizeof(PLASMA_Complex64_t),    A(id+1, id+1),           INOUT );
           if( id<(NT-2) )
              QUARK_Task_Pack_Arg(quark, MYTASK, sizeof(PLASMA_Complex64_t),    A(id+2, id+1),           INOUT );

       cur_id = id;
           for (ii = 1; ii < grsiz ; ii++) {
                cur_id = cur_id+1;
        if( id<(NT-1) )
                   QUARK_Task_Pack_Arg(quark, MYTASK, sizeof(PLASMA_Complex64_t),    A(cur_id+1, cur_id+1),           INOUT );
        if( id<(NT-2) )
                    QUARK_Task_Pack_Arg(quark, MYTASK, sizeof(PLASMA_Complex64_t),    A(cur_id+2, cur_id+1),           INOUT );
           }

           QUARK_Insert_Task_Packed(quark, MYTASK);
}
#undef A
/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ztrdalg_v2_quark = PCORE_ztrdalg_v2_quark
#define CORE_ztrdalg_v2_quark PCORE_ztrdalg_v2_quark
#endif
void CORE_ztrdalg_v2_quark(Quark *quark)
{
    PLASMA_desc *pA;
    PLASMA_Complex64_t *V;
    PLASMA_Complex64_t *TAU;
    int    uplo;
    int    grsiz, lcsweep, id, blksweep;

    quark_unpack_args_8(quark, uplo, pA, V, TAU, grsiz, lcsweep, id, blksweep);
    CORE_ztrdalg_v2(uplo, pA, V, TAU, grsiz, lcsweep, id, blksweep);
}
