/**
 *
 * @file pcpack.c
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
 * @generated c Thu Sep 15 12:09:23 2011
 *
 **/

#include <stdlib.h>
#include <sys/types.h>
#include "common.h"

/** ****************************************************************************
 *
 * @ingroup InPlaceTransformation
 *
 * plasma_pcpack pack all extra elements at the end of the matrix
 *
 *      +---------------+
 *      |               |
 *      |               |
 *      |     A11       |
 *      |               |
 *      |               |
 *      +---------------+
 *      |     A21       |
 *      +---------------+
 *
 * This matrix is initially stored as (example of Column Major, it's
 * the same for row major. We just consider the transpose matrix) :
 *  A11(:,0), A21(:,0), A11(:,1), A21(:,1), ... 
 *
 * On exit, it's stored as follow.
 *  A11(:,:), A12(:,:)
 *******************************************************************************
 *
 * @param[in] plasma
 *         Plasma context
 *
 * @param[in] m
 *         Number of rows in matrix A
 *
 * @param[in] n
 *         Number of columns in matrix A
 *
 * @param[in,out] A
 *         Matrix A to pack. (see above for entry and exit format)
 *
 * @param[in] m0
 *         Number of rows of A21
 *
 ******************************************************************************/
void plasma_pcpack(plasma_context_t *plasma)
{
    PLASMA_Complex32_t *A, *W, *Wl;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    int m, n, m0;
    int i, m1, size, rank, start, end, end2, bs, mod;

    plasma_unpack_args_6(m, n, A, m0, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;

    /* Quick return */
    if ( n <= 1 )
      return;

    m1 = m - m0;

    size = PLASMA_SIZE;
    rank = PLASMA_RANK;

    mod   = (n-1) % size;
    bs    = (n-1) / size;
    start = rank * bs;
    if ( rank < mod ) {
        bs++;
    }
    start += min( mod, rank );
    end    = start+bs;

    W  = (PLASMA_Complex32_t*)plasma_private_alloc(plasma, (m0*bs), PlasmaComplexFloat);
    Wl = (PLASMA_Complex32_t*)plasma_private_alloc(plasma, m1,      PlasmaComplexFloat);

    /* Save leftover pieces that are otherwise going to be overwritten */
    CORE_clacpy( PlasmaUpperLower, m0, bs, &(A[(int64_t)start*m+m1]), m, W, m0 );

    /* Pack A */
    end2 = ((n-1) / size) * size + 1;
    for(i=rank+1; i<end2; i+=size) {
        memcpy( Wl, &(A[i*m]), m1*sizeof(PLASMA_Complex32_t));
        plasma_barrier(plasma);
        memcpy( &(A[i*m1]), Wl, m1*sizeof(PLASMA_Complex32_t));
    }

    if ( rank < (n - end2)) {
        i = end2 + rank;
        memcpy( Wl, &(A[i*m]), m1*sizeof(PLASMA_Complex32_t));
        plasma_barrier(plasma);
        memcpy( &(A[i*m1]), Wl, m1*sizeof(PLASMA_Complex32_t));
    }
    else
        plasma_barrier(plasma);

    /* Restore leftover pieces */
    CORE_clacpy( PlasmaUpperLower, m0, bs, W, m0, &(A[(int64_t)m1*n+start*m0]), m0 );

    plasma_private_free(plasma, W);
    plasma_private_free(plasma, Wl);
}


/** ****************************************************************************
 *
 * @ingroup InPlaceTransformation
 *
 * plasma_pcunpack unpack all extra elements from the end of the matrix
 *
 *      +---------------+
 *      |               |
 *      |               |
 *      |     A11       |
 *      |               |
 *      |               |
 *      +---------------+
 *      |     A21       |
 *      +---------------+
 *
 * This matrix is initially stored as (example of Column Major, it's
 * the same for row major. We just consider the transpose matrix) :
 *  A11(:,:), A12(:,:)
 *
 * On exit, it's stored as follow.
 *  A11(:,0), A21(:,0), A11(:,1), A21(:,1), ... 
 *
 *******************************************************************************
 *
 * @param[in] plasma
 *         Plasma context
 *
 * @param[in] m
 *         Number of rows in matrix A
 *
 * @param[in] n
 *         Number of columns in matrix A
 *
 * @param[in,out] A
 *         Matrix A to pack. (see above for entry and exit format)
 *
 * @param[in] m0
 *         Number of rows of A21
 *
 ******************************************************************************/
void plasma_pcunpack(plasma_context_t *plasma)
{
    PLASMA_Complex32_t *A, *W, *Wl;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    int m, n, m0;
    int i, m1, size, rank, start, end, end2, bs, mod;

    plasma_unpack_args_6(m, n, A, m0, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;

    /* Quick return */
    if ( n <= 1 )
      return;

    m1 = m - m0;

    size = PLASMA_SIZE;
    rank = PLASMA_RANK;

    mod   = (n-1) % size;
    bs    = (n-1) / size;
    start = rank * bs;
    if ( rank < mod ) {
        bs++;
    }
    start += min( mod, rank );
    end    = start+bs;

    W  = (PLASMA_Complex32_t*)plasma_private_alloc(plasma, (m0*bs), PlasmaComplexFloat);
    Wl = (PLASMA_Complex32_t*)plasma_private_alloc(plasma, m1,      PlasmaComplexFloat);

    /* Save leftover pieces that are otherwise going to be overwritten */
    CORE_clacpy( PlasmaUpperLower, m0, bs, &(A[(int64_t)m1*n+start*m0]), m0, W, m0 );

    /* Unpack A */
    end2 = ((n-1) % size) ;
    for(i=n-1-rank; i>end2; i-=size) {
        memcpy( Wl, &(A[i*m1]), m1*sizeof(PLASMA_Complex32_t));
        plasma_barrier(plasma);
        memcpy( &(A[i*m]), Wl, m1*sizeof(PLASMA_Complex32_t));
    }

    if ( rank < end2 ) {
        i = rank+1;
        memcpy( Wl, &(A[i*m1]), m1*sizeof(PLASMA_Complex32_t));
        plasma_barrier(plasma);
        memcpy( &(A[i*m]), Wl, m1*sizeof(PLASMA_Complex32_t));
    }
    else
        plasma_barrier(plasma);

    /* Restore leftover pieces */
    CORE_clacpy( PlasmaUpperLower, m0, bs, W, m0, &(A[(int64_t)start*m+m1]), m );

    plasma_private_free(plasma, W);
    plasma_private_free(plasma, Wl);
}
