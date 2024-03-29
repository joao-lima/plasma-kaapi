/**
 *
 * @file splrnt.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Thu Sep 15 12:09:18 2011
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup float
 *
 *  PLASMA_splrnt - Generate a random matrix by tiles.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of A.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[out] A
 *          On exit, The random matrix A generated.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in] seed
 *          The seed used in the random generation.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_splrnt_Tile
 * @sa PLASMA_splrnt_Tile_Async
 * @sa PLASMA_cplrnt
 * @sa PLASMA_dplrnt
 * @sa PLASMA_splrnt
 * @sa PLASMA_splgsy
 * @sa PLASMA_splgsy
 *
 ******************************************************************************/
int PLASMA_splrnt( int M, int N,
                   float *A, int LDA,
                   unsigned long long int seed )
{
    int NB;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_splrnt", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        plasma_error("PLASMA_splrnt", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        plasma_error("PLASMA_splrnt", "illegal value of N");
        return -2;
    }
    if (LDA < max(1, M)) {
        plasma_error("PLASMA_splrnt", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (min(M, N) == 0)
        return PLASMA_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_SGEMM, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_splrnt", "plasma_tune() failed");
        return status;
    }
    
    /* Set NT */
    NB = PLASMA_NB;
    plasma_sequence_create(plasma, &sequence);
    descA = plasma_desc_init(
        PlasmaRealFloat, NB, NB, NB*NB,
        LDA, N, 0, 0, M, N);
    descA.mat = A;

    /* Call the tile interface */
    PLASMA_splrnt_Tile_Async( &descA, seed, sequence, &request );

    plasma_siptile2lap( descA, A, NB, NB, LDA, N );
    plasma_dynamic_sync();

    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);

    return status;
}

/***************************************************************************//**
 *
 * @ingroup float_Tile
 *
 *  PLASMA_splrnt_Tile - Generate a random matrix by tiles.
 *  Tile equivalent of PLASMA_splrnt().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          On exit, The random matrix A generated.
 *
 * @param[in] seed
 *          The seed used in the random generation.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_splrnt
 * @sa PLASMA_splrnt_Tile_Async
 * @sa PLASMA_cplrnt_Tile
 * @sa PLASMA_dplrnt_Tile
 * @sa PLASMA_splrnt_Tile
 * @sa PLASMA_splgsy_Tile
 * @sa PLASMA_splgsy_Tile
 *
 ******************************************************************************/
int PLASMA_splrnt_Tile( PLASMA_desc *A,
                        unsigned long long int seed )
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_splrnt_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_splrnt_Tile_Async( A, seed, sequence, &request );
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup float_Tile_Async
 *
 *  PLASMA_splrnt_Tile_Async - Generate a random matrix by tiles.
 *  Non-blocking equivalent of PLASMA_splrnt_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations ar runtime.
 *
 *******************************************************************************
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa PLASMA_splrnt
 * @sa PLASMA_splrnt_Tile
 * @sa PLASMA_cplrnt_Tile_Async
 * @sa PLASMA_dplrnt_Tile_Async
 * @sa PLASMA_splrnt_Tile_Async
 * @sa PLASMA_splgsy_Tile_Async
 * @sa PLASMA_splgsy_Tile_Async
 *
 ******************************************************************************/
int PLASMA_splrnt_Tile_Async( PLASMA_desc     *A,
                              unsigned long long int seed,
                              PLASMA_sequence *sequence, 
                              PLASMA_request  *request)
{
    PLASMA_desc descA = *A;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_splrnt_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_splrnt_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_splrnt_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_splrnt_Tile", "invalid descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_splrnt_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if (min( descA.m, descA.n ) == 0)
        return PLASMA_SUCCESS;

    plasma_parallel_call_4(plasma_psplrnt,
        PLASMA_desc,      descA,
        unsigned long long int, seed,
        PLASMA_sequence*, sequence,
        PLASMA_request*,  request);

    return PLASMA_SUCCESS;
}
