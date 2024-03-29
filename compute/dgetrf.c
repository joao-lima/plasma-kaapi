/**
 *
 * @file dgetrf.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @date 2009-11-15
 *
 * @generated d Thu Sep 15 12:09:15 2011
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup double
 *
 *  PLASMA_dgetrf - Computes an LU factorization of a general M-by-N matrix A
 *  using the tile LU algorithm with partial tile pivoting with row interchanges.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the tile factors L and U from the factorization.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] IPIV
 *          The pivot indices that define the permutations.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if i, U(i,i) is exactly zero. The factorization has been completed,
 *               but the factor U is exactly singular, and division by zero will occur
 *               if it is used to solve a system of equations.
 *
 *******************************************************************************
 *
 * @sa PLASMA_dgetrf_Tile
 * @sa PLASMA_dgetrf_Tile_Async
 * @sa PLASMA_cgetrf
 * @sa PLASMA_dgetrf
 * @sa PLASMA_sgetrf
 *
 ******************************************************************************/
int PLASMA_dgetrf(int M, int N,
                  double *A, int LDA,
                  int *IPIV)
{
    int NB, NBNB, minMN;
    int status;
    PLASMA_desc descA ;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dgetrf", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        plasma_error("PLASMA_dgetrf", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        plasma_error("PLASMA_dgetrf", "illegal value of N");
        return -2;
    }
    if (LDA < max(1, M)) {
        plasma_error("PLASMA_dgetrf", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (min(M, N) == 0)
        return PLASMA_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNBSIZE */
    status = plasma_tune(PLASMA_FUNC_DGESV, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dgetrf", "plasma_tune() failed");
        return status;
    }

    /* Set NT & NTRHS */
    NB   = PLASMA_NB;
    NBNB = NB*NB;

    plasma_sequence_create(plasma, &sequence);

    descA = plasma_desc_init(
        PlasmaRealDouble,
        NB, NB, NBNB,
        LDA, N, 0, 0, M, N);
    descA.mat = A;

    minMN = min(M, N);
    memset(IPIV, 0, minMN*sizeof(int));

    /* Call the tile interface */
    plasma_dynamic_call_4(plasma_pdgetrf_reclap,
        PLASMA_desc, descA,
        int*, IPIV,
        PLASMA_sequence*, sequence,
        PLASMA_request*, &request);

    plasma_dynamic_sync();

    /*
     * Generate the correct IPIV (Has to be move in a task)
     */
    { 
        int i, inc, tmp, j;
        for(i=1; i<descA.mt; i++) {
            inc = i*descA.mb;
            tmp = min( minMN - inc, descA.mb);
            if ( tmp < 1 )
              break;
            
            for (j=0; j<tmp; j++)
                IPIV[inc+j] = IPIV[inc+j] + inc;
        }
    }

    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);

    return status;
}

/***************************************************************************//**
 *
 * @ingroup double_Tile
 *
 *  PLASMA_dgetrf_Tile - Computes the tile LU factorization of a matrix.
 *  Tile equivalent of PLASMA_dgetrf().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the tile factors L and U from the factorization.
 *
 * @param[out] L
 *          On exit, auxiliary factorization data, related to the tile L factor,
 *          required by PLASMA_dgetrs to solve the system of equations.
 *
 * @param[out] IPIV
 *          The pivot indices that define the permutations (not equivalent to LAPACK).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval >0 if i, U(i,i) is exactly zero. The factorization has been completed,
 *               but the factor U is exactly singular, and division by zero will occur
 *               if it is used to solve a system of equations.
 *
 *******************************************************************************
 *
 * @sa PLASMA_dgetrf
 * @sa PLASMA_dgetrf_Tile_Async
 * @sa PLASMA_cgetrf_Tile
 * @sa PLASMA_dgetrf_Tile
 * @sa PLASMA_sgetrf_Tile
 * @sa PLASMA_dgetrs_Tile
 *
 ******************************************************************************/
int PLASMA_dgetrf_Tile(PLASMA_desc *A, int *IPIV)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dgetrf_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_dgetrf_Tile_Async(A, IPIV, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup double_Tile_Async
 *
 *  PLASMA_dgetrf_Tile_Async - Computes the tile LU factorization of a matrix.
 *  Non-blocking equivalent of PLASMA_dgetrf_Tile().
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
 * @sa PLASMA_dgetrf
 * @sa PLASMA_dgetrf_Tile
 * @sa PLASMA_cgetrf_Tile_Async
 * @sa PLASMA_dgetrf_Tile_Async
 * @sa PLASMA_sgetrf_Tile_Async
 * @sa PLASMA_dgetrs_Tile_Async
 *
 ******************************************************************************/
int PLASMA_dgetrf_Tile_Async(PLASMA_desc *A, int *IPIV,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA = *A;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dgetrf_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_dgetrf_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_dgetrf_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dgetrf_Tile", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_dgetrf_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
/*
    if (min(M, N) == 0)
        return PLASMA_SUCCESS;
*/

    plasma_dynamic_call_3(
        plasma_pdbarrier_tl2pnl,
        PLASMA_desc, descA,
        PLASMA_sequence*, sequence,
        PLASMA_request*,  request);
    
    plasma_dynamic_call_4(plasma_pdgetrf_rectil,
        PLASMA_desc, descA,
        int*, IPIV,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    plasma_dynamic_call_3(
        plasma_pdbarrier_pnl2tl,
        PLASMA_desc, descA,
        PLASMA_sequence*, sequence,
        PLASMA_request*,  request);
    
    return PLASMA_SUCCESS;
}
