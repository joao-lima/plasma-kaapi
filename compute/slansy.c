/**
 *
 * @file slansy.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Thu Sep 15 12:09:17 2011
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup float
 *
 *  PLASMA_slansy returns the value
 *
 *     slansy = ( max(abs(A(i,j))), NORM = PlasmaMaxNorm
 *              (
 *              ( norm1(A),         NORM = PlasmaOneNorm
 *              (
 *              ( normI(A),         NORM = PlasmaInfNorm
 *              (
 *              ( normF(A),         NORM = PlasmaFrobeniusNorm
 *
 *  where norm1 denotes the one norm of a matrix (maximum column sum),
 *  normI denotes the infinity norm of a matrix (maximum row sum) and
 *  normF denotes the Frobenius norm of a matrix (square root of sum
 *  of squares). Note that max(abs(A(i,j))) is not a consistent matrix
 *  norm.
 *
 *******************************************************************************
 *
 * @param[in] norm
 *          = PlasmaMaxNorm: Max norm
 *          = PlasmaOneNorm: One norm
 *          = PlasmaInfNorm: Infinity norm
 *          = PlasmaFrobeniusNorm: Frobenius norm
 *
 * @param[in] uplo
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] N
 *          The number of columns/rows of the matrix A. N >= 0. When N = 0,
 *          the returned value is set to zero.
 *
 * @param[in] A
 *          The N-by-N matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in] work
 *          real array of dimension PLASMA_SIZE is
 *          PLASMA_STATIC_SCHEDULING is used, and NULL otherwise.
 *
 *******************************************************************************
 *
 * @return
 *          \retval the norm described above.
 *
 *******************************************************************************
 *
 * @sa PLASMA_slansy_Tile
 * @sa PLASMA_slansy_Tile_Async
 * @sa PLASMA_clansy
 * @sa PLASMA_dlansy
 * @sa PLASMA_slansy
 *
 ******************************************************************************/
float PLASMA_slansy(PLASMA_enum norm, PLASMA_enum uplo, int N,
                     float *A, int LDA, float *work)
{
    int NB;
    int status;
    float value;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_slansy", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ( (norm != PlasmaMaxNorm) && (norm != PlasmaOneNorm)
        && (norm != PlasmaInfNorm) && (norm != PlasmaFrobeniusNorm) ) {
        plasma_error("PLASMA_slansy", "illegal value of norm");
        return -1;
    }
    if ( (uplo != PlasmaUpper) && (uplo != PlasmaLower) ) {
        plasma_error("PLASMA_slansy", "illegal value of uplo");
        return -2;
    }
    if (N < 0) {
        plasma_error("PLASMA_slansy", "illegal value of N");
        return -3;
    }
    if (LDA < max(1, N)) {
        plasma_error("PLASMA_slansy", "illegal value of LDA");
        return -5;
    }

    /* Quick return */
    if ( N == 0)
      return (float)0.0;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_SGEMM, N, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_slansy", "plasma_tune() failed");
        return status;
    }

    /* Set NT */
    NB   = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_sooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N, plasma_desc_mat_free(&(descA)) );
    } else {
        plasma_siplap2tile(  descA, A, NB, NB, LDA, N, 0, 0, N, N);
    }

    /* Call the tile interface */
    PLASMA_slansy_Tile_Async(norm, uplo, &descA, work, &value, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
    } else {
        plasma_siptile2lap( descA, A, NB, NB, LDA, N );
        plasma_dynamic_sync();
    }

    plasma_sequence_destroy(plasma, sequence);
    return value;
}

/***************************************************************************//**
 *
 * @ingroup float_Tile
 *
 *  PLASMA_slansy_Tile - Tile equivalent of PLASMA_slansy().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] norm
 *          = PlasmaMaxNorm: Max norm
 *          = PlasmaOneNorm: One norm
 *          = PlasmaInfNorm: Infinity norm
 *          = PlasmaFrobeniusNorm: Frobenius norm
 *
 * @param[in] uplo
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] A
 *          On entry, the triangular factor U or L.
 *          On exit, if UPLO = 'U', the upper triangle of A is
 *          overwritten with the upper triangle of the product U * U';
 *          if UPLO = 'L', the lower triangle of A is overwritten with
 *          the lower triangle of the product L' * L.
 *
 * @param[in] work
 *          real array of dimension PLASMA_SIZE is
 *          PLASMA_STATIC_SCHEDULING is used, and NULL otherwise.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_slansy
 * @sa PLASMA_slansy_Tile_Async
 * @sa PLASMA_clansy_Tile
 * @sa PLASMA_dlansy_Tile
 * @sa PLASMA_slansy_Tile
 *
 ******************************************************************************/
float PLASMA_slansy_Tile(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc *A, float *work)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    float value;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_slansy_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_slansy_Tile_Async(norm, uplo, A, work, &value, sequence, &request);
    plasma_dynamic_sync();
    plasma_sequence_destroy(plasma, sequence);
    return value;
}

/***************************************************************************//**
 *
 * @ingroup float_Tile_Async
 *
 *  PLASMA_slansy_Tile_Async - Non-blocking equivalent of PLASMA_slansy_Tile().
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
 * @sa PLASMA_slansy
 * @sa PLASMA_slansy_Tile
 * @sa PLASMA_clansy_Tile_Async
 * @sa PLASMA_dlansy_Tile_Async
 * @sa PLASMA_slansy_Tile_Async
 *
 ******************************************************************************/
int PLASMA_slansy_Tile_Async(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc *A, float *work, float *value,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA = *A;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_slansy_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_slansy_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_slansy_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_slansy_Tile", "invalid descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_slansy_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ( (norm != PlasmaMaxNorm) && (norm != PlasmaOneNorm)
         && (norm != PlasmaInfNorm) && (norm != PlasmaFrobeniusNorm) ) {
        plasma_error("PLASMA_slansy_Tile", "illegal value of norm");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ( (uplo != PlasmaUpper) && (uplo != PlasmaLower) ) {
        plasma_error("PLASMA_slansy_Tile", "illegal value of uplo");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
    if ( descA.m == 0) {
        *value = 0.0;
        return PLASMA_SUCCESS;
    }

    plasma_parallel_call_7(plasma_pslansy,
        PLASMA_enum, norm,
        PLASMA_enum, uplo,
        PLASMA_desc, descA,
        float*, work,
        float*, value,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}
