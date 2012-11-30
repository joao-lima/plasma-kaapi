/**
 *
 * @file dtrtri.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Thu Sep 15 12:09:16 2011
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup double
 *
 *  PLASMA_dtrtri - Computes the inverse of a complex upper or lower
 *  triangular matrix A.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] diag
 *          = PlasmaNonUnit: A is non-unit triangular;
 *          = PlasmaUnit:    A is unit triangular.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the triangular matrix A.  If UPLO = 'U', the
 *          leading N-by-N upper triangular part of the array A
 *          contains the upper triangular matrix, and the strictly
 *          lower triangular part of A is not referenced.  If UPLO =
 *          'L', the leading N-by-N lower triangular part of the array
 *          A contains the lower triangular matrix, and the strictly
 *          upper triangular part of A is not referenced.  If DIAG =
 *          'U', the diagonal elements of A are also not referenced and
 *          are assumed to be 1.  On exit, the (triangular) inverse of
 *          the original matrix, in the same storage format.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if i, A(i,i) is exactly zero.  The triangular
 *               matrix is singular and its inverse can not be computed.
 *
 *******************************************************************************
 *
 * @sa PLASMA_dtrtri_Tile
 * @sa PLASMA_dtrtri_Tile_Async
 * @sa PLASMA_ctrtri
 * @sa PLASMA_dtrtri
 * @sa PLASMA_strtri
 * @sa PLASMA_dpotri
 *
 ******************************************************************************/
int PLASMA_dtrtri(PLASMA_enum uplo, PLASMA_enum diag, int N,
                  double *A, int LDA)
{
    int NB;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dtrtri", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        plasma_error("PLASMA_dtrtri", "illegal value of uplo");
        return -1;
    }
    if (diag != PlasmaUnit && diag != PlasmaNonUnit) {
        plasma_error("PLASMA_dtrtri", "illegal value of diag");
        return -2;
    }
    if (N < 0) {
        plasma_error("PLASMA_dtrtri", "illegal value of N");
        return -3;
    }
    if (LDA < max(1, N)) {
        plasma_error("PLASMA_dtrtri", "illegal value of LDA");
        return -5;
    }
    /* Quick return */
    if (max(N, 0) == 0)
        return PLASMA_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_DPOSV, N, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dtrtri", "plasma_tune() failed");
        return status;
    }

    /* Set NT */
    NB   = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N, plasma_desc_mat_free(&(descA)) );
    } else {
        plasma_diplap2tile(  descA, A, NB, NB, LDA, N, 0, 0, N, N);
    }

    /* Call the tile interface */
    PLASMA_dtrtri_Tile_Async(uplo, diag, &descA, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dooptile2lap( descA, A, NB, NB, LDA, N );
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
    } else {
        plasma_diptile2lap( descA, A, NB, NB, LDA, N );
        plasma_dynamic_sync();
    }

    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup double_Tile
 *
 *  PLASMA_dtrtri_Tile - Computes the inverse of a complex upper or
 *  lower triangular matrix A.
 *  Tile equivalent of PLASMA_dtrtri().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] diag
 *          = PlasmaNonUnit: A is non-unit triangular;
 *          = PlasmaUnit:    A us unit triangular.
 *
 * @param[in] A
 *          On entry, the triangular matrix A.  If UPLO = 'U', the
 *          leading N-by-N upper triangular part of the array A
 *          contains the upper triangular matrix, and the strictly
 *          lower triangular part of A is not referenced.  If UPLO =
 *          'L', the leading N-by-N lower triangular part of the array
 *          A contains the lower triangular matrix, and the strictly
 *          upper triangular part of A is not referenced.  If DIAG =
 *          'U', the diagonal elements of A are also not referenced and
 *          are assumed to be 1.  On exit, the (triangular) inverse of
 *          the original matrix, in the same storage format.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval >0 if i, A(i,i) is exactly zero.  The triangular
 *               matrix is singular and its inverse can not be computed.
 *
 *******************************************************************************
 *
 * @sa PLASMA_dtrtri
 * @sa PLASMA_dtrtri_Tile_Async
 * @sa PLASMA_ctrtri_Tile
 * @sa PLASMA_dtrtri_Tile
 * @sa PLASMA_strtri_Tile
 * @sa PLASMA_dpotri_Tile
 *
 ******************************************************************************/
int PLASMA_dtrtri_Tile(PLASMA_enum uplo, PLASMA_enum diag, PLASMA_desc *A)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dtrtri_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_dtrtri_Tile_Async(uplo, diag, A, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup double_Tile_Async
 *
 *  PLASMA_dtrtri_Tile_Async - Computes the inverse of a complex upper or lower
 *  triangular matrix A.
 *  Non-blocking equivalent of PLASMA_dtrtri_Tile().
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
 * @sa PLASMA_dtrtri
 * @sa PLASMA_dtrtri_Tile
 * @sa PLASMA_ctrtri_Tile_Async
 * @sa PLASMA_dtrtri_Tile_Async
 * @sa PLASMA_strtri_Tile_Async
 * @sa PLASMA_dpotri_Tile_Async
 *
 ******************************************************************************/
int PLASMA_dtrtri_Tile_Async(PLASMA_enum uplo, PLASMA_enum diag, PLASMA_desc *A,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA = *A;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dtrtri_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_dtrtri_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_dtrtri_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dtrtri_Tile", "invalid descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_dtrtri_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        plasma_error("PLASMA_dtrtri_Tile", "illegal value of uplo");
        return plasma_request_fail(sequence, request, -1);
    }
    if (diag != PlasmaUnit && diag != PlasmaNonUnit) {
        plasma_error("PLASMA_dtrtri_Tile", "illegal value of diag");
        return plasma_request_fail(sequence, request, -2);
    }
    /* Quick return */
/*
    if (max(N, 0) == 0)
        return PLASMA_SUCCESS;
*/
    plasma_dynamic_call_5(plasma_pdtrtri,
        PLASMA_enum, uplo,
        PLASMA_enum, diag,
        PLASMA_desc, descA,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}