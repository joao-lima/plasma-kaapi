/**
 *
 * @file cpotrs.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Jakub Kurzak
 * @date 2010-11-15
 * @generated c Thu Sep 15 12:09:16 2011
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex32_t
 *
 *  PLASMA_cpotrs - Solves a system of linear equations A * X = B with a symmetric positive
 *  definite (or Hermitian positive definite in the complex case) matrix A using the Cholesky
 *  factorization A = U**H*U or A = L*L**H computed by PLASMA_cpotrf.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrix B. NRHS >= 0.
 *
 * @param[in] A
 *          The triangular factor U or L from the Cholesky factorization A = U**H*U or A = L*L**H,
 *          computed by PLASMA_cpotrf.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_cpotrs_Tile
 * @sa PLASMA_cpotrs_Tile_Async
 * @sa PLASMA_cpotrs
 * @sa PLASMA_dpotrs
 * @sa PLASMA_spotrs
 * @sa PLASMA_cpotrf
 *
 ******************************************************************************/
int PLASMA_cpotrs(PLASMA_enum uplo, int N, int NRHS,
                  PLASMA_Complex32_t *A, int LDA,
                  PLASMA_Complex32_t *B, int LDB)
{
    int NB;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA, descB;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_cpotrs", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        plasma_error("PLASMA_cpotrs", "illegal value of uplo");
        return -1;
    }
    if (N < 0) {
        plasma_error("PLASMA_cpotrs", "illegal value of N");
        return -2;
    }
    if (NRHS < 0) {
        plasma_error("PLASMA_cpotrs", "illegal value of NRHS");
        return -3;
    }
    if (LDA < max(1, N)) {
        plasma_error("PLASMA_cpotrs", "illegal value of LDA");
        return -5;
    }
    if (LDB < max(1, N)) {
        plasma_error("PLASMA_cpotrs", "illegal value of LDB");
        return -7;
    }
    /* Quick return */
    if (min(N, NRHS) == 0)
        return PLASMA_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_CPOSV, N, N, NRHS);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_cpotrs", "plasma_tune() failed");
        return status;
    }

    /* Set NT & NTRHS */
    NB    = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_cooplap2tile( descA, A, NB, NB, LDA, N,    0, 0, N, N   , plasma_desc_mat_free(&(descA)) );
        plasma_cooplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, N, NRHS, plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)));
    } else {
        plasma_ciplap2tile( descA, A, NB, NB, LDA, N,    0, 0, N, N   );
        plasma_ciplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, N, NRHS);
    }

    /* Call the tile interface */
    PLASMA_cpotrs_Tile_Async(uplo, &descA, &descB, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_cooptile2lap( descB, B, NB, NB, LDB, NRHS );
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
        plasma_desc_mat_free(&descB);
    } else {
        plasma_ciptile2lap( descA, A, NB, NB, LDA, N    );
        plasma_ciptile2lap( descB, B, NB, NB, LDB, NRHS );
        plasma_dynamic_sync();
    }
    
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex32_t_Tile
 *
 *  PLASMA_cpotrs_Tile - Solves a system of linear equations using previously
 *  computed Cholesky factorization.
 *  Tile equivalent of PLASMA_cpotrs().
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
 * @param[in] A
 *          The triangular factor U or L from the Cholesky factorization A = U**H*U or A = L*L**H,
 *          computed by PLASMA_cpotrf.
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_cpotrs
 * @sa PLASMA_cpotrs_Tile_Async
 * @sa PLASMA_cpotrs_Tile
 * @sa PLASMA_dpotrs_Tile
 * @sa PLASMA_spotrs_Tile
 * @sa PLASMA_cpotrf_Tile
 *
 ******************************************************************************/
int PLASMA_cpotrs_Tile(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_cpotrs_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_cpotrs_Tile_Async(uplo, A, B, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex32_t_Tile_Async
 *
 *  PLASMA_cpotrs_Tile_Async - Solves a system of linear equations using previously
 *  computed Cholesky factorization.
 *  Non-blocking equivalent of PLASMA_cpotrs_Tile().
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
 * @sa PLASMA_cpotrs
 * @sa PLASMA_cpotrs_Tile
 * @sa PLASMA_cpotrs_Tile_Async
 * @sa PLASMA_dpotrs_Tile_Async
 * @sa PLASMA_spotrs_Tile_Async
 * @sa PLASMA_cpotrf_Tile_Async
 *
 ******************************************************************************/
int PLASMA_cpotrs_Tile_Async(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA = *A;
    PLASMA_desc descB = *B;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_cpotrs_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_cpotrs_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_cpotrs_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_cpotrs_Tile", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (plasma_desc_check(&descB) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_cpotrs_Tile", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (descA.nb != descA.mb || descB.nb != descB.mb) {
        plasma_error("PLASMA_cpotrs_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        plasma_error("PLASMA_cpotrs_Tile", "illegal value of uplo");
        return plasma_request_fail(sequence, request, -1);
    }
    /* Quick return */
/*
    if (min(N, NRHS) == 0)
        return PLASMA_SUCCESS;
*/
    plasma_parallel_call_9(plasma_pctrsm,
        PLASMA_enum, PlasmaLeft,
        PLASMA_enum, uplo,
        PLASMA_enum, uplo == PlasmaUpper ? PlasmaConjTrans : PlasmaNoTrans,
        PLASMA_enum, PlasmaNonUnit,
        PLASMA_Complex32_t, 1.0,
        PLASMA_desc, descA,
        PLASMA_desc, descB,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    plasma_parallel_call_9(plasma_pctrsm,
        PLASMA_enum, PlasmaLeft,
        PLASMA_enum, uplo,
        PLASMA_enum, uplo == PlasmaUpper ? PlasmaNoTrans : PlasmaConjTrans,
        PLASMA_enum, PlasmaNonUnit,
        PLASMA_Complex32_t, 1.0,
        PLASMA_desc, descA,
        PLASMA_desc, descB,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}
