/**
 *
 * @file dsyrk.c
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
 *  PLASMA_dsyrk - Performs one of the hermitian rank k operations
 *
 *    \f[ C = \alpha [ op( A ) \times g( op( A )' )] + \beta C \f],
 *
 *  where op( X ) is one of
 *
 *    op( X ) = X  or op( X ) = g( X' )
 *
 *  where alpha and beta are real scalars, C is an n-by-n hermitian
 *  matrix and A is an n-by-k matrix in the first case and a k-by-n
 *  matrix in the second case.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = PlasmaUpper: Upper triangle of C is stored;
 *          = PlasmaLower: Lower triangle of C is stored.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed or ugate transposed:
 *          = PlasmaNoTrans:   A is not transposed;
 *          = PlasmaTrans  :   A is transposed.
 *
 * @param[in] N
 *          N specifies the order of the matrix C. N must be at least zero.
 *
 * @param[in] K
 *          K specifies the number of columns of the matrix op( A ).
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is K when trans = PlasmaNoTrans,
 *          and is N otherwise.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA must be at least
 *          max( 1, N ), otherwise LDA must be at least max( 1, K ).
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array uplo part of the matrix is overwritten
 *          by the uplo part of the updated matrix.
 *
 * @param[in] LDC
 *          The leading dimension of the array C. LDC >= max( 1, N ).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_dsyrk_Tile
 * @sa PLASMA_csyrk
 * @sa PLASMA_dsyrk
 * @sa PLASMA_ssyrk
 *
 ******************************************************************************/
int PLASMA_dsyrk(PLASMA_enum uplo, PLASMA_enum trans, int N, int K,
                 double alpha, double *A, int LDA,
                 double beta,  double *C, int LDC)
{
    int NB;
    int Am, An;
    int status;
    PLASMA_desc descA, descC;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dsyrk", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if ((uplo != PlasmaUpper) && (uplo != PlasmaLower)) {
        plasma_error("PLASMA_dsyrk", "illegal value of uplo");
        return -1;
    }
    if ((trans != PlasmaNoTrans) && (trans != PlasmaTrans)) {
        plasma_error("PLASMA_dsyrk", "illegal value of trans");
        return -2;
    }
    if ( trans == PlasmaNoTrans ) { 
        Am = N; An = K;
    } else {
        Am = K; An = N;
    }
    if (N < 0) {
        plasma_error("PLASMA_dsyrk", "illegal value of N");
        return -3;
    }
    if (K < 0) {
        plasma_error("PLASMA_dsyrk", "illegal value of K");
        return -4;
    }
    if (LDA < max(1, Am)) {
        plasma_error("PLASMA_dsyrk", "illegal value of LDA");
        return -7;
    }
    if (LDC < max(1, N)) {
        plasma_error("PLASMA_dsyrk", "illegal value of LDC");
        return -10;
    }

    /* Quick return */
    if (N == 0 ||
        ((alpha == (double)0.0 || K == 0.0) && beta == (double)1.0))
        return PLASMA_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNBSIZE */
    status = plasma_tune(PLASMA_FUNC_DSYRK, N, K, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dsyrk", "plasma_tune() failed");
        return status;
    }

    /* Set MT & NT & KT */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dooplap2tile( descA, A, NB, NB, LDA, An, 0, 0, Am, An, plasma_desc_mat_free(&(descA)) );
        plasma_dooplap2tile( descC, C, NB, NB, LDC, N,  0, 0, N,  N,  plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descC)));
    } else {
        plasma_diplap2tile( descA, A, NB, NB, LDA, An, 0, 0, Am, An );
        plasma_diplap2tile( descC, C, NB, NB, LDC, N,  0, 0, N,  N );
    }

    /* Call the tile interface */
    PLASMA_dsyrk_Tile_Async(uplo, trans, alpha, &descA, beta, &descC, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dooptile2lap( descC, C, NB, NB, LDC, N );
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
        plasma_desc_mat_free(&descC);
    } else {
        plasma_diptile2lap( descA, A, NB, NB, LDA, An );
        plasma_diptile2lap( descC, C, NB, NB, LDC, N );
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
 *  PLASMA_dsyrk_Tile - Performs rank k update.
 *  Tile equivalent of PLASMA_dsyrk().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = PlasmaUpper: Upper triangle of C is stored;
 *          = PlasmaLower: Lower triangle of C is stored.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed or ugate transposed:
 *          = PlasmaNoTrans:   A is not transposed;
 *          = PlasmaTrans:     A is transposed.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is K when trans = PlasmaNoTrans,
 *          and is N otherwise.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array uplo part of the matrix is overwritten
 *          by the uplo part of the updated matrix.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_dsyrk_Tile
 * @sa PLASMA_csyrk
 * @sa PLASMA_dsyrk
 * @sa PLASMA_ssyrk
 *
 ******************************************************************************/
int PLASMA_dsyrk_Tile(PLASMA_enum uplo, PLASMA_enum trans,
                      double alpha, PLASMA_desc *A,
                      double beta,  PLASMA_desc *C)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dsyrk_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_dsyrk_Tile_Async(uplo, trans, alpha, A, beta, C, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup double_Tile_Async
 *
 *  PLASMA_dsyrk_Tile_Async - Performs rank-k update.
 *  Non-blocking equivalent of PLASMA_dsyrk_Tile().
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
 * @sa PLASMA_dsyrk
 * @sa PLASMA_dsyrk_Tile
 * @sa PLASMA_csyrk_Tile_Async
 * @sa PLASMA_dsyrk_Tile_Async
 * @sa PLASMA_ssyrk_Tile_Async
 *
 ******************************************************************************/
int PLASMA_dsyrk_Tile_Async(PLASMA_enum uplo, PLASMA_enum trans,
                            double alpha, PLASMA_desc *A,
                            double beta,  PLASMA_desc *C,
                            PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    PLASMA_desc descA = *A;
    PLASMA_desc descC = *C;
    int N, K;
    int Am, An, Ai, Aj, Amb, Anb;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dsyrk_Tile_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_dsyrk_Tile_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_dsyrk_Tile_Async", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dsyrk_Tile_Async", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (plasma_desc_check(&descC) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dsyrk_Tile_Async", "invalid third descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ((uplo != PlasmaUpper) && (uplo != PlasmaLower)) {
        plasma_error("PLASMA_dsyrk", "illegal value of uplo");
        return plasma_request_fail(sequence, request, -1);
    }
    if ((trans != PlasmaNoTrans) && (trans != PlasmaTrans)) {
        plasma_error("PLASMA_dsyrk", "illegal value of transA");
        return plasma_request_fail(sequence, request, -2);
    }

    if ( trans == PlasmaNoTrans ) {
        Am  = descA.m;
        An  = descA.n;
        Amb = descA.mb;
        Anb = descA.nb;
        Ai  = descA.i;
        Aj  = descA.j;
    } else {
        Am  = descA.n;
        An  = descA.m;
        Amb = descA.nb;
        Anb = descA.mb;
        Ai  = descA.j;
        Aj  = descA.i;
    }

    if (descC.mb != descC.nb) {
        plasma_error("PLASMA_dsyrk_Tile_Async", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (Amb != descC.mb) {
        plasma_error("PLASMA_dsyrk_Tile_Async", "tile sizes have to match");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (descC.m != descC.n) {
        plasma_error("PLASMA_dsyrk_Tile_Async", "only square matrix C is supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (Am != descC.m) {
        plasma_error("PLASMA_dsyrk_Tile_Async", "sizes of matrices have to match");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    N = descC.m;
    K = An;

    /* Quick return */
    if ( N == 0 ||
        ((alpha == (double)0.0 || K == 0) && beta == (double)1.0))
        return PLASMA_SUCCESS;

    plasma_parallel_call_8(plasma_pdsyrk,
        PLASMA_enum, uplo,
        PLASMA_enum, trans,
        double, alpha,
        PLASMA_desc, descA,
        double, beta,
        PLASMA_desc, descC,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}
