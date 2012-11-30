/**
 *
 * @file sormqr.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @date 2010-11-15
 * @generated s Thu Sep 15 12:09:17 2011
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup float
 *
 *  PLASMA_sormqr - overwrites the general M-by-N matrix C with Q*C, where Q is an orthogonal
 *  matrix (unitary in the complex case) defined as the product of elementary reflectors returned
 *  by PLASMA_sgeqrf. Q is of order M.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Intended usage:
 *          = PlasmaLeft:  apply Q or Q**T from the left;
 *          = PlasmaRight: apply Q or Q**T from the right.
 *          Currently only PlasmaLeft is supported.
 *
 * @param[in] trans
 *          Intended usage:
 *          = PlasmaNoTrans:   no transpose, apply Q;
 *          = PlasmaTrans: ugate transpose, apply Q**T.
 *          Currently only PlasmaTrans is supported.
 *
 * @param[in] M
 *          The number of rows of the matrix C. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix C. N >= 0.
 *
 * @param[in] K
 *          The number of columns of elementary tile reflectors whose product defines the matrix Q.
 *          If side == PlasmaLeft,  M >= K >= 0.
 *          If side == PlasmaRight, N >= K >= 0.
 *
 * @param[in] A
 *          Details of the QR factorization of the original matrix A as returned by PLASMA_sgeqrf.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. 
 *          If side == PlasmaLeft,  LDA >= max(1,M).
 *          If side == PlasmaRight, LDA >= max(1,N).
 *
 * @param[in] T
 *          Auxiliary factorization data, computed by PLASMA_sgeqrf.
 *
 * @param[in,out] B
 *          On entry, the M-by-N matrix B.
 *          On exit, B is overwritten by Q*B or Q**T*B.
 *
 * @param[in] LDB
 *          The leading dimension of the array C. LDC >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_sormqr_Tile
 * @sa PLASMA_sormqr_Tile_Async
 * @sa PLASMA_cormqr
 * @sa PLASMA_dormqr
 * @sa PLASMA_sormqr
 * @sa PLASMA_sgeqrf
 *
 ******************************************************************************/
int PLASMA_sormqr(PLASMA_enum side, PLASMA_enum trans, int M, int N, int K,
                  float *A, int LDA,
                  float *T,
                  float *B, int LDB)
{
    int NB, IB, IBNB, Am, MT, KT;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA, descB, descT;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_sormqr", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    if ( side == PlasmaLeft ) {
        Am = M;
    } else {
        Am = N;
    }

    /* Check input arguments */
    if ((side != PlasmaLeft) && (side != PlasmaRight)) {
        plasma_error("PLASMA_sormqr", "illegal value of side");
        return -1;
    }
    if ((trans != PlasmaTrans) && (trans != PlasmaNoTrans)){
        plasma_error("PLASMA_sormqr", "illegal value of trans");
        return -2;
    }
    if (M < 0) {
        plasma_error("PLASMA_sormqr", "illegal value of M");
        return -3;
    }
    if (N < 0) {
        plasma_error("PLASMA_sormqr", "illegal value of N");
        return -4;
    }
    if ( (K < 0) || (K > Am) ) {
        plasma_error("PLASMA_sormqr", "illegal value of K");
        return -5;
    }
    if ( LDA < max(1, Am) ) {
        plasma_error("PLASMA_sormqr", "illegal value of LDA");
        return -7;
    }
    if (LDB < max(1, M)) {
        plasma_error("PLASMA_sormqr", "illegal value of LDB");
        return -10;
    }
    /* Quick return - currently NOT equivalent to LAPACK's:
     * CALL DLASET( 'Full', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB ) */
    if (min(M, min(N, K)) == 0)
        return PLASMA_SUCCESS;

    /* Tune NB & IB depending on M, K & N; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_SGELS, M, K, N);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sormqr", "plasma_tune() failed");
        return status;
    }

    /* Set MT, NT & NTRHS */
    NB   = PLASMA_NB;
    IB   = PLASMA_IB;
    IBNB = IB*NB;
    MT   = (Am%NB==0) ? (Am/NB) : (Am/NB+1);
    KT   = (K%NB==0)  ? (K /NB) : (K /NB+1);

    plasma_sequence_create(plasma, &sequence);

    if (plasma->householder == PLASMA_FLAT_HOUSEHOLDER) {
        descT = plasma_desc_init(
            PlasmaRealFloat,
            IB, NB, IBNB,
            MT*IB, KT*NB, 0, 0, MT*IB, KT*NB);
    }
    else {
        /* Double the size of T to accomodate the tree reduction phase */
        descT = plasma_desc_init(
            PlasmaRealFloat,
            IB, NB, IBNB,
            MT*IB, 2*KT*NB, 0, 0, MT*IB, 2*KT*NB);
    }
    descT.mat = T;

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_sooplap2tile( descA, A, NB, NB, LDA, K, 0, 0, Am, K, plasma_desc_mat_free(&(descA)) );
        plasma_sooplap2tile( descB, B, NB, NB, LDB, N, 0, 0, M,  N, plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)));
    } else {
        plasma_siplap2tile( descA, A, NB, NB, LDA, K, 0, 0, Am, K);
        plasma_siplap2tile( descB, B, NB, NB, LDB, N, 0, 0, M,  N);
    }

    /* Call the tile interface */
    PLASMA_sormqr_Tile_Async(
        side, trans, &descA, &descT, &descB, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_sooptile2lap( descB, B, NB, NB, LDB, N );
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
        plasma_desc_mat_free(&descB);
    } else {
        plasma_siptile2lap( descA, A, NB, NB, LDA, K );
        plasma_siptile2lap( descB, B, NB, NB, LDB, N );
        plasma_dynamic_sync();
    }

    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup float_Tile
 *
 *  PLASMA_sormqr_Tile - overwrites the general M-by-N matrix C with Q*C, where Q is an orthogonal
 *  matrix (unitary in the complex case) defined as the product of elementary reflectors returned
 *  by PLASMA_sgeqrf_Tile Q is of order M.
 *  All matrices are passed through descriptors. All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Intended usage:
 *          = PlasmaLeft:  apply Q or Q**T from the left;
 *          = PlasmaRight: apply Q or Q**T from the right.
 *          Currently only PlasmaLeft is supported.
 *
 * @param[in] trans
 *          Intended usage:
 *          = PlasmaNoTrans:   no transpose, apply Q;
 *          = PlasmaTrans: ugate transpose, apply Q**T.
 *          Currently only PlasmaTrans is supported.
 *
 * @param[in] A
 *          Details of the QR factorization of the original matrix A as returned by PLASMA_sgeqrf.
 *
 * @param[in] T
 *          Auxiliary factorization data, computed by PLASMA_sgeqrf.
 *
 * @param[in,out] B
 *          On entry, the M-by-N matrix B.
 *          On exit, B is overwritten by Q*B or Q**T*B.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_sormqr
 * @sa PLASMA_sormqr_Tile_Async
 * @sa PLASMA_cormqr_Tile
 * @sa PLASMA_dormqr_Tile
 * @sa PLASMA_sormqr_Tile
 * @sa PLASMA_sgeqrf_Tile
 *
 ******************************************************************************/
int PLASMA_sormqr_Tile(PLASMA_enum side, PLASMA_enum trans,
                       PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_sormqr_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_sormqr_Tile_Async(side, trans, A, T, B, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup float_Tile_Async
 *
 *  Non-blocking equivalent of PLASMA_sormqr_Tile().
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
 * @sa PLASMA_sormqr
 * @sa PLASMA_sormqr_Tile
 * @sa PLASMA_cormqr_Tile_Async
 * @sa PLASMA_dormqr_Tile_Async
 * @sa PLASMA_sormqr_Tile_Async
 * @sa PLASMA_sgeqrf_Tile_Async
 *
 ******************************************************************************/
int PLASMA_sormqr_Tile_Async(PLASMA_enum side, PLASMA_enum trans,
                             PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA = *A;
    PLASMA_desc descT = *T;
    PLASMA_desc descB = *B;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_sormqr_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_sormqr_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_sormqr_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sormqr_Tile", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (plasma_desc_check(&descT) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sormqr_Tile", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (plasma_desc_check(&descB) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sormqr_Tile", "invalid third descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (descA.nb != descA.mb || descB.nb != descB.mb) {
        plasma_error("PLASMA_sormqr_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ((side != PlasmaLeft) && (side != PlasmaRight)) {
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ((trans != PlasmaTrans) && (trans != PlasmaNoTrans)){
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Quick return - currently NOT equivalent to LAPACK's:
     * CALL DLASET( 'Full', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB ) */
/*
    if (min(M, min(N, K)) == 0)
        return PLASMA_SUCCESS;
*/
    if (plasma->householder == PLASMA_FLAT_HOUSEHOLDER) {
        if ( (trans == PlasmaTrans) &&
             (side == PlasmaLeft) ) {
            plasma_parallel_call_7(plasma_psormqr,
                PLASMA_enum, side,
                PLASMA_enum, trans,
                PLASMA_desc, descA,
                PLASMA_desc, descB,
                PLASMA_desc, descT,
                PLASMA_sequence*, sequence,
                PLASMA_request*, request);
        }
        else {
            plasma_dynamic_call_7(plasma_psormqr,
                PLASMA_enum, side,
                PLASMA_enum, trans,
                PLASMA_desc, descA,
                PLASMA_desc, descB,
                PLASMA_desc, descT,
                PLASMA_sequence*, sequence,
                PLASMA_request*, request);
        }
    }
    else {
        plasma_dynamic_call_8(plasma_psormqrrh,
            PLASMA_enum, side,
            PLASMA_enum, trans,
            PLASMA_desc, descA,
            PLASMA_desc, descB,
            PLASMA_desc, descT,
            PLASMA_enum, PLASMA_RHBLK,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);
    }

    return PLASMA_SUCCESS;
}
