/**
 *
 * @file dormlq.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Dulceneia Becker
 * @date 2010-11-15
 * @generated d Thu Sep 15 12:09:17 2011
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup double
 *
 *  PLASMA_dormlq - overwrites the general M-by-N matrix C with Q*C, where Q is an orthogonal
 *  matrix (unitary in the complex case) defined as the product of elementary reflectors returned
 *  by PLASMA_dgelqf. Q is of order M.
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
 *          The number of rows of elementary tile reflectors whose product defines the matrix Q.
 *          M >= K >= 0.
 *
 * @param[in] A
 *          Details of the LQ factorization of the original matrix A as returned by PLASMA_dgelqf.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,K).
 *
 * @param[in] T
 *          Auxiliary factorization data, computed by PLASMA_dgelqf.
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
 * @sa PLASMA_dormlq_Tile
 * @sa PLASMA_dormlq_Tile_Async
 * @sa PLASMA_cunmlq
 * @sa PLASMA_dormlq
 * @sa PLASMA_sormlq
 * @sa PLASMA_dgelqf
 *
 ******************************************************************************/
int PLASMA_dormlq(PLASMA_enum side, PLASMA_enum trans, int M, int N, int K,
                  double *A, int LDA,
                  double *T,
                  double *B, int LDB)
{
    int NB, IB, IBNB, KT, NT, An;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA, descB, descT;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dormlq", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    if (side == PlasmaLeft)
        An = M;
    else 
        An = N;

    /* Check input arguments */
    if ( (side != PlasmaLeft) && (side != PlasmaRight) ) {
        plasma_error("PLASMA_dormlq", "illegal value of side");
        return -1;
    }
    if ( (trans != PlasmaTrans) && (trans != PlasmaNoTrans) ){
        plasma_error("PLASMA_dormlq", "illegal value of trans");
        return -2;
    }
    if (M < 0) {
        plasma_error("PLASMA_dormlq", "illegal value of M");
        return -3;
    }
    if (N < 0) {
        plasma_error("PLASMA_dormlq", "illegal value of N");
        return -4;
    }
    if ((K < 0) || (K > An)) {
        plasma_error("PLASMA_dormlq", "illegal value of K");
        return -5;
    }
    if (LDA < max(1, K)) {
        plasma_error("PLASMA_dormlq", "illegal value of LDA");
        return -7;
    }
    if (LDB < max(1, M)) {
        plasma_error("PLASMA_dormlq", "illegal value of LDB");
        return -10;
    }
    /* Quick return - currently NOT equivalent to LAPACK's:
     * CALL DLASET( 'Full', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB ) */
    if (min(M, min(N, K)) == 0)
        return PLASMA_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_DGELS, M, K, N);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dormlq", "plasma_tune() failed");
        return status;
    }

    /* Set MT, NT & NTRHS */
    NB   = PLASMA_NB;
    IB   = PLASMA_IB;
    IBNB = IB*NB;
    KT   = ( K%NB==0) ? (K /NB) : (K /NB+1);
    NT   = (An%NB==0) ? (An/NB) : (An/NB+1);

    plasma_sequence_create(plasma, &sequence);

    if (plasma->householder == PLASMA_FLAT_HOUSEHOLDER) {
        descT = plasma_desc_init(
            PlasmaRealDouble,
            IB, NB, IBNB,
            KT*IB, NT*NB, 0, 0, KT*IB, NT*NB);
    }
    else {
        /* Double the size of T to accomodate the tree reduction phase */
        descT = plasma_desc_init(
            PlasmaRealDouble,
            IB, NB, IBNB,
            KT*IB, 2*NT*NB, 0, 0, KT*IB, 2*NT*NB);
    }
    descT.mat = T;

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dooplap2tile( descA, A, NB, NB, LDA, An, 0, 0, K, An, plasma_desc_mat_free(&(descA)) );
        plasma_dooplap2tile( descB, B, NB, NB, LDB, N,  0, 0, M, N,  plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)));
    } else {
        plasma_diplap2tile( descA, A, NB, NB, LDA, An, 0, 0, K, An);
        plasma_diplap2tile( descB, B, NB, NB, LDB, N,  0, 0, M, N);
    }

    /* Call the tile interface */
    PLASMA_dormlq_Tile_Async(
        side, trans, &descA, &descT, &descB, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dooptile2lap( descB, B, NB, NB, LDB, N );
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
        plasma_desc_mat_free(&descB);
    } else {
        plasma_diptile2lap( descA, A, NB, NB, LDA, An );
        plasma_diptile2lap( descB, B, NB, NB, LDB, N );
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
 *  PLASMA_dormlq_Tile - overwrites the general M-by-N matrix C with Q*C, where Q is an orthogonal
 *  matrix (unitary in the complex case) defined as the product of elementary reflectors returned
 *  by PLASMA_dgelqf_Tile Q is of order M.
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
 *          Details of the LQ factorization of the original matrix A as returned by PLASMA_dgelqf.
 *
 * @param[in] T
 *          Auxiliary factorization data, computed by PLASMA_dgelqf.
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
 * @sa PLASMA_dormlq
 * @sa PLASMA_dormlq_Tile_Async
 * @sa PLASMA_cunmlq_Tile
 * @sa PLASMA_dormlq_Tile
 * @sa PLASMA_sormlq_Tile
 * @sa PLASMA_dgelqf_Tile
 *
 ******************************************************************************/
int PLASMA_dormlq_Tile(PLASMA_enum side, PLASMA_enum trans,
                       PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dormlq_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_dormlq_Tile_Async(side, trans, A, T, B, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup double_Tile_Async
 *
 *  Non-blocking equivalent of PLASMA_dormlq_Tile().
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
 * @sa PLASMA_dormlq
 * @sa PLASMA_dormlq_Tile
 * @sa PLASMA_cunmlq_Tile_Async
 * @sa PLASMA_dormlq_Tile_Async
 * @sa PLASMA_sormlq_Tile_Async
 * @sa PLASMA_dgelqf_Tile_Async
 *
 ******************************************************************************/
int PLASMA_dormlq_Tile_Async(PLASMA_enum side, PLASMA_enum trans,
                             PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA = *A;
    PLASMA_desc descT = *T;
    PLASMA_desc descB = *B;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dormlq_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_dormlq_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_dormlq_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dormlq_Tile", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (plasma_desc_check(&descT) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dormlq_Tile", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (plasma_desc_check(&descB) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dormlq_Tile", "invalid third descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (descA.nb != descA.mb || descB.nb != descB.mb) {
        plasma_error("PLASMA_dormlq_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ( (side != PlasmaLeft) && (side != PlasmaRight) ) {
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ( (trans != PlasmaTrans) && (trans != PlasmaNoTrans) ){
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
            plasma_parallel_call_7(plasma_pdormlq,
                PLASMA_enum, side,
                PLASMA_enum, trans,
                PLASMA_desc, descA,
                PLASMA_desc, descB,
                PLASMA_desc, descT,
                PLASMA_sequence*, sequence,
                PLASMA_request*, request);
        } else {
            plasma_dynamic_call_7(plasma_pdormlq,
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
        plasma_dynamic_call_8(plasma_pdormlqrh,
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
