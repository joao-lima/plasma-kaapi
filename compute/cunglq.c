/**
 *
 * @file cunglq.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @date 2010-11-15
 * @generated c Thu Sep 15 12:09:18 2011
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex32_t
 *
 *  PLASMA_cunglq - Generates an M-by-N matrix Q with orthonormal rows, which is defined as the
 *  first M rows of a product of the elementary reflectors returned by PLASMA_cgelqf.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix Q. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix Q. N >= M.
 *
 * @param[in] K
 *          The number of rows of elementary tile reflectors whose product defines the matrix Q.
 *          M >= K >= 0.
 *
 * @param[in] A
 *          Details of the LQ factorization of the original matrix A as returned by PLASMA_cgelqf.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in] T
 *          Auxiliary factorization data, computed by PLASMA_cgelqf.
 *
 * @param[out] Q
 *          On exit, the M-by-N matrix Q.
 *
 * @param[in] LDQ
 *          The leading dimension of the array Q. LDQ >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval PLASMA_SUCCESS <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_cunglq_Tile
 * @sa PLASMA_cunglq_Tile_Async
 * @sa PLASMA_cunglq
 * @sa PLASMA_dorglq
 * @sa PLASMA_sorglq
 * @sa PLASMA_cgelqf
 *
 ******************************************************************************/
int PLASMA_cunglq(int M, int N, int K,
                  PLASMA_Complex32_t *A, int LDA,
                  PLASMA_Complex32_t *T,
                  PLASMA_Complex32_t *Q, int LDQ)
{
    int NB, IB, IBNB, KT, NT;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA, descQ, descT;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_cunglq", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (M < 0) {
        plasma_error("PLASMA_cunglq", "illegal value of M");
        return -1;
    }
    if (N < M) {
        plasma_error("PLASMA_cunglq", "illegal value of N");
        return -2;
    }
    if (K < 0 || K > M) {
        plasma_error("PLASMA_cunglq", "illegal value of K");
        return -3;
    }
    if (LDA < max(1, M)) {
        plasma_error("PLASMA_cunglq", "illegal value of LDA");
        return -5;
    }
    if (LDQ < max(1, M)) {
        plasma_error("PLASMA_cunglq", "illegal value of LDQ");
        return -8;
    }
    /* Quick return - currently NOT equivalent to LAPACK's:
     * CALL DLASET( 'Full', MAX( M, N ), NRHS, ZERO, ZERO, B, LDQ ) */
    if (min(M, min(N, K)) == 0)
        return PLASMA_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_CGELS, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_cunglq", "plasma_tune() failed");
        return status;
    }

    /* Set MT & NT */
    NB   = PLASMA_NB;
    IB   = PLASMA_IB;
    IBNB = IB*NB;
    NT   = (N%NB==0) ? (N/NB) : (N/NB+1);
    KT   = (K%NB==0) ? (K/NB) : (K/NB+1);

    plasma_sequence_create(plasma, &sequence);

    if (plasma->householder == PLASMA_FLAT_HOUSEHOLDER) {
        descT = plasma_desc_init(
            PlasmaComplexFloat,
            IB, NB, IBNB,
            KT*IB, NT*NB, 0, 0, KT*IB, NT*NB);
    }
    else {
        /* Double the size of T to accomodate the tree reduction phase */
        descT = plasma_desc_init(
            PlasmaComplexFloat,
            IB, NB, IBNB,
            KT*IB, 2*NT*NB, 0, 0, KT*IB, 2*NT*NB);
    }
    descT.mat = T;

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_cooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, K, N, plasma_desc_mat_free(&(descA)) );
        plasma_cooplap2tile( descQ, Q, NB, NB, LDQ, N, 0, 0, M, N, plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descQ)));
    } else {
        plasma_ciplap2tile( descA, A, NB, NB, LDA, N, 0, 0, K, N);
        plasma_ciplap2tile( descQ, Q, NB, NB, LDQ, N, 0, 0, M, N);
    }

    /* Call the tile interface */
    PLASMA_cunglq_Tile_Async(&descA, &descT, &descQ, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_cooptile2lap( descQ, Q, NB, NB, LDQ, N );
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
        plasma_desc_mat_free(&descQ);
    } else {
        plasma_ciptile2lap( descA, A, NB, NB, LDA, N );
        plasma_ciptile2lap( descQ, Q, NB, NB, LDQ, N );
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
 * PLASMA_cunglq_Tile - Generates an M-by-N matrix Q with orthonormal rows, which is defined as the
 * first M rows of a product of the elementary reflectors returned by PLASMA_cgelqf.
 * All matrices are passed through descriptors. All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          Details of the LQ factorization of the original matrix A as returned by PLASMA_cgelqf.
 *
 * @param[in] T
 *          Auxiliary factorization data, computed by PLASMA_cgelqf.
 *
 * @param[out] B
 *          On exit, the M-by-N matrix Q.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_cunglq
 * @sa PLASMA_cunglq_Tile_Async
 * @sa PLASMA_cunglq_Tile
 * @sa PLASMA_dorglq_Tile
 * @sa PLASMA_sorglq_Tile
 * @sa PLASMA_cgelqf_Tile
 *
 ******************************************************************************/
int PLASMA_cunglq_Tile(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_cunglq_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_cunglq_Tile_Async(A, T, B, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex32_t_Tile_Async
 *
 *  Non-blocking equivalent of PLASMA_cunglq_Tile().
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
 * @sa PLASMA_cunglq
 * @sa PLASMA_cunglq_Tile
 * @sa PLASMA_cunglq_Tile_Async
 * @sa PLASMA_dorglq_Tile_Async
 * @sa PLASMA_sorglq_Tile_Async
 * @sa PLASMA_cgelqf_Tile_Async
 *
 ******************************************************************************/
int PLASMA_cunglq_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *Q,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA = *A;
    PLASMA_desc descT = *T;
    PLASMA_desc descQ = *Q;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_cunglq_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_cunglq_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_cunglq_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_cunglq_Tile", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (plasma_desc_check(&descT) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_cunglq_Tile", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (plasma_desc_check(&descQ) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_cunglq_Tile", "invalid third descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (descA.nb != descA.mb || descQ.nb != descQ.mb) {
        plasma_error("PLASMA_cunglq_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Quick return - currently NOT equivalent to LAPACK's:
     * CALL DLASET( 'Full', MAX( M, N ), NRHS, ZERO, ZERO, Q, LDQ ) */
/*
    if (min(M, N) == 0)
        return PLASMA_SUCCESS;
*/
    if (plasma->householder == PLASMA_FLAT_HOUSEHOLDER) {
        plasma_dynamic_call_5(plasma_pcunglq,
            PLASMA_desc, descA,
            PLASMA_desc, descQ,
            PLASMA_desc, descT,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);
    }
    else {
        plasma_dynamic_call_6(plasma_pcunglqrh,
            PLASMA_desc, descA,
            PLASMA_desc, descQ,
            PLASMA_desc, descT,
            PLASMA_enum, PLASMA_RHBLK,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);
    }

    return PLASMA_SUCCESS;
}
