/**
 *
 * @file dgelqf.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Jakub Kurzak
 * @author Dulceneia Becker
 * @date 2010-11-15
 * @generated d Thu Sep 15 12:09:14 2011
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup double
 *
 *  PLASMA_dgelqf - Computes the tile LQ factorization of a complex M-by-N matrix A: A = L * Q.
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
 *          On entry, the M-by-N matrix A.
 *          On exit, the elements on and below the diagonal of the array contain the m-by-min(M,N)
 *          lower trapezoidal matrix L (L is lower triangular if M <= N); the elements above the
 *          diagonal represent the unitary matrix Q as a product of elementary reflectors, stored
 *          by tiles.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] T
 *          On exit, auxiliary factorization data, required by PLASMA_dgelqs to solve the system
 *          of equations.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_dgelqf_Tile
 * @sa PLASMA_dgelqf_Tile_Async
 * @sa PLASMA_cgelqf
 * @sa PLASMA_dgelqf
 * @sa PLASMA_sgelqf
 * @sa PLASMA_dgelqs
 *
 ******************************************************************************/
int PLASMA_dgelqf(int M, int N,
                  double *A, int LDA,
                  double *T)
{
    int NB, IB, IBNB, MT, NT;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA, descT;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dgelqf", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        plasma_error("PLASMA_dgelqf", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        plasma_error("PLASMA_dgelqf", "illegal value of N");
        return -2;
    }
    if (LDA < max(1, M)) {
        plasma_error("PLASMA_dgelqf", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (min(M, N) == 0)
        return PLASMA_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNBSIZE */
    status = plasma_tune(PLASMA_FUNC_DGELS, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dgelqf", "plasma_tune() failed");
        return status;
    }

    /* Set MT & NT */
    NB   = PLASMA_NB;
    IB   = PLASMA_IB;
    IBNB = IB*NB;
    MT   = (M%NB==0) ? (M/NB) : (M/NB+1);
    NT   = (N%NB==0) ? (N/NB) : (N/NB+1);

    plasma_sequence_create(plasma, &sequence);

     if (plasma->householder == PLASMA_FLAT_HOUSEHOLDER) {
        descT = plasma_desc_init(
            PlasmaRealDouble,
            IB, NB, IBNB,
            MT*IB, NT*NB, 0, 0, MT*IB, NT*NB);
    }
    else {
        /* Double the size of T to accomodate the tree reduction phase */
        descT = plasma_desc_init(
            PlasmaRealDouble,
            IB, NB, IBNB,
            MT*IB, 2*NT*NB, 0, 0, MT*IB, 2*NT*NB);
    }
    descT.mat = T;

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, M, N, plasma_desc_mat_free(&(descA)) );
    } else {
        plasma_diplap2tile( descA, A, NB, NB, LDA, N, 0, 0, M, N);
    }

    /* Call the tile interface */
    PLASMA_dgelqf_Tile_Async(&descA, &descT, sequence, &request);

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
 *  PLASMA_dgelqf_Tile - Computes the tile LQ factorization of a matrix.
 *  Tile equivalent of PLASMA_dgelqf().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit, the elements on and below the diagonal of the array contain the m-by-min(M,N)
 *          lower trapezoidal matrix L (L is lower triangular if M <= N); the elements above the
 *          diagonal represent the unitary matrix Q as a product of elementary reflectors, stored
 *          by tiles.
 *
 * @param[out] T
 *          On exit, auxiliary factorization data, required by PLASMA_dgelqs to solve the system
 *          of equations.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_dgelqf
 * @sa PLASMA_dgelqf_Tile_Async
 * @sa PLASMA_cgelqf_Tile
 * @sa PLASMA_dgelqf_Tile
 * @sa PLASMA_sgelqf_Tile
 * @sa PLASMA_dgelqs_Tile
 *
 ******************************************************************************/
int PLASMA_dgelqf_Tile(PLASMA_desc *A, PLASMA_desc *T)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dgelqf_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_dgelqf_Tile_Async(A, T, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup double_Tile_Async
 *
 *  PLASMA_dgelqf_Tile_Async - Computes the tile LQ factorization of a matrix.
 *  Non-blocking equivalent of PLASMA_dgelqf_Tile().
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
 * @sa PLASMA_dgelqf
 * @sa PLASMA_dgelqf_Tile
 * @sa PLASMA_cgelqf_Tile_Async
 * @sa PLASMA_dgelqf_Tile_Async
 * @sa PLASMA_sgelqf_Tile_Async
 * @sa PLASMA_dgelqs_Tile_Async
 *
 ******************************************************************************/
int PLASMA_dgelqf_Tile_Async(PLASMA_desc *A, PLASMA_desc *T,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA = *A;
    PLASMA_desc descT = *T;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dgelqf_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_dgelqf_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_dgelqf_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dgelqf_Tile", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (plasma_desc_check(&descT) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dgelqf_Tile", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_dgelqf_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
/*
    if (min(M, N) == 0)
        return PLASMA_SUCCESS;
*/
    if (plasma->householder == PLASMA_FLAT_HOUSEHOLDER) {
        plasma_parallel_call_4(plasma_pdgelqf,
            PLASMA_desc, descA,
            PLASMA_desc, descT,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);
    }
    else {
        plasma_dynamic_call_5(plasma_pdgelqfrh,
            PLASMA_desc, descA,
            PLASMA_desc, descT,
            PLASMA_enum, PLASMA_RHBLK,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);
    }

    return PLASMA_SUCCESS;
}
