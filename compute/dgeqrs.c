/**
 *
 * @file dgeqrs.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Jakub Kurzak
 * @date 2010-11-15
 * @generated d Thu Sep 15 12:09:14 2011
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup double
 *
 *  PLASMA_dgeqrs - Compute a minimum-norm solution min || A*X - B || using the RQ factorization
 *  A = R*Q computed by PLASMA_dgeqrf.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= M >= 0.
 *
 * @param[in] NRHS
 *          The number of columns of B. NRHS >= 0.
 *
 * @param[in,out] A
 *          Details of the QR factorization of the original matrix A as returned by PLASMA_dgeqrf.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= M.
 *
 * @param[in] T
 *          Auxiliary factorization data, computed by PLASMA_dgeqrf.
 *
 * @param[in,out] B
 *          On entry, the m-by-nrhs right hand side matrix B.
 *          On exit, the n-by-nrhs solution matrix X.
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
 * @sa PLASMA_dgeqrs_Tile
 * @sa PLASMA_dgeqrs_Tile_Async
 * @sa PLASMA_cgeqrs
 * @sa PLASMA_dgeqrs
 * @sa PLASMA_sgeqrs
 * @sa PLASMA_dgeqrf
 *
 *******************************************************************************/
int PLASMA_dgeqrs(int M, int N, int NRHS,
                  double *A, int LDA,
                  double *T,
                  double *B, int LDB)
{
    int NB, IB, IBNB, MT, NT;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA, descB, descT;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dgeqrs", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        plasma_error("PLASMA_dgeqrs", "illegal value of M");
        return -1;
    }
    if (N < 0 || N > M) {
        plasma_error("PLASMA_dgeqrs", "illegal value of N");
        return -2;
    }
    if (NRHS < 0) {
        plasma_error("PLASMA_dgeqrs", "illegal value of N");
        return -3;
    }
    if (LDA < max(1, M)) {
        plasma_error("PLASMA_dgeqrs", "illegal value of LDA");
        return -5;
    }
    if (LDB < max(1, max(1, M))) {
        plasma_error("PLASMA_dgeqrs", "illegal value of LDB");
        return -8;
    }
    /* Quick return */
    if (min(M, min(N, NRHS)) == 0) {
        return PLASMA_SUCCESS;
    }

    /* Tune NB & IB depending on M, N & NRHS; Set NBNBSIZE */
    status = plasma_tune(PLASMA_FUNC_DGELS, M, N, NRHS);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dgeqrs", "plasma_tune() failed");
        return status;
    }

    /* Set MT, NT & NTRHS */
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
        plasma_dooplap2tile( descA, A, NB, NB, LDA, N,    0, 0, M, N   , plasma_desc_mat_free(&(descA)) );
        plasma_dooplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, M, NRHS, plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)));
    } else {
        plasma_diplap2tile( descA, A, NB, NB, LDA, N,    0, 0, M, N   );
        plasma_diplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, M, NRHS);
    }

    /* Call the tile interface */
    PLASMA_dgeqrs_Tile_Async(&descA, &descT, &descB, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dooptile2lap( descA, A, NB, NB, LDA, N    );
        plasma_dooptile2lap( descB, B, NB, NB, LDB, NRHS );
    plasma_dynamic_sync();
    plasma_desc_mat_free(&descA);
    plasma_desc_mat_free(&descB);
    } else {
        plasma_diptile2lap( descA, A, NB, NB, LDA, N    );
        plasma_diptile2lap( descB, B, NB, NB, LDB, NRHS );
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
 *  PLASMA_dgeqrs_Tile - Computes a minimum-norm solution using the tile QR factorization.
 *  Tile equivalent of PLASMA_dgeqrf().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          Details of the QR factorization of the original matrix A as returned by PLASMA_dgeqrf.
 *
 * @param[in] T
 *          Auxiliary factorization data, computed by PLASMA_dgeqrf.
 *
 * @param[in,out] B
 *          On entry, the m-by-nrhs right hand side matrix B.
 *          On exit, the n-by-nrhs solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_dgeqrs
 * @sa PLASMA_dgeqrs_Tile_Async
 * @sa PLASMA_cgeqrs_Tile
 * @sa PLASMA_dgeqrs_Tile
 * @sa PLASMA_sgeqrs_Tile
 * @sa PLASMA_dgeqrf_Tile
 *
 ******************************************************************************/
int PLASMA_dgeqrs_Tile(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dgeqrs_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_dgeqrs_Tile_Async(A, T, B, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup double_Tile_Async
 *
 *  PLASMA_dgeqrs_Tile_Async - Computes a minimum-norm solution using the tile
 *  QR factorization.
 *  Non-blocking equivalent of PLASMA_dgeqrs_Tile().
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
 * @sa PLASMA_dgeqrs
 * @sa PLASMA_dgeqrs_Tile
 * @sa PLASMA_cgeqrs_Tile_Async
 * @sa PLASMA_dgeqrs_Tile_Async
 * @sa PLASMA_sgeqrs_Tile_Async
 * @sa PLASMA_dgeqrf_Tile_Async
 *
 ******************************************************************************/
int PLASMA_dgeqrs_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA = *A;
    PLASMA_desc descT = *T;
    PLASMA_desc descB = *B;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dgeqrs_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_dgeqrs_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_dgeqrs_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dgeqrs_Tile", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (plasma_desc_check(&descT) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dgeqrs_Tile", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (plasma_desc_check(&descB) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dgeqrs_Tile", "invalid third descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (descA.nb != descA.mb || descB.nb != descB.mb) {
        plasma_error("PLASMA_dgeqrs_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
/*
    if (min(M, min(N, NRHS)) == 0) {
        return PLASMA_SUCCESS;
    }
*/
    if (plasma->householder == PLASMA_FLAT_HOUSEHOLDER) {
        plasma_parallel_call_7(plasma_pdormqr,
            PLASMA_enum, PlasmaLeft,
            PLASMA_enum, PlasmaTrans,
            PLASMA_desc, descA,
            PLASMA_desc, descB,
            PLASMA_desc, descT,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);
    }
    else {
        plasma_dynamic_call_8(plasma_pdormqrrh,
            PLASMA_enum, PlasmaLeft,
            PLASMA_enum, PlasmaTrans,
            PLASMA_desc, descA,
            PLASMA_desc, descB,
            PLASMA_desc, descT,
            PLASMA_enum, PLASMA_RHBLK,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);
    }

    plasma_parallel_call_9(plasma_pdtrsm,
        PLASMA_enum, PlasmaLeft,
        PLASMA_enum, PlasmaUpper,
        PLASMA_enum, PlasmaNoTrans,
        PLASMA_enum, PlasmaNonUnit,
        double, 1.0,
        PLASMA_desc, plasma_desc_submatrix(descA, 0, 0, descA.n, descA.n),
        PLASMA_desc, plasma_desc_submatrix(descB, 0, 0, descA.n, descB.n),
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}
