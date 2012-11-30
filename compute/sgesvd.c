/**
 *
 * @file sgesvd.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Azzam Haidar
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @generated s Thu Sep 15 12:09:19 2011
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup float
 *
 *  PLASMA_sgesvd - computes the singular value decomposition (SVD) of a complex
 *  M-by-N matrix A, optionally computing the left and/or right singular
 *  vectors. The SVD is written
 *
 *       A = U * SIGMA * transpose(V)
 *
 *  where SIGMA is an M-by-N matrix which is zero except for its
 *  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
 *  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
 *  are the singular values of A; they are real and non-negative, and
 *  are returned in descending order.  The first min(m,n) columns of
 *  U and V are the left and right singular vectors of A.
 *
 *  Note that the routine returns V**T, not V.
 *  Not LAPACK Compliant for now!
 *  Note: Only PlasmaNoVec supported!
 *******************************************************************************
 *
 * @param[in] jobu
 *          Specifies options for computing all or part of the matrix U.
 *          Intended usage:
 *          = PlasmaVec: all M columns of U are returned in array U;
 *          = PlasmaNoVec: no columns of U (no left singular vectors) are
 *                     computed.
 *          Note: Only PlasmaNoVec supported!
 *
 * @param[in] jobvt
 *          Specifies options for computing all or part of the matrix V**T.
 *          Intended usage:
 *          = PlasmaVec: all M columns of U are returned in array U;
 *          = PlasmaNoVec: no columns of U (no left singular vectors) are
 *                     computed.
 *          Note: Only PlasmaNoVec supported!
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit,
 *          if JOBU = 'O',  A is overwritten with the first min(m,n)
 *                          columns of U (the left singular vectors,
 *                          stored columnwise);
 *          if JOBVT = 'O', A is overwritten with the first min(m,n)
 *                          rows of V**T (the right singular vectors,
 *                          stored rowwise);
 *          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
 *                          are destroyed.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] S
 *          The real singular values of A, sorted so that S(i) >= S(i+1).
 *
 * @param[out] U
 *          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
 *          If JOBU = 'A', U contains the M-by-M unitary matrix U;
 *          if JOBU = 'S', U contains the first min(m,n) columns of U
 *          (the left singular vectors, stored columnwise);
 *          if JOBU = 'N' or 'O', U is not referenced.
 *
 * @param[in] LDU
 *          The leading dimension of the array U.  LDU >= 1; if
 *          JOBU = 'S' or 'A', LDU >= M.
 *
 * @param[out] VT
 *         If JOBVT = 'A', VT contains the N-by-N unitary matrix
 *         V**T;
 *         if JOBVT = 'S', VT contains the first min(m,n) rows of
 *         V**T (the right singular vectors, stored rowwise);
 *         if JOBVT = 'N' or 'O', VT is not referenced.
 *
 * @param[in] LDVT
 *         The leading dimension of the array VT.  LDVT >= 1; if
 *         JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
 *
 * @param[in, out] descT
 *          On entry, descriptor as return by PLASMA_Alloc_Workspace_sgesvd
 *          On exit, contains auxiliary factorization data.
 *
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_sgesvd_Tile
 * @sa PLASMA_sgesvd_Tile_Async
 * @sa PLASMA_cgesvd
 * @sa PLASMA_dgesvd
 * @sa PLASMA_sgesvd
 *
 ******************************************************************************/
int PLASMA_sgesvd(PLASMA_enum jobu, PLASMA_enum jobvt, int M, int N,
                  float *A, int LDA,
                  float *S,
                  float *U, int LDU,
                  float *VT, int LDVT,
                  PLASMA_desc *descT)
{
    int NB, IB, IBNB, minMN, MT, NT, minMTNT;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA, descU, descVT;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_sgesvd", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    
    /* Tune NB & IB depending on M & N; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_SGESVD, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sgesvd", "plasma_tune() failed");
        return status;
    }

    /* Set MT, NT */
    NB    = PLASMA_NB;
    IB    = PLASMA_IB;
    IBNB  = IB*NB;
    MT    = (M%NB==0) ? (M/NB) : (M/NB+1);
    NT    = (N%NB==0) ? (N/NB) : (N/NB+1);
    minMN = min(M,N);
    minMTNT = min(MT,NT);

    /* Check input arguments */
    if (jobu != PlasmaNoVec  && jobu !=PlasmaVec) {
        plasma_error("PLASMA_sgesvd", "illegal value of jobu");
        return -1;
    }
    if (jobvt != PlasmaNoVec && jobvt != PlasmaVec) {
        plasma_error("PLASMA_sgesvd", "illegal value of jobvt");
        return -2;
    }
    if (M < 0) {
        plasma_error("PLASMA_sgesvd", "illegal value of M");
        return -3;
    }
    if (N < 0) {
        plasma_error("PLASMA_sgesvd", "illegal value of N");
        return -4;
    }
    if (LDA < max(1, M)) {
        plasma_error("PLASMA_sgesvd", "illegal value of LDA");
        return -6;
    }
    if (LDU < 1) {
        plasma_error("PLASMA_sgesvd", "illegal value of LDU");
        return -9;
    }
    if (LDVT < 1) {
        plasma_error("PLASMA_sgesvd", "illegal value of LDVT");
        return -11;
    }
    if ( (plasma_desc_check(descT) != PLASMA_SUCCESS) || 
         ( descT->m != MT*IB ) || (descT->n != NT*NB) ) {
        plasma_error("PLASMA_sgesvd", "invalid T descriptor");
        return -12;
    }
    /* Quick return */
    if (min(M, N) == 0) {
        return PLASMA_SUCCESS;
    }

    if (jobu == PlasmaVec) {
        plasma_error("PLASMA_sgesvd", "computing the singular vectors is not supported in this version");
        return -1;
    }
    if (jobvt == PlasmaVec) {
        plasma_error("PLASMA_sgesvd", "computing the singular vectors is not supported in this version");
        return -2;
    }

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_sooplap2tile( descA,   A, NB, NB,  LDA, N, 0, 0, M, N, plasma_desc_mat_free(&(descA)) );
        if (jobu == PlasmaVec){
            plasma_sooplap2tile( descU,   U, NB, NB,  LDU, M, 0, 0, M, M, plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descU)));
        }
        if (jobvt == PlasmaVec){
            plasma_sooplap2tile( descVT, VT, NB, NB, LDVT, N, 0, 0, N, N, plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descU)); plasma_desc_mat_free(&(descVT)));
        }
    } else {
        plasma_siplap2tile( descA,   A, NB, NB,  LDA, N, 0, 0, M, N);
        if (jobu == PlasmaVec){
            plasma_siplap2tile( descU,   U, NB, NB,  LDU, M, 0, 0, M, M);
        }
        if (jobvt == PlasmaVec){
            plasma_siplap2tile( descVT, VT, NB, NB, LDVT, N, 0, 0, N, N);
        }
    }

    /* Call the tile interface */
    PLASMA_sgesvd_Tile_Async(jobu, jobvt, &descA, S, &descU, &descVT, descT, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_sooptile2lap( descA,   A, NB, NB,  LDA, N );
        if (jobu == PlasmaVec){
            plasma_sooptile2lap( descU,   U, NB, NB,  LDU, M );
        }
        if (jobvt == PlasmaVec){
            plasma_sooptile2lap( descVT, VT, NB, NB, LDVT, N );
        }
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
        if (jobu == PlasmaVec){
            plasma_desc_mat_free(&descU);
        }
        if (jobvt == PlasmaVec){
            plasma_desc_mat_free(&descVT);
        }
    } else {
        plasma_siptile2lap( descA,   A, NB, NB,  LDA, N );
        if (jobu == PlasmaVec){
            plasma_siptile2lap( descU,   U, NB, NB,  LDU, M );
        }
        if (jobvt == PlasmaVec){
            plasma_siptile2lap( descVT, VT, NB, NB, LDVT, N );
        }
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
 *  PLASMA_sgesvd_Tile - computes the singular value decomposition (SVD) of a complex
 *  M-by-N matrix A, optionally computing the left and/or right singular
 *  vectors.
 *  Tile equivalent of PLASMA_sgesvd().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] jobu
 *          Specifies options for computing all or part of the matrix U.
 *          Intended usage:
 *          = PlasmaVec: all M columns of U are returned in array U;
 *          = PlasmaNoVec: no columns of U (no left singular vectors) are
 *                     computed.
 *
 * @param[in] jobvt
 *          Specifies options for computing all or part of the matrix V**T.
 *          Intended usage:
 *          = PlasmaVec: all M columns of U are returned in array U;
 *          = PlasmaNoVec: no columns of U (no left singular vectors) are
 *                     computed.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit,
 *          if JOBU = 'O',  A is overwritten with the first min(m,n)
 *                          columns of U (the left singular vectors,
 *                          stored columnwise);
 *          if JOBVT = 'O', A is overwritten with the first min(m,n)
 *                          rows of V**T (the right singular vectors,
 *                          stored rowwise);
 *          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
 *                          are destroyed.
 *
 * @param[out] S
 *          The singular values of A, sorted so that S(i) >= S(i+1).
 *
 * @param[out] U
 *          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
 *          If JOBU = 'A', U contains the M-by-M unitary matrix U;
 *          if JOBU = 'S', U contains the first min(m,n) columns of U
 *          (the left singular vectors, stored columnwise);
 *          if JOBU = 'N' or 'O', U is not referenced.
 *
 * @param[out] VT
 *         If JOBVT = 'A', VT contains the N-by-N unitary matrix
 *         V**T;
 *         if JOBVT = 'S', VT contains the first min(m,n) rows of
 *         V**T (the right singular vectors, stored rowwise);
 *         if JOBVT = 'N' or 'O', VT is not referenced.
 *
 * @param[out] T
 *         On exit, auxiliary factorization data.
 *
 *******************************************************************************
 *
 * @return
 *          \return PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_sgesvd
 * @sa PLASMA_sgesvd_Tile_Async
 * @sa PLASMA_cgesvd_Tile
 * @sa PLASMA_dgesvd_Tile
 * @sa PLASMA_sgesvd_Tile
 *
 ******************************************************************************/
int PLASMA_sgesvd_Tile(PLASMA_enum jobu, PLASMA_enum jobvt, PLASMA_desc *A,
                      float *S, PLASMA_desc *U, PLASMA_desc *VT, PLASMA_desc *T)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_sgesvd_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_sgesvd_Tile_Async(jobu, jobvt, A, S, U, VT, T, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup float_Tile_Async
 *
 *  PLASMA_sgesvd_Tile_Async - computes the singular value decomposition (SVD) of a complex
 *  M-by-N matrix A, optionally computing the left and/or right singular
 *  vectors.
 *  Non-blocking equivalent of PLASMA_sgesvd_Tile().
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
 * @sa PLASMA_sgesvd
 * @sa PLASMA_sgesvd_Tile
 * @sa PLASMA_cgesvd_Tile_Async
 * @sa PLASMA_dgesvd_Tile_Async
 * @sa PLASMA_sgesvd_Tile_Async
 *
 ******************************************************************************/
int PLASMA_sgesvd_Tile_Async(PLASMA_enum jobu, PLASMA_enum jobvt, PLASMA_desc *A,
                             float *S, PLASMA_desc *U, PLASMA_desc *VT, PLASMA_desc *T,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA  = *A;
    PLASMA_desc descT  = *T;
    float *E;
    int minMN = min(descA.m, descA.n);
    int NCVT = 0;
    int NRU = 0;
    int NCC = 0;

    plasma_context_t *plasma;

    plasma = plasma_context_self();

    if (jobu != PlasmaNoVec  && jobu !=PlasmaVec) {
        plasma_error("PLASMA_sgesvd_Tile_Async", "illegal value of jobu");
        return PLASMA_ERR_NOT_SUPPORTED;
    }
    if (jobvt != PlasmaNoVec && jobvt != PlasmaVec) {
        plasma_error("PLASMA_sgesvd_Tile_Async", "illegal value of jobvt");
        return PLASMA_ERR_NOT_SUPPORTED;
    }
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_sgesvd_Tile_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_sgesvd_Tile_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_sgesvd_Tile_Async", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sgesvd_Tile_Async", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ((jobu != PlasmaNoVec) && (plasma_desc_check(U) != PLASMA_SUCCESS)) {
        plasma_error("PLASMA_sgesvd_Tile_Async", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ((jobvt != PlasmaNoVec) && (plasma_desc_check(VT) != PLASMA_SUCCESS) ) {
        plasma_error("PLASMA_sgesvd_Tile_Async", "invalid third descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (plasma_desc_check(&descT) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sgesvd_Tile_Async", "invalid fourth descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_sgesvd_Tile_Async", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (( (jobu != PlasmaNoVec) && (U->nb != U->mb) )  || ( (jobvt != PlasmaNoVec) && (VT->nb != VT->mb) )) {
        plasma_error("PLASMA_sgesvd_Tile_Async", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ((jobu == PlasmaVec) || (jobvt == PlasmaVec) ){
        plasma_error("PLASMA_sgesvd_Tile_Async", "computing the singular vectors is not supported in this version");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    
    E = (float *)plasma_shared_alloc(plasma, minMN-1, PlasmaRealDouble);
 
    /* 
     * Reduction to bidiagonal form with a two-stage approach.
     */

    /* 
     * 1: Reduction to BAND bidiagonal form
     * May be further optimized using the algo described in Trefethen
     */
    /* if (plasma->householder == PLASMA_FLAT_HOUSEHOLDER) { */
        plasma_dynamic_call_4(plasma_psgerbb,
            PLASMA_desc, descA,
            PLASMA_desc, descT,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);
    /* } */
    /* else { */
    /*     plasma_dynamic_call_4(plasma_psgerbbrh, */
    /*         PLASMA_desc, descA, */
    /*         PLASMA_desc, descT, */
    /*         PLASMA_sequence*, sequence, */
    /*         PLASMA_request*, request); */
    /* } */

    /* Build the U of the first stage */
    /* if (jobu == PlasmaVec){ */
    /*    /\* Initialize U to Identity *\/ */
    /*    plasma_dynamic_call_6(plasma_pslaset, */
    /*        PLASMA_enum, PlasmaUpperLower, */
    /*        float, 0.0, */
    /*        float, 1.0, */
    /*        PLASMA_desc, descU, */
    /*        PLASMA_sequence*, sequence, */
    /*        PLASMA_request*, request); */
    /*    /\* Accumulate the transformations from the first stage *\/ */
    /*    if (plasma->householder == PLASMA_FLAT_HOUSEHOLDER) { */
    /*        plasma_dynamic_call_6(plasma_psorgbr, */
    /*            PLASMA_enum, PlasmaLeft, */
    /*            PLASMA_desc, descA, */
    /*            PLASMA_desc, descU, */
    /*            PLASMA_desc, descT, */
    /*            PLASMA_sequence*, sequence, */
    /*            PLASMA_request*, request); */
    /*    } */
    /*    else { */
    /*        plasma_dynamic_call_6(plasma_psorgbrrh, */
    /*            PLASMA_enum, PlasmaLeft, */
    /*            PLASMA_desc, descA, */
    /*            PLASMA_desc, descU, */
    /*            PLASMA_desc, descT, */
    /*            PLASMA_sequence*, sequence, */
    /*            PLASMA_request*, request); */
    /*    } */
    /* } */
    
    /* Build the VT of the first stage */
    /* if (jobvt == PlasmaVec){ */
    /*     /\* Initialize VT to Identity *\/ */
    /*    plasma_dynamic_call_6(plasma_pslaset, */
    /*        PLASMA_enum, PlasmaUpperLower, */
    /*        float, 0.0, */
    /*        float, 1.0, */
    /*        PLASMA_desc, descVT, */
    /*        PLASMA_sequence*, sequence, */
    /*        PLASMA_request*, request); */

    /*    /\* Accumulate the transformations from the first stage *\/ */
    /*    if (plasma->householder == PLASMA_FLAT_HOUSEHOLDER) { */
    /*        plasma_dynamic_call_6(plasma_psorgbr, */
    /*            PLASMA_enum, PlasmaRight, */
    /*            PLASMA_desc, descA, */
    /*            PLASMA_desc, descVT, */
    /*            PLASMA_desc, descT, */
    /*            PLASMA_sequence*, sequence, */
    /*            PLASMA_request*, request); */
    /*    } */
    /*    else { */
    /*        plasma_dynamic_call_6(plasma_psorgbrrh, */
    /*            PLASMA_enum, PlasmaRight, */
    /*            PLASMA_desc, descA, */
    /*            PLASMA_desc, descVT, */
    /*            PLASMA_desc, descT, */
    /*            PLASMA_sequence*, sequence, */
    /*            PLASMA_request*, request); */
    /*    } */
    /* } */

    /*
     * Set the V's to zero before the 2nd stage i.e., bulge chasing
     */
    plasma_dynamic_call_5(plasma_pslaset2,
        PLASMA_enum, PlasmaLower,
        float, 0.0,
        PLASMA_desc, descA.m >= descA.n ? descA : plasma_desc_submatrix(descA, descA.mb, 0, descA.m-descA.mb, descA.n),
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    plasma_dynamic_call_5(plasma_pslaset2,
        PLASMA_enum, PlasmaUpper,
        float, 0.0,
        PLASMA_desc, descA.m >= descA.n ? plasma_desc_submatrix(descA, 0, descA.nb, descA.m, descA.n-descA.nb) : descA,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    /* 
     * 2: Reduction from BAND bidiagonal to the final condensed form
     */
    plasma_dynamic_call_7(plasma_psgbrdb,
        PLASMA_enum, descA.m >= descA.n ? PlasmaUpper : PlasmaLower,
        PLASMA_desc, descA,
        float*, S,
        float*, E,
        PLASMA_desc, descT,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    /*
     * Compute the singular values ONLY for now
     */
    plasma_dynamic_sync();
    if (descA.m >= descA.n)
       LAPACKE_sbdsqr(
              LAPACK_COL_MAJOR, lapack_const(PlasmaUpper),
              minMN, NCVT, NRU, NCC,
              S, E,
              NULL, 1, NULL, 1, NULL, 1 );
    else {
       LAPACKE_sbdsqr(
              LAPACK_COL_MAJOR, lapack_const(PlasmaLower),
              minMN, NCVT, NRU, NCC,
              S, E,
              NULL, 1, NULL, 1, NULL, 1 );
    }

    plasma_shared_free(plasma, E);

    return PLASMA_SUCCESS;
}
