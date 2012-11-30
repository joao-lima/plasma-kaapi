/**
 * @file dsytrd.c
 *
 *  PLASMA computational routines
 *  Release Date: November, 15th 2009
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Azzam Haidar
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @generated d Thu Sep 15 12:09:19 2011
 *
 **/
#include "common.h"
/***************************************************************************//**
 *
 * @ingroup double
 *
 *  PLASMA_dsytrd - reduces a complex Hermitian matrix A to real symmetric
 *  tridiagonal form S using a two-stage approach
 *  First stage: reduction to band tridiagonal form (unitary Q1);
 *  Second stage: reduction from band to tridiagonal form (unitary
 *  Q2).  Let Q = Q1 * Q2 be the global unitary transformation; Q**T *
 *  A * Q = S.
 *  Not LAPACK compliant as A does not contain the T elements
 *  Note: Only PlasmaNoVec supported!
 *
 *******************************************************************************
 *
 * @param[in] jobz
 *          Intended usage:
 *          = PlasmaNoVec: computes eigenvalues only;
 *          = PlasmaVec: computes eigenvalues and eigenvectors.
 *          Note: Only PlasmaNoVec supported!
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or
 *          lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the symmetric (or Hermitian) matrix A.
 *          If uplo = PlasmaUpper, the leading N-by-N upper triangular
 *          part of A contains the upper triangular part of the matrix
 *          A, and the strictly lower triangular part of A is not
 *          referenced.
 *          If uplo = PlasmaLower, the leading N-by-N lower triangular
 *          part of A contains the lower triangular part of the matrix
 *          A, and the strictly upper triangular part of A is not
 *          referenced.
 *          On exit, the lower triangle (if uplo = PlasmaLower) or the
 *          upper triangle (if uplo = PlasmaUpper) of A, including the
 *          diagonal, is destroyed.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[out] D
 *          On exit, the diagonal elements of the tridiagonal matrix: 
 *          D(i) = A(i,i).
 *
 * @param[out] E
 *          On exit, he off-diagonal elements of the tridiagonal matrix: 
 *          E(i) = A(i,i+1) if uplo = PlasmaUpper, E(i) = A(i+1,i) if uplo = PlasmaLower.
 *
 * @param[in, out] descT
 *          On entry, descriptor as return by PLASMA_Alloc_Workspace_dsyev
 *          On exit, contains auxiliary factorization data.
 *
 * @param[out] Q
 *          On exit, if jobz = PlasmaVec and info = 0, the eigenvectors.
 *
 * @param[in] LDQ
 *          The leading dimension of the array Q. LDQ >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if INFO = i, the algorithm failed to converge; i
 *               off-diagonal elements of an intermediate tridiagonal
 *               form did not converge to zero.
 *
 *******************************************************************************
 *
 * @sa PLASMA_dsytrd_Tile
 * @sa PLASMA_dsytrd_Tile_Async
 * @sa PLASMA_chetrd
 * @sa PLASMA_dsytrd
 * @sa PLASMA_ssytrd
 *
 ******************************************************************************/
int PLASMA_dsytrd(PLASMA_enum jobz, PLASMA_enum uplo, int N,
                 double *A, int LDA,
                 double *D,
                 double *E,
                 PLASMA_desc *descT,
                 double *Q, int LDQ)
{
    int NB, IB, IBNB, NT;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA, descQ;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_error("PLASMA_dsytrd", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    /* Tune NB & IB depending on N; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_DSYTRD, N, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dsytrd", "plasma_tune() failed");
        return status;
    }

    /* Set NT */
    NB    = PLASMA_NB;
    IB    = PLASMA_IB;
    IBNB  = IB*NB;
    NT    = (N%NB==0) ? (N/NB) : (N/NB+1);

    /* Check input arguments */
    if (jobz != PlasmaNoVec && jobz != PlasmaVec) {
        plasma_error("PLASMA_dsytrd", "illegal value of jobz");
        return -1;
    }
    if (uplo != PlasmaLower && uplo != PlasmaUpper) {
        plasma_error("PLASMA_dsytrd", "illegal value of uplo");
        return -2;
    }
    if (N < 0) {
        plasma_error("PLASMA_dsytrd", "illegal value of N");
        return -3;
    }
    if (LDA < max(1, N)) {
        plasma_error("PLASMA_dsytrd", "illegal value of LDA");
        return -5;
    }
    if ( (plasma_desc_check(descT) != PLASMA_SUCCESS) || 
         ( descT->m != NT*IB ) || (descT->n != NT*NB) ) {
        plasma_error("PLASMA_dsytrd", "invalid T descriptor");
        return -8;
    }
    if (LDQ < max(1, N)) {
        plasma_error("PLASMA_dsytrd", "illegal value of LDQ");
        return -10;
    }
    /* Quick return */
    if (N == 0)
        return PLASMA_SUCCESS;

    if (jobz == PlasmaVec) {
        plasma_error("PLASMA_dsytrd", "computing the eigenvectors is not supported in this version");
        return -1;
    }

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dooplap2tile( descA, A, NB, NB, LDA, N,    0, 0, N, N   , plasma_desc_mat_free(&(descA)) );
        if (jobz == PlasmaVec) {
           plasma_dooplap2tile( descQ, Q, NB, NB, LDQ, N,    0, 0, N, N   , plasma_desc_mat_free(&(descQ)) );
        }
    } else {
        plasma_diplap2tile( descA, A, NB, NB, LDA, N,    0, 0, N, N   );
        if (jobz == PlasmaVec)
           plasma_diplap2tile( descQ, Q, NB, NB, LDQ, N,    0, 0, N, N   );
    }

    /* Call the tile interface */
    PLASMA_dsytrd_Tile_Async(jobz, uplo, &descA, D, E, descT, &descQ, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dooptile2lap( descA, A, NB, NB, LDA, N    );
        if (jobz == PlasmaVec) {
           plasma_dooptile2lap( descQ, Q, NB, NB, LDQ, N    );
        }
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
        if (jobz == PlasmaVec)
           plasma_desc_mat_free(&descQ);
    } else {
        plasma_diptile2lap( descA, A, NB, NB, LDA, N    );
        if (jobz == PlasmaVec)
           plasma_diptile2lap( descQ, Q, NB, NB, LDQ, N    );
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
 *  PLASMA_dsytrd_Tile - reduces a complex Hermitian matrix A to real symmetric
 *  tridiagonal form S using a two-stage approach
 *  First stage: reduction to band tridiagonal form (unitary Q1);
 *  Second stage: reduction from band to tridiagonal form (unitary Q2).
 *  Let Q = Q1 * Q2 be the global unitary transformation;
 *  Q**T * A * Q = S.
 *  Tile equivalent of PLASMA_dsytrd().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] jobz
 *          Intended usage:
 *          = PlasmaNoVec: computes eigenvalues only;
 *          = PlasmaVec: computes eigenvalues and eigenvectors.
 *          PlasmaVec is NOT supported.
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in,out] A
 *          On entry, the symmetric (or Hermitian) matrix A.  If uplo
 *          = PlasmaUpper, the leading N-by-N upper triangular part of
 *          A contains the upper triangular part of the matrix A, and
 *          the strictly lower triangular part of A is not referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of
 *          A contains the lower triangular part of the matrix A, and
 *          the strictly upper triangular part of A is not referenced.
 *          On exit, if jobz = PlasmaVec, then if return value = 0, A
 *          contains the orthonormal eigenvectors of the matrix A.
 *          If jobz = PlasmaNoVec, then on exit the lower triangle (if
 *          uplo = PlasmaLower) or the upper triangle (if uplo =
 *          PlasmaUpper) of A, including the diagonal, is destroyed.*
 *
 * @param[out] D
 *          On exit, the diagonal elements of the tridiagonal matrix: 
 *          D(i) = A(i,i).
 *
 * @param[out] E
 *          On exit, he off-diagonal elements of the tridiagonal matrix: 
 *          E(i) = A(i,i+1) if uplo = PlasmaUpper, 
 *          E(i) = A(i+1,i) if uplo = PlasmaLower.
 *
 * @param[out] T
 *          On exit, auxiliary factorization data.
 *
 * @param[out] Q
 *          On exit, if jobz = PlasmaVec and info = 0, the eigenvectors.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if INFO = i, the algorithm failed to converge; i
 *               off-diagonal elements of an intermediate tridiagonal
 *               form did not converge to zero.
 *
 *******************************************************************************
 *
 * @sa PLASMA_dsytrd
 * @sa PLASMA_dsytrd_Tile_Async
 * @sa PLASMA_chetrd_Tile
 * @sa PLASMA_dsytrd_Tile
 * @sa PLASMA_ssytrd_Tile
 * @sa PLASMA_dsyev_Tile
 *
 ******************************************************************************/
int PLASMA_dsytrd_Tile(PLASMA_enum jobz, PLASMA_enum uplo,
                      PLASMA_desc *A, double *D, double *E, PLASMA_desc *T, PLASMA_desc *Q)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dsytrd_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_dsytrd_Tile_Async(jobz, uplo, A, D, E, T, Q, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup double_Tile_Async
 *
 *  PLASMA_dsytrd_Tile_Async - Computes all eigenvalues and,
 *  optionally, eigenvectors of a complex Hermitian matrix A using a
 *  two-stage approach:
 *  First stage: reduction to band tridiagonal form;
 *  Second stage: reduction from band to tridiagonal form.
 *
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
 * @sa PLASMA_dsytrd
 * @sa PLASMA_dsytrd_Tile
 * @sa PLASMA_chetrd_Tile_Async
 * @sa PLASMA_dsytrd_Tile_Async
 * @sa PLASMA_ssytrd_Tile_Async
 *
 ******************************************************************************/
int PLASMA_dsytrd_Tile_Async(PLASMA_enum jobz, PLASMA_enum uplo,
                            PLASMA_desc *A,
                            double *D,
                            double *E,
                            PLASMA_desc *T,
                            PLASMA_desc *Q,
                            PLASMA_sequence *sequence, PLASMA_request *request)
{
    int NB, IB, IBNB, NT;
    PLASMA_desc descA = *A;
    PLASMA_desc descT = *T;

    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dsytrd_Tile_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_dsytrd_Tile_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_dsytrd_Tile_Async", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Set NT & NTRHS */
    NB   = PLASMA_NB;
    IB   = PLASMA_IB;
    IBNB = IB*NB;
    NT   = (descA.ln%NB==0) ? (descA.ln/NB) : (descA.ln/NB+1);

    /* Check descriptors for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dsytrd_Tile_Async", "invalid descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (plasma_desc_check(&descT) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dsytrd_Tile_Async", "invalid descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ( (jobz == PlasmaVec) && (plasma_desc_check(Q) != PLASMA_SUCCESS) ) {
        plasma_error("PLASMA_dsytrd_Tile_Async", "invalid descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    /* Check input arguments */
    if (jobz != PlasmaNoVec && jobz != PlasmaVec) {
        plasma_error("PLASMA_dsytrd_Tile_Async", "illegal value of jobz");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (descA.m != descA.n) {
        plasma_error("PLASMA_dsytrd_Tile_Async", "matrix need to be square");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
   if (descA.nb != descA.mb) {
        plasma_error("PLASMA_dsytrd_Tile_Async", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (jobz == PlasmaVec) {
        plasma_error("PLASMA_dsytrd_Tile_Async", "computing the eigenvectors is not supported in this version");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ( (jobz == PlasmaVec) && (Q->nb != Q->mb) ) {
        plasma_error("PLASMA_dsytrd_Tile_Async", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    /* Reduction to tridiagonal form
     * with a two-stage approach.
     */

    /* Reduction to BAND tridiagonal form
     */
    plasma_dynamic_call_5(plasma_pdsyrbt,
        PLASMA_enum, uplo,
        PLASMA_desc, descA,
        PLASMA_desc, descT,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    /* 
     * Build the Q of the first stage
     */
    /* if (jobz == PlasmaVec){ */
    /*    /\* Initialize Q to Identity *\/ */
    /*    plasma_dynamic_call_6(plasma_pdlaset, */
    /*        PLASMA_enum, PlasmaUpperLower, */
    /*        double, 0.0, */
    /*        double, 1.0, */
    /*        PLASMA_desc, descQ, */
    /*        PLASMA_sequence*, sequence, */
    /*        PLASMA_request*, request); */
    /*    /\* Accumulate the transformations from the first stage*\/ */
    /*    plasma_dynamic_call_6(plasma_pdorgtr, */
    /*        PLASMA_enum, uplo, */
    /*        PLASMA_desc, descA, */
    /*        PLASMA_desc, descQ, */
    /*        PLASMA_desc, descT, */
    /*        PLASMA_sequence*, sequence, */
    /*        PLASMA_request*, request); */
    /* } */

    /* Set the V's to zero before the 2nd stage (bulge chasing) */
    /*
    */
    plasma_dynamic_call_5(plasma_pdlaset2,
        PLASMA_enum, uplo,
        double, 0.0,
        PLASMA_desc, uplo == PlasmaLower ? plasma_desc_submatrix(descA, descA.mb, 0, descA.m-descA.mb, descA.n-descA.nb) 
                                         : plasma_desc_submatrix(descA, 0, descA.nb, descA.m-descA.mb, descA.n-descA.nb),
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    /* Reduction from BAND tridiagonal to the final condensed form
     */
    plasma_dynamic_call_7(plasma_pdsbrdt,
        PLASMA_enum, uplo,
        PLASMA_desc, descA,
        double*, D,
        double*, E,
        PLASMA_desc, descT,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);


    return PLASMA_SUCCESS;
}
