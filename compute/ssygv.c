/**
 * @file ssygv.c
 *
 *  PLASMA computational routines
 *  Release Date: November, 15th 2009
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @generated s Thu Sep 15 12:09:18 2011
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup float
 *
 *  PLASMA_ssygv - Computes all eigenvalues and, optionally,
 *  eigenvectors of a complex generalized Hermitian-definite
 *  eigenproblem of the form: A*x=(lambda)*B*x, A*Bx=(lambda)*x, or
 *  B*A*x=(lambda)*x.
 *  Here A and B are assumed to be Hermitian and B is also positive
 *  definite.
 *  Note: Only PlasmaNoVec supported!
 *
 *******************************************************************************
 *
 * @param[in] PlasmaItype
 *          Intended usage:
 *          = 1: A*x=(lambda)*B*x
 *          = 2: A*Bx=(lambda)*x 
 *          = 3: B*A*x=(lambda)*x 
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
 *          = PlasmaUpper: Upper triangle of A and B are stored;
 *          = PlasmaLower: Lower triangle of A and B are stored.
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
 *          On exit, if jobz = PlasmaVec, then if return value = 0, A
 *          contains the matrix Z of eigenvectors.
 *          The eigenvectors are normalized as follows:
 *          if ITYPE = 1 or 2, Z**T*B*Z = I;
 *          if ITYPE = 3,      Z**T*inv(B)*Z = I.
 *          If jobz = PlasmaNoVec, then on exit the lower triangle (if
 *          uplo = PlasmaLower) or the upper triangle (if uplo =
 *          PlasmaUpper) of A, including the diagonal, is destroyed.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in,out] B
 *          On entry, the symmetric (or Hermitian) positive definite
 *          matrix B.
 *          If uplo = PlasmaUpper, the leading N-by-N upper triangular
 *          part of B contains the upper triangular part of the matrix
 *          B, and the strictly lower triangular part of B is not
 *          referenced.
 *          If uplo = PlasmaLower, the leading N-by-N lower triangular
 *          part of B contains the lower triangular part of the matrix
 *          B, and the strictly upper triangular part of B is not
 *          referenced.
 *          On exit, if return value <= N, the part of B containing
 *          the matrix is overwritten by the triangular factor U or L
 *          from the Cholesky factorization B = U**T*U or B = L*L**T.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDA >= max(1,N).
 *
 * @param[out] W
 *          On exit, if info = 0, the eigenvalues.
 *
 * @param[in, out] descT
 *          On entry, descriptor as return by PLASMA_Alloc_Workspace_ssygv
 *          On exit, contains auxiliary factorization data.
 *
 * @param[out] Q
 *          On exit, if jobz = PlasmaVec and info = 0, the eigenvectors.
 *
 * @param[in] LDQ
 *          The leading dimension of Q.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval <=N if INFO = i, plasma_ssygv failed to converge; i
 *               off-diagonal elements of an intermediate tridiagonal
 *               form did not converge to zero.
 *          \retval >N if INFO = N + i, for 1 <= i <= N, then the leading
 *                     minor of order i of B is not positive definite.
 *                     The factorization of B could not be completed and
 *                     no eigenvalues or eigenvectors were computed.
 *
 *******************************************************************************
 *
 * @sa PLASMA_ssygv_Tile
 * @sa PLASMA_ssygv_Tile_Async
 * @sa PLASMA_chegv
 * @sa PLASMA_dsygv
 * @sa PLASMA_ssygv
 *
 ******************************************************************************/
int PLASMA_ssygv(PLASMA_enum itype, PLASMA_enum jobz, PLASMA_enum uplo, int N,
                 float *A, int LDA,
                 float *B, int LDB,
                 float *W,
                 PLASMA_desc *descT,
                 float *Q, int LDQ)
{
    int NB, IB, IBNB, NT;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA, descB, descQ;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_error("PLASMA_ssygv", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    /* Tune NB & IB depending on N; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_SSYGV, N, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_ssygv", "plasma_tune() failed");
        return status;
    }

    /* Set NT */
    NB   = PLASMA_NB;
    IB   = PLASMA_IB;
    IBNB = IB*NB;
    NT   = (N%NB==0) ? (N/NB) : (N/NB+1);

    /* Check input arguments */
    if (itype != 1 && itype != 2 && itype != 3) {
        plasma_error("PLASMA_ssygv", "Illegal value of itype");
        return -1;
    }
    if (jobz != PlasmaNoVec && jobz != PlasmaVec) {
        plasma_error("PLASMA_ssygv", "illegal value of jobz");
        return -2;
    }
    if (uplo != PlasmaLower && uplo!= PlasmaUpper) {
        plasma_error("PLASMA_ssygv", "only PlasmaLower supported");
        return -3;
    }
    if (N < 0) {
        plasma_error("PLASMA_ssygv", "illegal value of N");
        return -4;
    }
    if (LDA < max(1, N)) {
        plasma_error("PLASMA_ssygv", "illegal value of LDA");
        return -6;
    }
    if (LDB < max(1, N)) {
        plasma_error("PLASMA_ssygv", "illegal value of LDB");
        return -8;
    }
    if ( (plasma_desc_check(descT) != PLASMA_SUCCESS) || 
         ( descT->m != NT*IB ) || (descT->n != NT*NB) ) {
        plasma_error("PLASMA_ssygv", "invalid T descriptor");
        return -10;
    }
    if (LDQ < max(1, N)) {
        plasma_error("PLASMA_ssygv", "illegal value of LDQ");
        return -12;
    }
    
    /* Quick return */
    if (N == 0)
        return PLASMA_SUCCESS;

    if (jobz == PlasmaVec) {
        plasma_error("PLASMA_ssygv", "computing the eigenvectors is not supported in this version");
        return -1;
    }

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_sooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N, 
                             plasma_desc_mat_free(&(descA)) );
        plasma_sooplap2tile( descB, B, NB, NB, LDB, N, 0, 0, N, N, 
                             plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)) );
        if (jobz == PlasmaVec) {
            /* No need for conversion, it's just output */
            plasma_sdesc_alloc( descQ, NB, NB, LDQ, N, 0, 0, N, N, 
                                plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)); plasma_desc_mat_free(&(descQ)) );
        }
    } else {
        plasma_siplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N );
        plasma_siplap2tile( descB, B, NB, NB, LDB, N, 0, 0, N, N );
        if (jobz == PlasmaVec) {
            /* No need for conversion, it's just output */
            descQ = plasma_desc_init(
                PlasmaRealFloat, NB, NB, NB*NB,
                LDQ, N, 0, 0, N, N);
            descQ.mat = Q;
        }
    }

    /* Call the tile interface */
    PLASMA_ssygv_Tile_Async(itype, PlasmaNoVec, uplo, 
                            &descA, &descB, W, 
                            descT, &descQ, 
                            sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_sooptile2lap( descA, A, NB, NB, LDA, N );
        plasma_sooptile2lap( descB, B, NB, NB, LDB, N );
        if (jobz == PlasmaVec) {
           plasma_sooptile2lap( descQ, Q, NB, NB, LDQ, N );
        }
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
        plasma_desc_mat_free(&descB);
        if (jobz == PlasmaVec)
           plasma_desc_mat_free(&descQ);
    } else {
        plasma_siptile2lap( descA, A, NB, NB, LDA, N );
        plasma_siptile2lap( descB, B, NB, NB, LDB, N );
        if (jobz == PlasmaVec)
           plasma_siptile2lap( descQ, Q, NB, NB, LDQ, N );
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
 *  PLASMA_ssygv_Tile - Computes all eigenvalues and, optionally,
 *  eigenvectors of a complex generalized Hermitian-definite
 *  eigenproblem of the form: A*x=(lambda)*B*x, A*Bx=(lambda)*x, or
 *  B*A*x=(lambda)*x.
 *  Here A and B are assumed to be Hermitian and B is also
 *  positive definite.
 *
 *  Tile equivalent of PLASMA_ssygv().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] PlasmaItype
 *          Intended usage:
 *          = 1: A*x=(lambda)*B*x
 *          = 2: A*Bx=(lambda)*x 
 *          = 3: B*A*x=(lambda)*x 
 *
 * @param[in] jobz
 *          Intended usage:
 *          = PlasmaNoVec: computes eigenvalues only;
 *          = PlasmaVec: computes eigenvalues and eigenvectors.
 *          Currently only PlasmaVec is NOT supported.
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or
 *          lower triangular:
 *          = PlasmaUpper: Upper triangle of A and B are stored;
 *          = PlasmaLower: Lower triangle of A and B are stored.
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
 *          On exit, if jobz = PlasmaVec, then if return value = 0, A
 *          contains the matrix Z of eigenvectors.
 *          The eigenvectors are normalized as follows:
 *          if ITYPE = 1 or 2, Z**T*B*Z = I;
 *          if ITYPE = 3,      Z**T*inv(B)*Z = I.
 *          If jobz = PlasmaNoVec, then on exit the lower triangle (if
 *          uplo = PlasmaLower) or the upper triangle (if uplo =
 *          PlasmaUpper) of A, including the diagonal, is destroyed.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in,out] B
 *          On entry, the symmetric (or Hermitian) positive definite
 *          matrix B.
 *          If uplo = PlasmaUpper, the leading N-by-N upper triangular
 *          part of B contains the upper triangular part of the matrix
 *          B, and the strictly lower triangular part of B is not
 *          referenced.
 *          If uplo = PlasmaLower, the leading N-by-N lower triangular
 *          part of B contains the lower triangular part of the matrix
 *          B, and the strictly upper triangular part of B is not
 *          referenced.
 *          On exit, if return value <= N, the part of B containing
 *          the matrix is overwritten by the triangular factor U or L
 *          from the Cholesky factorization B = U**T*U or B = L*L**T.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDA >= max(1,N).
 *
 * @param[in,out] T
 *          On entry, descriptor as return by
 *          PLASMA_Alloc_Workspace_ssygv
 *          On exit, contains auxiliary factorization data.
 *
 * @param[out] W
 *          On exit, if info = 0, the eigenvalues.
 *
 * @param[out] Q
 *          On exit, if jobz = PlasmaVec and info = 0, the eigenvectors.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval <=N if INFO = i, plasma_ssygv failed to converge; i
 *               off-diagonal elements of an intermediate tridiagonal
 *               form did not converge to zero.
 *          \retval >N if INFO = N + i, for 1 <= i <= N, then the leading
 *                     minor of order i of B is not positive definite.
 *                     The factorization of B could not be completed and
 *                     no eigenvalues or eigenvectors were computed.
 *
 *******************************************************************************
 *
 * @sa PLASMA_ssygv
 * @sa PLASMA_ssygv_Tile_Async
 * @sa PLASMA_chegv_Tile
 * @sa PLASMA_dsygv_Tile
 * @sa PLASMA_ssygv_Tile
 *
 ******************************************************************************/
int PLASMA_ssygv_Tile(PLASMA_enum itype, PLASMA_enum jobz, PLASMA_enum uplo, 
                      PLASMA_desc *A, PLASMA_desc *B, float *W, 
                      PLASMA_desc *T, PLASMA_desc *Q)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_ssygv_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_ssygv_Tile_Async(itype, jobz, uplo, A, B, W, T, Q, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup float_Tile_Async
 *
 *  PLASMA_ssygv_Tile - Computes all eigenvalues and, optionally,
 *  eigenvectors of a complex generalized Hermitian-definite
 *  eigenproblem of the form: A*x=(lambda)*B*x, A*Bx=(lambda)*x, or
 *  B*A*x=(lambda)*x.
 *  Here A and B are assumed to be Hermitian and B is also
 *  positive definite.
 * 
 *  Non-blocking equivalent of PLASMA_ssygv_Tile().
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
 * @sa PLASMA_ssygv
 * @sa PLASMA_ssygv_Tile
 * @sa PLASMA_chegv_Tile_Async
 * @sa PLASMA_dsygv_Tile_Async
 * @sa PLASMA_ssygv_Tile_Async
 *
 ******************************************************************************/
int PLASMA_ssygv_Tile_Async(PLASMA_enum itype, PLASMA_enum jobz, PLASMA_enum uplo, 
                            PLASMA_desc *A,
                            PLASMA_desc *B,
                            float *W,
                            PLASMA_desc *T,
                            PLASMA_desc *Q,
                            PLASMA_sequence *sequence, PLASMA_request *request)
{
    int NB, IB, IBNB, NT;
    int status;
    PLASMA_desc descA = *A;
    PLASMA_desc descB = *B;
    PLASMA_desc descT = *T;
    float *E;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_ssygv_Tile_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_ssygv_Tile_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_ssygv_Tile_Async", "NULL request");
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
        plasma_error("PLASMA_ssygv_Tile_Async", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (plasma_desc_check(&descB) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_ssygv_Tile_Async", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (plasma_desc_check(&descT) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_ssygv_Tile_Async", "invalid descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (jobz == PlasmaVec){
       if (plasma_desc_check(Q) != PLASMA_SUCCESS) {
           plasma_error("PLASMA_ssygv_Tile_Async", "invalid descriptor");
           return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
       }
    }
    /* Check input arguments */
    if (itype != 1 && itype != 2 && itype != 3) {
        plasma_error("PLASMA_ssygv_Tile_Async", "Illegal value of itype");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (jobz != PlasmaNoVec && jobz != PlasmaVec) {
        plasma_error("PLASMA_ssygv_Tile_Async", "illegal value of jobz");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (uplo != PlasmaLower && uplo != PlasmaUpper) {
        plasma_error("PLASMA_ssyev_Tile_Async", "illegal value of uplo");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_ssygv_Tile_Async", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (jobz == PlasmaVec) {
        plasma_error("PLASMA_ssygv_Tile_Async", "computing the eigenvectors is not supported in this version");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ( (jobz == PlasmaVec) && (Q->nb != Q->mb) ) {
        plasma_error("PLASMA_ssygv_Tile_Async", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    E = (float *)plasma_shared_alloc(plasma, descA.n-1, PlasmaRealDouble);

    /* Currently NOT equivalent to LAPACK's
     */
    plasma_parallel_call_4(plasma_pspotrf,
        PLASMA_enum, uplo,
        PLASMA_desc, descB,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    status = sequence->status;
    if (status != 0){
       status = descA.ln + status;
       return status;
    }

    /* 
     * Transform problem to standard eigenvalue problem and solve
     */
    plasma_dynamic_call_6(plasma_pssygst,
        PLASMA_enum, itype,
        PLASMA_enum, uplo,
        PLASMA_desc, descA,
        PLASMA_desc, descB,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    /* Reduction to tridiagonal form
     * with a two-stage approach.
     */

    /* 
     *Reduction to BAND tridiagonal form
     */
    plasma_dynamic_call_5(plasma_pssyrbt,
        PLASMA_enum, uplo,
        PLASMA_desc, descA,
        PLASMA_desc, descT,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    /*
     * Build the Q of the first stage
     */
    /* if (jobz == PlasmaVec){ */
    /*     /\* Initialize Q to Identity *\/ */
    /*     plasma_dynamic_call_6(plasma_pslaset, */
    /*                           PLASMA_enum, PlasmaUpperLower, */
    /*                           float, 0.0, */
    /*                           float, 1.0, */
    /*                           PLASMA_desc, descQ, */
    /*                           PLASMA_sequence*, sequence, */
    /*                           PLASMA_request*, request); */

    /*     /\* Accumulate the transformations from the first stage *\/ */
    /*     plasma_dynamic_call_6(plasma_psorgtr, */
    /*                           PLASMA_enum, uplo, */
    /*                           PLASMA_desc, descA, */
    /*                           PLASMA_desc, descQ, */
    /*                           PLASMA_desc, descT, */
    /*                           PLASMA_sequence*, sequence, */
    /*                           PLASMA_request*, request); */
    /* } */
    
    /* Set the V's to zero before the 2nd stage (bulge chasing) */
    plasma_dynamic_call_5(
        plasma_pslaset2,
        PLASMA_enum, uplo,
        float, 0.0,
        PLASMA_desc, uplo == PlasmaLower ? plasma_desc_submatrix(descA, descA.mb, 0, descA.m-descA.mb, descA.n-descA.nb) : 
                                           plasma_desc_submatrix(descA, 0, descA.nb, descA.m-descA.mb, descA.n-descA.nb),
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);
    /* 
     * Reduction from BAND tridiagonal to the final condensed form
     */
    plasma_dynamic_call_7(plasma_pssbrdt,
        PLASMA_enum, uplo,
        PLASMA_desc, descA,
        float*, W,
        float*, E,
        PLASMA_desc, descT,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    /* For eigenvalues only, call DSTERF.  
     * For eigenvectors, first call SORGTR to generate the unitary matrix, 
     * then call ZSTEQR.
     */

    if (jobz == PlasmaNoVec)
        status = LAPACKE_ssterf(descA.n, W, E);
    else {
        status = LAPACKE_ssterf(descA.n, W, E);
        /* Accumulate the transformations from the second stage */
        /*
        plasma_dynamic_call_6(plasma_psorgtr,
            PLASMA_enum, uplo,
            PLASMA_desc, descA,
            PLASMA_desc, descQ,
            PLASMA_desc, descT,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);
        LAPACKE_ssteqr(jobz, descA.n, W, E, descQ.mat, descQ.lm);
        */
    }

    /* If matrix was scaled, then rescale eigenvalues appropriately.
     */

    plasma_shared_free(plasma, E);

    return PLASMA_SUCCESS;
}
