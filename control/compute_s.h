/**
 *
 * @file compute_s.h
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Thu Sep 15 12:09:27 2011
 *
 **/

/***************************************************************************//**
 *  Macro for matrix conversion / Lapack interface
 **/
#define plasma_sdesc_alloc( descA, mb, nb, lm, ln, i, j, m, n, free)   \
    descA = plasma_desc_init(                                          \
        PlasmaRealFloat, (mb), (nb), ((mb)*(nb)),                  \
        (m), (n), (i), (j), (m), (n));                                 \
    if ( plasma_desc_mat_alloc( &(descA) ) ) {                         \
        plasma_error( __func__, "plasma_shared_alloc() failed");       \
        {free;};                                                       \
        return PLASMA_ERR_OUT_OF_RESOURCES;                            \
    }

#define plasma_sooplap2tile( descA, A, mb, nb, lm, ln, i, j, m, n, free) \
    descA = plasma_desc_init(                                           \
        PlasmaRealFloat, (mb), (nb), ((mb)*(nb)),                   \
        (lm), (ln), (i), (j), (m), (n));                                  \
    if ( plasma_desc_mat_alloc( &(descA) ) ) {                          \
        plasma_error( __func__, "plasma_shared_alloc() failed");        \
        {free;};                                                        \
        return PLASMA_ERR_OUT_OF_RESOURCES;                             \
    }                                                                   \
    plasma_parallel_call_5(                                             \
        plasma_pslapack_to_tile,                                        \
        float*, (A),                                       \
        int,                 (lm),                                      \
        PLASMA_desc,         (descA),                                   \
        PLASMA_sequence*,    sequence,                                  \
        PLASMA_request*,     &request);

#define plasma_siplap2tile( descA, A, mb, nb, lm, ln, i, j, m, n)     \
    descA = plasma_desc_init(                                         \
        PlasmaRealFloat, (mb), (nb), ((mb)*(nb)),                 \
        (lm), (ln), (i), (j), (m), (n));                              \
    descA.mat = A;                                                    \
    PLASMA_sgecfi_Async((lm), (ln), (A), PlasmaCM, (mb), (nb),        \
                        PlasmaCCRB, (mb), (nb), sequence, &request);



#define plasma_sooptile2lap( descA, A, mb, nb, lm, ln)  \
    plasma_parallel_call_5(plasma_pstile_to_lapack,     \
        PLASMA_desc,         (descA),                   \
        float*, (A),                       \
        int,                 (lm),                      \
        PLASMA_sequence*,    sequence,                  \
        PLASMA_request*,     &request);

#define plasma_siptile2lap( descA, A, mb, nb, lm, ln)                   \
    PLASMA_sgecfi_Async((lm), (ln), (A), PlasmaCCRB, (mb), (nb),        \
                        PlasmaCM, (mb), (nb), sequence, &request);

/***************************************************************************//**
 *  Declarations of parallel functions (static scheduling) - alphabetical order
 **/
void plasma_psaxpy  (plasma_context_t *plasma);
void plasma_psgelqf (plasma_context_t *plasma);
void plasma_psgemm  (plasma_context_t *plasma);
void plasma_psgeqrf (plasma_context_t *plasma);
void plasma_psgerbb (plasma_context_t *plasma);
void plasma_psgetmi2(plasma_context_t *plasma);
void plasma_psgetrf_incpiv(plasma_context_t *plasma);
#ifdef COMPLEX
void plasma_pssymm  (plasma_context_t *plasma);
void plasma_pssyrk  (plasma_context_t *plasma);
void plasma_pssyr2k (plasma_context_t *plasma);
#endif
void plasma_pslacpy (plasma_context_t *plasma);
void plasma_pslag2d (plasma_context_t *plasma);
void plasma_pslange (plasma_context_t *plasma);
#ifdef COMPLEX
void plasma_pslansy (plasma_context_t *plasma);
#endif
void plasma_pslansy (plasma_context_t *plasma);
void plasma_pspack  (plasma_context_t *plasma);
void plasma_psplgsy (plasma_context_t *plasma);
void plasma_psplgsy (plasma_context_t *plasma);
void plasma_psplrnt (plasma_context_t *plasma);
void plasma_pspotrf (plasma_context_t *plasma);
void plasma_psshift (plasma_context_t *plasma);
void plasma_pssymm  (plasma_context_t *plasma);
void plasma_pssyrk  (plasma_context_t *plasma);
void plasma_pssyr2k (plasma_context_t *plasma);
void plasma_pstrmm  (plasma_context_t *plasma);
void plasma_pstrsm  (plasma_context_t *plasma);
void plasma_pstrsmpl(plasma_context_t *plasma);
void plasma_psorglq (plasma_context_t *plasma);
void plasma_psorgqr (plasma_context_t *plasma);
void plasma_psorgqrrh(plasma_context_t *plasma);
void plasma_psormlq (plasma_context_t *plasma);
void plasma_psormqr (plasma_context_t *plasma);
void plasma_psunpack(plasma_context_t *plasma);

/***************************************************************************//**
 *  Declarations of internal sequential functions 
 **/
int plasma_sshift(plasma_context_t *plasma, int m, int n, float *A,
                  int nprob, int me, int ne, int L,
                  PLASMA_sequence *sequence, PLASMA_request *request);

/***************************************************************************//**
 *  Declarations of parallel functions (dynamic scheduling) - alphabetical order
 **/
void plasma_psaxpy_quark(float alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psbarrier_tl2pnl_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psbarrier_pnl2tl_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgelqf_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgelqfrh_quark(PLASMA_desc A, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgemm_quark(PLASMA_enum transA, PLASMA_enum transB, float alpha, PLASMA_desc A, PLASMA_desc B, float beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgeqrf_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgeqrfrh_quark(PLASMA_desc A, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgerbh_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgerbb_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgerbbrh_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgetmi2_quark(PLASMA_enum idep, PLASMA_enum odep, PLASMA_enum storev, int m, int n, int mb, int nb, float *A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgetrf_incpiv_quark(PLASMA_desc A, PLASMA_desc L, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgetrf_reclap_quark(PLASMA_desc A, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgetrf_rectil_quark(PLASMA_desc A, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
#ifdef COMPLEX
void plasma_pssygst_quark(PLASMA_enum itype, PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pssymm_quark(PLASMA_enum side, PLASMA_enum uplo, float alpha, PLASMA_desc A, PLASMA_desc B, float beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pssyrk_quark(PLASMA_enum uplo, PLASMA_enum trans, float alpha, PLASMA_desc A, float beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pssyr2k_quark(PLASMA_enum uplo, PLASMA_enum trans, float alpha, PLASMA_desc A, PLASMA_desc B, float beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pssyrbt_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgbrdb_quark(PLASMA_enum uplo, PLASMA_desc A, float *D, float *E, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pssbrdt_quark(PLASMA_enum uplo, PLASMA_desc A, float *D, float *E, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
void plasma_pslacpy_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pslag2d_quark(PLASMA_desc A, PLASMA_desc SB, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pslange_quark(PLASMA_enum norm, PLASMA_desc A, float *work, float *result, PLASMA_sequence *sequence, PLASMA_request *request);
#ifdef COMPLEX
void plasma_pslansy_quark(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc A, float *work, float *result, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
void plasma_pslansy_quark(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc A, float *work, float *result, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pslaset_quark( PLASMA_enum uplo, float alpha, float beta, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pslaset2_quark(PLASMA_enum uplo, float alpha,                          PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pslaswp_quark(PLASMA_desc B, int *IPIV, int inc, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pslauum_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psplgsy_quark(float bump, PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psplgsy_quark(float bump, PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_psplrnt_quark(PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pspotrf_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psshift_quark(int, int, int, float *, int *, int, int, PLASMA_sequence*, PLASMA_request*);
void plasma_pssymm_quark(PLASMA_enum side, PLASMA_enum uplo, float alpha, PLASMA_desc A, PLASMA_desc B, float beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pssyrk_quark(PLASMA_enum uplo, PLASMA_enum trans, float alpha, PLASMA_desc A, float beta,  PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pssyr2k_quark(PLASMA_enum uplo, PLASMA_enum trans, float alpha, PLASMA_desc A, PLASMA_desc B, float beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pstrmm_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, float alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pstrsm_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, float alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pstrsmpl_quark(PLASMA_desc A, PLASMA_desc B, PLASMA_desc L, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pstrtri_quark(PLASMA_enum uplo, PLASMA_enum diag, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psorgbr_quark(PLASMA_enum side, PLASMA_desc A, PLASMA_desc O, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psorgbrrh_quark(PLASMA_enum side, PLASMA_desc A, PLASMA_desc O, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psorgqr_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psorgqrrh_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psorglq_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psorglqrh_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psorgtr_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psormqr_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psormqrrh_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psormlq_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psormlqrh_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
