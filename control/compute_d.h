/**
 *
 * @file compute_d.h
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Thu Sep 15 12:09:27 2011
 *
 **/

/***************************************************************************//**
 *  Macro for matrix conversion / Lapack interface
 **/
#define plasma_ddesc_alloc( descA, mb, nb, lm, ln, i, j, m, n, free)   \
    descA = plasma_desc_init(                                          \
        PlasmaRealDouble, (mb), (nb), ((mb)*(nb)),                  \
        (m), (n), (i), (j), (m), (n));                                 \
    if ( plasma_desc_mat_alloc( &(descA) ) ) {                         \
        plasma_error( __func__, "plasma_shared_alloc() failed");       \
        {free;};                                                       \
        return PLASMA_ERR_OUT_OF_RESOURCES;                            \
    }

#define plasma_dooplap2tile( descA, A, mb, nb, lm, ln, i, j, m, n, free) \
    descA = plasma_desc_init(                                           \
        PlasmaRealDouble, (mb), (nb), ((mb)*(nb)),                   \
        (lm), (ln), (i), (j), (m), (n));                                  \
    if ( plasma_desc_mat_alloc( &(descA) ) ) {                          \
        plasma_error( __func__, "plasma_shared_alloc() failed");        \
        {free;};                                                        \
        return PLASMA_ERR_OUT_OF_RESOURCES;                             \
    }                                                                   \
    plasma_parallel_call_5(                                             \
        plasma_pdlapack_to_tile,                                        \
        double*, (A),                                       \
        int,                 (lm),                                      \
        PLASMA_desc,         (descA),                                   \
        PLASMA_sequence*,    sequence,                                  \
        PLASMA_request*,     &request);

#define plasma_diplap2tile( descA, A, mb, nb, lm, ln, i, j, m, n)     \
    descA = plasma_desc_init(                                         \
        PlasmaRealDouble, (mb), (nb), ((mb)*(nb)),                 \
        (lm), (ln), (i), (j), (m), (n));                              \
    descA.mat = A;                                                    \
    PLASMA_dgecfi_Async((lm), (ln), (A), PlasmaCM, (mb), (nb),        \
                        PlasmaCCRB, (mb), (nb), sequence, &request);



#define plasma_dooptile2lap( descA, A, mb, nb, lm, ln)  \
    plasma_parallel_call_5(plasma_pdtile_to_lapack,     \
        PLASMA_desc,         (descA),                   \
        double*, (A),                       \
        int,                 (lm),                      \
        PLASMA_sequence*,    sequence,                  \
        PLASMA_request*,     &request);

#define plasma_diptile2lap( descA, A, mb, nb, lm, ln)                   \
    PLASMA_dgecfi_Async((lm), (ln), (A), PlasmaCCRB, (mb), (nb),        \
                        PlasmaCM, (mb), (nb), sequence, &request);

/***************************************************************************//**
 *  Declarations of parallel functions (static scheduling) - alphabetical order
 **/
void plasma_pdaxpy  (plasma_context_t *plasma);
void plasma_pdgelqf (plasma_context_t *plasma);
void plasma_pdgemm  (plasma_context_t *plasma);
void plasma_pdgeqrf (plasma_context_t *plasma);
void plasma_pdgerbb (plasma_context_t *plasma);
void plasma_pdgetmi2(plasma_context_t *plasma);
void plasma_pdgetrf_incpiv(plasma_context_t *plasma);
#ifdef COMPLEX
void plasma_pdsymm  (plasma_context_t *plasma);
void plasma_pdsyrk  (plasma_context_t *plasma);
void plasma_pdsyr2k (plasma_context_t *plasma);
#endif
void plasma_pdlacpy (plasma_context_t *plasma);
void plasma_pdlag2s (plasma_context_t *plasma);
void plasma_pdlange (plasma_context_t *plasma);
#ifdef COMPLEX
void plasma_pdlansy (plasma_context_t *plasma);
#endif
void plasma_pdlansy (plasma_context_t *plasma);
void plasma_pdpack  (plasma_context_t *plasma);
void plasma_pdplgsy (plasma_context_t *plasma);
void plasma_pdplgsy (plasma_context_t *plasma);
void plasma_pdplrnt (plasma_context_t *plasma);
void plasma_pdpotrf (plasma_context_t *plasma);
void plasma_pdshift (plasma_context_t *plasma);
void plasma_pdsymm  (plasma_context_t *plasma);
void plasma_pdsyrk  (plasma_context_t *plasma);
void plasma_pdsyr2k (plasma_context_t *plasma);
void plasma_pdtrmm  (plasma_context_t *plasma);
void plasma_pdtrsm  (plasma_context_t *plasma);
void plasma_pdtrsmpl(plasma_context_t *plasma);
void plasma_pdorglq (plasma_context_t *plasma);
void plasma_pdorgqr (plasma_context_t *plasma);
void plasma_pdorgqrrh(plasma_context_t *plasma);
void plasma_pdormlq (plasma_context_t *plasma);
void plasma_pdormqr (plasma_context_t *plasma);
void plasma_pdunpack(plasma_context_t *plasma);

/***************************************************************************//**
 *  Declarations of internal sequential functions 
 **/
int plasma_dshift(plasma_context_t *plasma, int m, int n, double *A,
                  int nprob, int me, int ne, int L,
                  PLASMA_sequence *sequence, PLASMA_request *request);

/***************************************************************************//**
 *  Declarations of parallel functions (dynamic scheduling) - alphabetical order
 **/
void plasma_pdaxpy_quark(double alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdbarrier_tl2pnl_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdbarrier_pnl2tl_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgelqf_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgelqfrh_quark(PLASMA_desc A, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgemm_quark(PLASMA_enum transA, PLASMA_enum transB, double alpha, PLASMA_desc A, PLASMA_desc B, double beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgeqrf_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgeqrfrh_quark(PLASMA_desc A, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgerbh_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgerbb_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgerbbrh_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgetmi2_quark(PLASMA_enum idep, PLASMA_enum odep, PLASMA_enum storev, int m, int n, int mb, int nb, double *A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgetrf_incpiv_quark(PLASMA_desc A, PLASMA_desc L, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgetrf_reclap_quark(PLASMA_desc A, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgetrf_rectil_quark(PLASMA_desc A, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
#ifdef COMPLEX
void plasma_pdsygst_quark(PLASMA_enum itype, PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdsymm_quark(PLASMA_enum side, PLASMA_enum uplo, double alpha, PLASMA_desc A, PLASMA_desc B, double beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdsyrk_quark(PLASMA_enum uplo, PLASMA_enum trans, double alpha, PLASMA_desc A, double beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdsyr2k_quark(PLASMA_enum uplo, PLASMA_enum trans, double alpha, PLASMA_desc A, PLASMA_desc B, double beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdsyrbt_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgbrdb_quark(PLASMA_enum uplo, PLASMA_desc A, double *D, double *E, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdsbrdt_quark(PLASMA_enum uplo, PLASMA_desc A, double *D, double *E, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
void plasma_pdlacpy_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdlag2s_quark(PLASMA_desc A, PLASMA_desc SB, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdlange_quark(PLASMA_enum norm, PLASMA_desc A, double *work, double *result, PLASMA_sequence *sequence, PLASMA_request *request);
#ifdef COMPLEX
void plasma_pdlansy_quark(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc A, double *work, double *result, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
void plasma_pdlansy_quark(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc A, double *work, double *result, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdlaset_quark( PLASMA_enum uplo, double alpha, double beta, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdlaset2_quark(PLASMA_enum uplo, double alpha,                          PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdlaswp_quark(PLASMA_desc B, int *IPIV, int inc, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdlauum_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdplgsy_quark(double bump, PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdplgsy_quark(double bump, PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pdplrnt_quark(PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pdpotrf_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdshift_quark(int, int, int, double *, int *, int, int, PLASMA_sequence*, PLASMA_request*);
void plasma_pdsymm_quark(PLASMA_enum side, PLASMA_enum uplo, double alpha, PLASMA_desc A, PLASMA_desc B, double beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdsyrk_quark(PLASMA_enum uplo, PLASMA_enum trans, double alpha, PLASMA_desc A, double beta,  PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdsyr2k_quark(PLASMA_enum uplo, PLASMA_enum trans, double alpha, PLASMA_desc A, PLASMA_desc B, double beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdtrmm_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, double alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdtrsm_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, double alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdtrsmpl_quark(PLASMA_desc A, PLASMA_desc B, PLASMA_desc L, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdtrtri_quark(PLASMA_enum uplo, PLASMA_enum diag, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdorgbr_quark(PLASMA_enum side, PLASMA_desc A, PLASMA_desc O, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdorgbrrh_quark(PLASMA_enum side, PLASMA_desc A, PLASMA_desc O, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdorgqr_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdorgqrrh_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdorglq_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdorglqrh_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdorgtr_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdormqr_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdormqrrh_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdormlq_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdormlqrh_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
