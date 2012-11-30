/**
 *
 * @file compute_c.h
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated c Thu Sep 15 12:09:27 2011
 *
 **/

/***************************************************************************//**
 *  Macro for matrix conversion / Lapack interface
 **/
#define plasma_cdesc_alloc( descA, mb, nb, lm, ln, i, j, m, n, free)   \
    descA = plasma_desc_init(                                          \
        PlasmaComplexFloat, (mb), (nb), ((mb)*(nb)),                  \
        (m), (n), (i), (j), (m), (n));                                 \
    if ( plasma_desc_mat_alloc( &(descA) ) ) {                         \
        plasma_error( __func__, "plasma_shared_alloc() failed");       \
        {free;};                                                       \
        return PLASMA_ERR_OUT_OF_RESOURCES;                            \
    }

#define plasma_cooplap2tile( descA, A, mb, nb, lm, ln, i, j, m, n, free) \
    descA = plasma_desc_init(                                           \
        PlasmaComplexFloat, (mb), (nb), ((mb)*(nb)),                   \
        (lm), (ln), (i), (j), (m), (n));                                  \
    if ( plasma_desc_mat_alloc( &(descA) ) ) {                          \
        plasma_error( __func__, "plasma_shared_alloc() failed");        \
        {free;};                                                        \
        return PLASMA_ERR_OUT_OF_RESOURCES;                             \
    }                                                                   \
    plasma_parallel_call_5(                                             \
        plasma_pclapack_to_tile,                                        \
        PLASMA_Complex32_t*, (A),                                       \
        int,                 (lm),                                      \
        PLASMA_desc,         (descA),                                   \
        PLASMA_sequence*,    sequence,                                  \
        PLASMA_request*,     &request);

#define plasma_ciplap2tile( descA, A, mb, nb, lm, ln, i, j, m, n)     \
    descA = plasma_desc_init(                                         \
        PlasmaComplexFloat, (mb), (nb), ((mb)*(nb)),                 \
        (lm), (ln), (i), (j), (m), (n));                              \
    descA.mat = A;                                                    \
    PLASMA_cgecfi_Async((lm), (ln), (A), PlasmaCM, (mb), (nb),        \
                        PlasmaCCRB, (mb), (nb), sequence, &request);



#define plasma_cooptile2lap( descA, A, mb, nb, lm, ln)  \
    plasma_parallel_call_5(plasma_pctile_to_lapack,     \
        PLASMA_desc,         (descA),                   \
        PLASMA_Complex32_t*, (A),                       \
        int,                 (lm),                      \
        PLASMA_sequence*,    sequence,                  \
        PLASMA_request*,     &request);

#define plasma_ciptile2lap( descA, A, mb, nb, lm, ln)                   \
    PLASMA_cgecfi_Async((lm), (ln), (A), PlasmaCCRB, (mb), (nb),        \
                        PlasmaCM, (mb), (nb), sequence, &request);

/***************************************************************************//**
 *  Declarations of parallel functions (static scheduling) - alphabetical order
 **/
void plasma_pcaxpy  (plasma_context_t *plasma);
void plasma_pcgelqf (plasma_context_t *plasma);
void plasma_pcgemm  (plasma_context_t *plasma);
void plasma_pcgeqrf (plasma_context_t *plasma);
void plasma_pcgerbb (plasma_context_t *plasma);
void plasma_pcgetmi2(plasma_context_t *plasma);
void plasma_pcgetrf_incpiv(plasma_context_t *plasma);
#ifdef COMPLEX
void plasma_pchemm  (plasma_context_t *plasma);
void plasma_pcherk  (plasma_context_t *plasma);
void plasma_pcher2k (plasma_context_t *plasma);
#endif
void plasma_pclacpy (plasma_context_t *plasma);
void plasma_pclag2z (plasma_context_t *plasma);
void plasma_pclange (plasma_context_t *plasma);
#ifdef COMPLEX
void plasma_pclanhe (plasma_context_t *plasma);
#endif
void plasma_pclansy (plasma_context_t *plasma);
void plasma_pcpack  (plasma_context_t *plasma);
void plasma_pcplghe (plasma_context_t *plasma);
void plasma_pcplgsy (plasma_context_t *plasma);
void plasma_pcplrnt (plasma_context_t *plasma);
void plasma_pcpotrf (plasma_context_t *plasma);
void plasma_pcshift (plasma_context_t *plasma);
void plasma_pcsymm  (plasma_context_t *plasma);
void plasma_pcsyrk  (plasma_context_t *plasma);
void plasma_pcsyr2k (plasma_context_t *plasma);
void plasma_pctrmm  (plasma_context_t *plasma);
void plasma_pctrsm  (plasma_context_t *plasma);
void plasma_pctrsmpl(plasma_context_t *plasma);
void plasma_pcunglq (plasma_context_t *plasma);
void plasma_pcungqr (plasma_context_t *plasma);
void plasma_pcungqrrh(plasma_context_t *plasma);
void plasma_pcunmlq (plasma_context_t *plasma);
void plasma_pcunmqr (plasma_context_t *plasma);
void plasma_pcunpack(plasma_context_t *plasma);

/***************************************************************************//**
 *  Declarations of internal sequential functions 
 **/
int plasma_cshift(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A,
                  int nprob, int me, int ne, int L,
                  PLASMA_sequence *sequence, PLASMA_request *request);

/***************************************************************************//**
 *  Declarations of parallel functions (dynamic scheduling) - alphabetical order
 **/
void plasma_pcaxpy_quark(PLASMA_Complex32_t alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcbarrier_tl2pnl_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcbarrier_pnl2tl_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgelqf_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgelqfrh_quark(PLASMA_desc A, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgemm_quark(PLASMA_enum transA, PLASMA_enum transB, PLASMA_Complex32_t alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_Complex32_t beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgeqrf_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgeqrfrh_quark(PLASMA_desc A, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgerbh_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgerbb_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgerbbrh_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgetmi2_quark(PLASMA_enum idep, PLASMA_enum odep, PLASMA_enum storev, int m, int n, int mb, int nb, PLASMA_Complex32_t *A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgetrf_incpiv_quark(PLASMA_desc A, PLASMA_desc L, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgetrf_reclap_quark(PLASMA_desc A, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgetrf_rectil_quark(PLASMA_desc A, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
#ifdef COMPLEX
void plasma_pchegst_quark(PLASMA_enum itype, PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pchemm_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_Complex32_t alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_Complex32_t beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcherk_quark(PLASMA_enum uplo, PLASMA_enum trans, float alpha, PLASMA_desc A, float beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcher2k_quark(PLASMA_enum uplo, PLASMA_enum trans, PLASMA_Complex32_t alpha, PLASMA_desc A, PLASMA_desc B, float beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcherbt_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgbrdb_quark(PLASMA_enum uplo, PLASMA_desc A, float *D, float *E, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pchbrdt_quark(PLASMA_enum uplo, PLASMA_desc A, float *D, float *E, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
void plasma_pclacpy_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pclag2z_quark(PLASMA_desc A, PLASMA_desc SB, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pclange_quark(PLASMA_enum norm, PLASMA_desc A, float *work, float *result, PLASMA_sequence *sequence, PLASMA_request *request);
#ifdef COMPLEX
void plasma_pclanhe_quark(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc A, float *work, float *result, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
void plasma_pclansy_quark(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc A, float *work, float *result, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pclaset_quark( PLASMA_enum uplo, PLASMA_Complex32_t alpha, PLASMA_Complex32_t beta, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pclaset2_quark(PLASMA_enum uplo, PLASMA_Complex32_t alpha,                          PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pclaswp_quark(PLASMA_desc B, int *IPIV, int inc, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pclauum_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcplghe_quark(float bump, PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcplgsy_quark(PLASMA_Complex32_t bump, PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pcplrnt_quark(PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pcpotrf_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcshift_quark(int, int, int, PLASMA_Complex32_t *, int *, int, int, PLASMA_sequence*, PLASMA_request*);
void plasma_pcsymm_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_Complex32_t alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_Complex32_t beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcsyrk_quark(PLASMA_enum uplo, PLASMA_enum trans, PLASMA_Complex32_t alpha, PLASMA_desc A, PLASMA_Complex32_t beta,  PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcsyr2k_quark(PLASMA_enum uplo, PLASMA_enum trans, PLASMA_Complex32_t alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_Complex32_t beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pctrmm_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, PLASMA_Complex32_t alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pctrsm_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, PLASMA_Complex32_t alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pctrsmpl_quark(PLASMA_desc A, PLASMA_desc B, PLASMA_desc L, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pctrtri_quark(PLASMA_enum uplo, PLASMA_enum diag, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcungbr_quark(PLASMA_enum side, PLASMA_desc A, PLASMA_desc O, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcungbrrh_quark(PLASMA_enum side, PLASMA_desc A, PLASMA_desc O, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcungqr_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcungqrrh_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcunglq_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcunglqrh_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcungtr_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcunmqr_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcunmqrrh_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcunmlq_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcunmlqrh_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
