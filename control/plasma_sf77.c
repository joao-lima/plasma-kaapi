/**
 *
 * @file plasma_sf77.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Bilel Hadri
 * @date 2010-11-15
 * @generated s Thu Sep 15 12:09:27 2011
 *
 **/
#include <stdlib.h>
#include "common.h"
#undef COMPLEX
#define REAL

/*
 * Lapack interface
 */
#define PLASMA_SGEBRD        PLASMA_FNAME(sgebrd,  SGEBRD )
#define PLASMA_SGEEV         PLASMA_FNAME(sgeev,   SGEEV  )
#define PLASMA_SGEHRD        PLASMA_FNAME(sgehrd,  SGEHRD )
#define PLASMA_SGELQF        PLASMA_FNAME(sgelqf,  SGELQF )
#define PLASMA_SGELQS        PLASMA_FNAME(sgelqs,  SGELQS )
#define PLASMA_SGELS         PLASMA_FNAME(sgels,   SGELS  )
#define PLASMA_SGEQRF        PLASMA_FNAME(sgeqrf,  SGEQRF )
#define PLASMA_SGEQRS        PLASMA_FNAME(sgeqrs,  SGEQRS )
#define PLASMA_SGESV         PLASMA_FNAME(sgesv,   SGESV  )
#define PLASMA_SGESVD        PLASMA_FNAME(sgesvd,  SGESVD )
#define PLASMA_SGETRF        PLASMA_FNAME(sgetrf,  SGETRF )
#define PLASMA_SGETRS        PLASMA_FNAME(sgetrs,  SGETRS )
#define PLASMA_SGESV_INCPIV  PLASMA_FNAME(sgesv_incpiv,  SGESV_INCPIV  )
#define PLASMA_SGETRF_INCPIV PLASMA_FNAME(sgetrf_incpiv, SGETRF_INCPIV )
#define PLASMA_SGETRS_INCPIV PLASMA_FNAME(sgetrs_incpiv, SGETRS_INCPIV )
#define PLASMA_SSYEV         PLASMA_FNAME(ssyev,   SSYEV  )
#define PLASMA_SSYGV         PLASMA_FNAME(ssygv,   SSYGV  )
#define PLASMA_SSYGST        PLASMA_FNAME(ssygst,  SSYGST )
#define PLASMA_SSYTRD        PLASMA_FNAME(ssytrd,  SSYTRD )
#define PLASMA_SPOSV         PLASMA_FNAME(sposv,   SPOSV  )
#define PLASMA_SPOTRF        PLASMA_FNAME(spotrf,  SPOTRF )
#define PLASMA_SPOTRI        PLASMA_FNAME(spotri,  SPOTRI )
#define PLASMA_SPOTRS        PLASMA_FNAME(spotrs,  SPOTRS )
#define PLASMA_STRSMPL       PLASMA_FNAME(strsmpl, STRSMPL)
#define PLASMA_SORGBR        PLASMA_FNAME(sorgbr,  SORGBR )
#define PLASMA_SORGHR        PLASMA_FNAME(sorghr,  SORGHR )
#define PLASMA_SORGLQ        PLASMA_FNAME(sorglq,  SORGLQ )
#define PLASMA_SORGQR        PLASMA_FNAME(sorgqr,  SORGQR )
#define PLASMA_SORGTR        PLASMA_FNAME(sorgtr,  SORGTR )
#define PLASMA_SORMLQ        PLASMA_FNAME(sormlq,  SORMLQ )
#define PLASMA_SORMQR        PLASMA_FNAME(sormqr,  SORMQR )
#define PLASMA_STRSM         PLASMA_FNAME(strsm,   STRSM  )
#define PLASMA_SGEMM         PLASMA_FNAME(sgemm,   SGEMM  )
#define PLASMA_SSYMM         PLASMA_FNAME(ssymm,   SSYMM  )
#define PLASMA_SSYRK         PLASMA_FNAME(ssyrk,   SSYRK  )
#ifdef COMPLEX
#define PLASMA_SSYMM         PLASMA_FNAME(ssymm,   SSYMM  )
#define PLASMA_SSYRK         PLASMA_FNAME(ssyrk,   SSYRK  )
#endif

/*
 * Tile interface
 */
#define PLASMA_SGEBRD_TILE        PLASMA_TILE_FNAME(sgebrd,  SGEBRD )
#define PLASMA_SGEEV_TILE         PLASMA_TILE_FNAME(sgeev,   SGEEV  )
#define PLASMA_SGEHRD_TILE        PLASMA_TILE_FNAME(sgehrd,  SGEHRD )
#define PLASMA_SGELQF_TILE        PLASMA_TILE_FNAME(sgelqf,  SGELQF )
#define PLASMA_SGELQS_TILE        PLASMA_TILE_FNAME(sgelqs,  SGELQS )
#define PLASMA_SGELS_TILE         PLASMA_TILE_FNAME(sgels,   SGELS  )
#define PLASMA_SGEQRF_TILE        PLASMA_TILE_FNAME(sgeqrf,  SGEQRF )
#define PLASMA_SGEQRS_TILE        PLASMA_TILE_FNAME(sgeqrs,  SGEQRS )
#define PLASMA_SGESV_TILE         PLASMA_TILE_FNAME(sgesv,   SGESV  )
#define PLASMA_SGESVD_TILE        PLASMA_TILE_FNAME(sgesvd,  SGESVD )
#define PLASMA_SGETRF_TILE        PLASMA_TILE_FNAME(sgetrf,  SGETRF )
#define PLASMA_SGETRS_TILE        PLASMA_TILE_FNAME(sgetrs,  SGETRS ) 
#define PLASMA_SGESV_INCPIV_TILE  PLASMA_TILE_FNAME(sgesv_incpiv,  SGESV_INCPIV  )
#define PLASMA_SGETRF_INCPIV_TILE PLASMA_TILE_FNAME(sgetrf_incpiv, SGETRF_INCPIV )
#define PLASMA_SGETRS_INCPIV_TILE PLASMA_TILE_FNAME(sgetrs_incpiv, SGETRS_INCPIV ) 
#define PLASMA_SSYEV_TILE         PLASMA_TILE_FNAME(ssyev,   SSYEV  )
#define PLASMA_SSYGV_TILE         PLASMA_TILE_FNAME(ssygv,   SSYGV  )
#define PLASMA_SSYGST_TILE        PLASMA_TILE_FNAME(ssygst,  SSYGST )
#define PLASMA_SSYTRD_TILE        PLASMA_TILE_FNAME(ssytrd,  SSYTRD )
#define PLASMA_SPOSV_TILE         PLASMA_TILE_FNAME(sposv,   SPOSV  )
#define PLASMA_SPOTRF_TILE        PLASMA_TILE_FNAME(spotrf,  SPOTRF )
#define PLASMA_SPOTRI_TILE        PLASMA_TILE_FNAME(spotri,  SPOTRI )
#define PLASMA_SPOTRS_TILE        PLASMA_TILE_FNAME(spotrs,  SPOTRS )
#define PLASMA_STRSM_TILE         PLASMA_TILE_FNAME(strsm,   STRSM  )
#define PLASMA_STRSMPL_TILE       PLASMA_TILE_FNAME(strsmpl, STRSMPL)
#define PLASMA_SORGBR_TILE        PLASMA_TILE_FNAME(sorgbr,  SORGBR )
#define PLASMA_SORGHR_TILE        PLASMA_TILE_FNAME(sorghr,  SORGHR )
#define PLASMA_SORGLQ_TILE        PLASMA_TILE_FNAME(sorglq,  SORGLQ )
#define PLASMA_SORGQR_TILE        PLASMA_TILE_FNAME(sorgqr,  SORGQR )
#define PLASMA_SORGTR_TILE        PLASMA_TILE_FNAME(sorgtr,  SORGTR )
#define PLASMA_SORMLQ_TILE        PLASMA_TILE_FNAME(sormlq,  SORMLQ )
#define PLASMA_SORMQR_TILE        PLASMA_TILE_FNAME(sormqr,  SORMQR )
#define PLASMA_SGEMM_TILE         PLASMA_TILE_FNAME(sgemm,   SGEMM  )
#define PLASMA_SSYMM_TILE         PLASMA_TILE_FNAME(ssymm,   SSYMM  )
#define PLASMA_SSYRK_TILE         PLASMA_TILE_FNAME(ssyrk,   SSYRK  )
#ifdef COMPLEX                    
#define PLASMA_SSYMM_TILE         PLASMA_TILE_FNAME(ssymm,   SSYMM  )
#define PLASMA_SSYRK_TILE         PLASMA_TILE_FNAME(ssyrk,   SSYRK  )
#endif

/*
 * Workspaces
 */
#define PLASMA_ALLOC_WORKSPACE_SGEBRD        PLASMA_WS_FNAME(sgehrd, SGEBRD)
#define PLASMA_ALLOC_WORKSPACE_SGEEV         PLASMA_WS_FNAME(sgeev,  SGEEV )
#define PLASMA_ALLOC_WORKSPACE_SGEHRD        PLASMA_WS_FNAME(sgehrd, SGEHRD)
#define PLASMA_ALLOC_WORKSPACE_SGELQF        PLASMA_WS_FNAME(sgelqf, SGELQF) 
#define PLASMA_ALLOC_WORKSPACE_SGELS         PLASMA_WS_FNAME(sgels,  SGELS )  
#define PLASMA_ALLOC_WORKSPACE_SGEQRF        PLASMA_WS_FNAME(sgeqrf, SGEQRF) 
#define PLASMA_ALLOC_WORKSPACE_SGESV_INCPIV  PLASMA_WS_FNAME(sgesv_incpiv,  SGESV_INCPIV ) 
#define PLASMA_ALLOC_WORKSPACE_SGETRF_INCPIV PLASMA_WS_FNAME(sgetrf_incpiv, SGETRF_INCPIV) 
#define PLASMA_ALLOC_WORKSPACE_SGESVD        PLASMA_WS_FNAME(sgesvd, SGESVD)
#define PLASMA_ALLOC_WORKSPACE_SSYEV         PLASMA_WS_FNAME(ssyev,  SSYEV )
#define PLASMA_ALLOC_WORKSPACE_SSYGV         PLASMA_WS_FNAME(ssygv,  SSYGV )
#define PLASMA_ALLOC_WORKSPACE_SSYTRD        PLASMA_WS_FNAME(ssytrd, SSYTRD)

/*
 * Tile Workspaces
 */
#define PLASMA_ALLOC_WORKSPACE_SGELQF_TILE        PLASMA_WST_FNAME(sgelqf,        SGELQF       )
#define PLASMA_ALLOC_WORKSPACE_SGELS_TILE         PLASMA_WST_FNAME(sgels,         SGELS        )
#define PLASMA_ALLOC_WORKSPACE_SGEQRF_TILE        PLASMA_WST_FNAME(sgeqrf,        SGEQRF       )
#define PLASMA_ALLOC_WORKSPACE_SGESV_INCPIV_TILE  PLASMA_WST_FNAME(sgesv_incpiv,  SGESV_INCPIV )
#define PLASMA_ALLOC_WORKSPACE_SGETRF_INCPIV_TILE PLASMA_WST_FNAME(sgetrf_incpiv, SGETRF_INCPIV)

#define PLASMA_SLAPACK_TO_TILE   PLASMA_FNAME(slapack_to_tile, SLAPACK_TO_TILE)
#define PLASMA_STILE_TO_LAPACK   PLASMA_FNAME(stile_to_lapack, STILE_TO_LAPACK)

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  FORTRAN API - math functions (simple interface)
 **/
void PLASMA_SGEBRD(PLASMA_enum *jobu, PLASMA_enum *jobvt, int *M, int *N, float *A, int *LDA, float *D, float *E, float *U, int *LDU, float *VT, int *LDVT, intptr_t *descT, int *INFO)
{   *INFO = PLASMA_sgebrd(*jobu, *jobvt, *M, *N, A, *LDA, D, E, U, *LDU, VT, *LDVT, (PLASMA_desc *)(*descT)); }

void PLASMA_SGELQF(int *M, int *N, float *A, int *LDA, float **T, int *INFO)
{   *INFO = PLASMA_sgelqf(*M, *N, A, *LDA, *T); }

void PLASMA_SGELQS(int *M, int *N, int *NRHS, float *A, int *LDA, float **T, float *B, int *LDB, int *INFO)
{   *INFO = PLASMA_sgelqs(*M, *N, *NRHS, A, *LDA, *T, B, *LDB); }

void PLASMA_SGELS(PLASMA_enum *trans, int *M, int *N, int *NRHS, float *A, int *LDA, float **T, float *B, int *LDB, int *INFO)
{   *INFO = PLASMA_sgels(*trans, *M, *N, *NRHS, A, *LDA, *T, B, *LDB); }

void PLASMA_SGEQRF(int *M, int *N, float *A, int *LDA, float **T, int *INFO)
{   *INFO = PLASMA_sgeqrf(*M, *N, A, *LDA, *T); }

void PLASMA_SGEQRS(int *M, int *N, int *NRHS, float *A, int *LDA, float **T, float *B, int *LDB, int *INFO)
{   *INFO = PLASMA_sgeqrs(*M, *N, *NRHS, A, *LDA, *T, B, *LDB); }

void PLASMA_SGESV(int *N, int *NRHS, float *A, int *LDA, int *IPIV, float *B, int *LDB, int *INFO)
{   *INFO = PLASMA_sgesv(*N, *NRHS, A, *LDA, IPIV, B, *LDB); }

void PLASMA_SGESVD(PLASMA_enum *jobu, PLASMA_enum *jobvt, int *M, int *N, float *A, int *LDA, float *S, float *U, int *LDU, float *VT, int *LDVT, intptr_t *T, int *INFO)
{   *INFO = PLASMA_sgesvd(*jobu, *jobvt, *M, *N, A, *LDA, S, U, *LDU, VT, *LDVT, (PLASMA_desc *)(*T)); }

void PLASMA_SGETRF(int *M, int *N, float *A, int *LDA, int *IPIV, int *INFO)
{   *INFO = PLASMA_sgetrf(*M, *N, A, *LDA, IPIV); }

void PLASMA_SGETRS(PLASMA_enum *trans, int *N, int *NRHS, float *A, int *LDA, int *IPIV, float *B, int *LDB, int *INFO)
{   *INFO = PLASMA_sgetrs(*trans, *N, *NRHS, A, *LDA, IPIV, B, *LDB); }

void PLASMA_SGESV_INCPIV(int *N, int *NRHS, float *A, int *LDA, float **LH, int **IPIVH, float *B, int *LDB, int *INFO)
{   *INFO = PLASMA_sgesv_incpiv(*N, *NRHS, A, *LDA, *LH, *IPIVH, B, *LDB); }

void PLASMA_SGETRF_INCPIV(int *M, int *N, float *A, int *LDA, float **LH, int **IPIVH, int *INFO)
{   *INFO = PLASMA_sgetrf_incpiv(*M, *N, A, *LDA, *LH, *IPIVH); }

void PLASMA_SGETRS_INCPIV(PLASMA_enum *uplo, int *N, int *NRHS, float *A, int *LDA, float **LH, int **IPIVH, float *B, int *LDB, int *INFO)
{   *INFO = PLASMA_sgetrs_incpiv(*uplo, *N, *NRHS, A, *LDA, *LH, *IPIVH, B, *LDB); }

void PLASMA_SSYEV(PLASMA_enum *jobz, PLASMA_enum *uplo, int *N, float *A, int *LDA, float *W, intptr_t *T, float *Q, int *LDQ, int *INFO)
{   *INFO = PLASMA_ssyev(*jobz, *uplo, *N, A, *LDA, W, (PLASMA_desc *)(*T), Q, *LDQ); }

void PLASMA_SSYGV(PLASMA_enum *itype, PLASMA_enum *jobz, PLASMA_enum *uplo, int *N, float *A, int *LDA, float *B, int *LDB, float *W, intptr_t *T, float *Q, int *LDQ, int *INFO)
{   *INFO = PLASMA_ssygv(*itype, *jobz, *uplo, *N, A, *LDA, B, *LDB, W, (PLASMA_desc*)(*T), Q, *LDQ); }

void PLASMA_SSYGST(PLASMA_enum *itype, PLASMA_enum *uplo, int *N, float *A, int *LDA, float *B, int *LDB, int *INFO)
{   *INFO = PLASMA_ssygst(*itype, *uplo, *N, A, *LDA, B, *LDB); }

void PLASMA_SSYTRD(PLASMA_enum *jobz, PLASMA_enum *uplo, int *N, float *A, int *LDA, float *D, float *E, intptr_t *T, float *Q, int *LDQ, int *INFO)
{   *INFO = PLASMA_ssytrd(*jobz, *uplo, *N, A, *LDA, D, E, (PLASMA_desc *)(*T), Q, *LDQ); }

void PLASMA_SPOSV(PLASMA_enum *uplo, int *N, int *NRHS, float *A, int *LDA, float *B, int *LDB, int *INFO)
{   *INFO = PLASMA_sposv(*uplo, *N, *NRHS, A, *LDA, B, *LDB); }

void PLASMA_SPOTRF(PLASMA_enum *uplo, int *N, float *A, int *LDA, int *INFO)
{   *INFO = PLASMA_spotrf(*uplo, *N, A, *LDA); }

void PLASMA_SPOTRI(PLASMA_enum *uplo, int *N, float *A, int *LDA, int *INFO)
{   *INFO = PLASMA_spotri(*uplo, *N, A, *LDA); }

void PLASMA_SPOTRS(PLASMA_enum *uplo, int *N, int *NRHS, float *A, int *LDA, float *B, int* LDB, int * INFO)
{   *INFO = PLASMA_spotrs(*uplo, *N, *NRHS, A, *LDA, B, *LDB); }

void PLASMA_STRSMPL(int *N, int *NRHS, float *A, int *LDA, float **LH, int **IPIVH, float *B, int *LDB, int *INFO)
{   *INFO = PLASMA_strsmpl(*N, *NRHS, A, *LDA, *LH, *IPIVH, B, *LDB); }

void PLASMA_SORGLQ(int *M, int *N, int *K, float *A, int *LDA, float **T, float *B, int *LDB, int *INFO)
{   *INFO = PLASMA_sorglq(*M, *N, *K, A, *LDA, *T, B, *LDB); }

void PLASMA_SORGQR(int *M, int *N, int *K, float *A, int *LDA, float **T, float *B, int *LDB, int *INFO)
{   *INFO = PLASMA_sorgqr(*M, *N, *K, A, *LDA, *T, B, *LDB); }

void PLASMA_SORMLQ(PLASMA_enum *side, PLASMA_enum *trans, int *M, int *N, int *K, float *A, int *LDA, float **T, float *B, int *LDB, int *INFO)
{   *INFO = PLASMA_sormlq(*side, *trans, *M, *N, *K, A, *LDA, *T, B, *LDB); }

void PLASMA_SORMQR(PLASMA_enum *side, PLASMA_enum *trans, int *M, int *N, int *K, float *A, int *LDA, float **T, float *B, int *LDB, int *INFO)
{   *INFO = PLASMA_sormqr(*side, *trans, *M, *N, *K, A, *LDA, *T, B, *LDB); }

void PLASMA_STRSM(PLASMA_enum *side, PLASMA_enum *uplo, PLASMA_enum *transA, PLASMA_enum *diag, int *N, int *NRHS, float *alpha, float *A, int *LDA, float *B, int *LDB, int *INFO)
{   *INFO = PLASMA_strsm(*side, *uplo, *transA, *diag, *N, *NRHS, *alpha, A, *LDA, B, *LDB); }

void PLASMA_SGEMM(PLASMA_enum *transA, PLASMA_enum *transB, int *M, int *N, int *K, float *alpha, float *A, int *LDA, float *B, int *LDB, float *beta, float *C, int *LDC, int *INFO)
{   *INFO = PLASMA_sgemm(*transA, *transB, *M, *N, *K, *alpha, A, *LDA, B, *LDB, *beta, C, *LDC); }

void PLASMA_SSYMM(PLASMA_enum *side, PLASMA_enum *uplo, int *M, int *N, float *alpha, float *A, int *LDA, float *B, int *LDB, float *beta, float *C, int *LDC, int *INFO)
{   *INFO = PLASMA_ssymm(*side, *uplo, *M, *N, *alpha, A, *LDA, B, *LDB, *beta, C, *LDC); }

void PLASMA_SSYRK(PLASMA_enum *uplo, PLASMA_enum *trans, int *N, int *K, float *alpha, float *A, int *LDA, float *beta, float *C, int *LDC, int *INFO)
{   *INFO = PLASMA_ssyrk(*uplo, *trans, *N, *K, *alpha, A, *LDA, *beta, C, *LDC); }

#ifdef COMPLEX
void PLASMA_SSYMM(PLASMA_enum *side, PLASMA_enum *uplo, int *M, int *N, float *alpha, float *A, int *LDA, float *B, int *LDB, float *beta, float *C, int *LDC, int *INFO)
{   *INFO = PLASMA_ssymm(*side, *uplo, *M, *N, *alpha, A, *LDA, B, *LDB, *beta, C, *LDC); }

void PLASMA_SSYRK(PLASMA_enum *uplo, PLASMA_enum *trans, int *N, int *K, float *alpha, float *A, int *LDA, float *beta, float *C, int *LDC, int *INFO)
{   *INFO = PLASMA_ssyrk(*uplo, *trans, *N, *K, *alpha, A, *LDA, *beta, C, *LDC); }
#endif

/***************************************************************************//**
 *  FORTRAN API - math functions (native interface)
 **/
void PLASMA_SGEBRD_TILE(PLASMA_enum *jobu, PLASMA_enum *jobvt, intptr_t *A, float *D, float *E, intptr_t *U, intptr_t *VT, intptr_t *T, int *INFO)
{   *INFO = PLASMA_sgebrd_Tile(*jobu, *jobvt, (PLASMA_desc *)(*A), D, E, (PLASMA_desc *)(*U), (PLASMA_desc *)(*VT), (PLASMA_desc *)(*T)); }

void PLASMA_SGELQF_TILE(intptr_t *A, intptr_t *T, int *INFO)
{   *INFO = PLASMA_sgelqf_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*T)); }

void PLASMA_SGELQS_TILE(intptr_t *A, intptr_t *B, intptr_t *T, int *INFO)
{   *INFO = PLASMA_sgelqs_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*B), (PLASMA_desc *)(*T)); }

void PLASMA_SGELS_TILE(PLASMA_enum *trans, intptr_t *A, intptr_t *B, intptr_t *T, int *INFO)
{   *INFO = PLASMA_sgels_Tile(*trans, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B), (PLASMA_desc *)(*T)); }

void PLASMA_SGEQRF_TILE(intptr_t *A, intptr_t *T, int *INFO)
{   *INFO = PLASMA_sgeqrf_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*T)); }

void PLASMA_SGEQRS_TILE(intptr_t *A, intptr_t *B, intptr_t *T, int *INFO)
{   *INFO = PLASMA_sgeqrs_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*B), (PLASMA_desc *)(*T)); }

void PLASMA_SGESV_TILE(intptr_t *A, int *IPIV, intptr_t *B, int *INFO)
{   *INFO = PLASMA_sgesv_Tile((PLASMA_desc *)(*A), IPIV, (PLASMA_desc *)(*B)); }

void PLASMA_SGESVD_TILE(PLASMA_enum *jobu, PLASMA_enum *jobvt, intptr_t *A, float *S, intptr_t *U, intptr_t *VT, intptr_t *T, int *INFO)
{   *INFO = PLASMA_sgesvd_Tile(*jobu, *jobvt, (PLASMA_desc *)(*A), S, (PLASMA_desc *)(*U), (PLASMA_desc *)(*VT), (PLASMA_desc *)(*T)); }

void PLASMA_SGETRF_TILE(intptr_t *A, int *IPIV, int *INFO)
{   *INFO = PLASMA_sgetrf_Tile((PLASMA_desc *)(*A), IPIV); }

void PLASMA_SGETRS_TILE(PLASMA_enum *trans, intptr_t *A, int *IPIV, intptr_t *B, int *INFO)
{   *INFO = PLASMA_sgetrs_Tile(*trans, (PLASMA_desc *)(*A), IPIV, (PLASMA_desc *)(*B)); }

void PLASMA_SGESV_INCPIV_TILE(intptr_t *A, intptr_t *L, int **IPIVH, intptr_t *B, int *INFO)
{   *INFO = PLASMA_sgesv_incpiv_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*L), *IPIVH, (PLASMA_desc *)(*B)); }

void PLASMA_SGETRF_INCPIV_TILE(intptr_t *A, intptr_t *L, int **IPIVH, int *INFO)
{   *INFO = PLASMA_sgetrf_incpiv_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*L), *IPIVH); }

void PLASMA_SGETRS_INCPIV_TILE(intptr_t *A, intptr_t *L, int **IPIVH, intptr_t *B, int *INFO)
{   *INFO = PLASMA_sgetrs_incpiv_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*L), *IPIVH, (PLASMA_desc *)(*B)); }

void PLASMA_SSYEV_TILE(PLASMA_enum *jobz, PLASMA_enum *uplo, intptr_t *A, float *W, intptr_t *T, intptr_t *Q, int *INFO)
{   *INFO = PLASMA_ssyev_Tile(*jobz, *uplo, (PLASMA_desc *)(*A), W, (PLASMA_desc *)(*T), (PLASMA_desc *)(*Q)); }

void PLASMA_SSYGV_TILE(PLASMA_enum *itype, PLASMA_enum *jobz, PLASMA_enum *uplo, intptr_t *A, intptr_t *B, float *W, intptr_t *T, intptr_t *Q, int *INFO)
{   *INFO = PLASMA_ssygv_Tile(*itype, *jobz, *uplo, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B), W, (PLASMA_desc *)(*T), (PLASMA_desc *)(*Q)); }

void PLASMA_SSYGST_TILE(PLASMA_enum *itype, PLASMA_enum *uplo, intptr_t *A, intptr_t *B, int *INFO)
{   *INFO = PLASMA_ssygst_Tile(*itype, *uplo, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B)); }

void PLASMA_SSYTRD_TILE(PLASMA_enum *jobz, PLASMA_enum *uplo, intptr_t *A, float *D, float *E, intptr_t *T, intptr_t *Q, int *INFO)
{   *INFO = PLASMA_ssytrd_Tile(*jobz, *uplo, (PLASMA_desc *)(*A), D, E, (PLASMA_desc *)(*T), (PLASMA_desc *)(*Q)); }

void PLASMA_SPOSV_TILE(PLASMA_enum *uplo, intptr_t *A, intptr_t *B, int *INFO)
{   *INFO = PLASMA_sposv_Tile(*uplo, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B)); }

void PLASMA_SPOTRF_TILE(PLASMA_enum *uplo, intptr_t *A, int *INFO)
{   *INFO = PLASMA_spotrf_Tile(*uplo, (PLASMA_desc *)(*A)); }

void PLASMA_SPOTRS_TILE(PLASMA_enum *uplo, intptr_t *A, intptr_t *B, int *INFO)
{   *INFO = PLASMA_spotrs_Tile(*uplo, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B)); }

void PLASMA_STRSMPL_TILE(intptr_t *A, intptr_t *L, int **IPIVH, intptr_t *B, int *INFO)
{   *INFO = PLASMA_strsmpl_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*L), *IPIVH, (PLASMA_desc *)(*B)); }

void PLASMA_SORGLQ_TILE(intptr_t *A, intptr_t *T, intptr_t *B, int *INFO)
{   *INFO = PLASMA_sorglq_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*T), (PLASMA_desc *)(*B)); }

void PLASMA_SORGQR_TILE(intptr_t *A, intptr_t *T, intptr_t *B, int *INFO)
{   *INFO = PLASMA_sorgqr_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*T), (PLASMA_desc *)(*B)); }

void PLASMA_SORMLQ_TILE(PLASMA_enum *side, PLASMA_enum *trans, intptr_t *A, intptr_t *T, intptr_t *B, int *INFO)
{   *INFO = PLASMA_sormlq_Tile(*side, *trans, (PLASMA_desc *)(*A), (PLASMA_desc *)(*T), (PLASMA_desc *)(*B)); }

void PLASMA_SORMQR_TILE(PLASMA_enum *side, PLASMA_enum *trans, intptr_t *A, intptr_t *T, intptr_t *B, int *INFO)
{   *INFO = PLASMA_sormqr_Tile(*side, *trans, (PLASMA_desc *)(*A), (PLASMA_desc *)(*T), (PLASMA_desc *)(*B)); }

void PLASMA_STRSM_TILE(PLASMA_enum *side, PLASMA_enum *uplo, PLASMA_enum *transA, PLASMA_enum *diag, float *alpha, intptr_t *A, intptr_t *B, int *INFO)
{   *INFO = PLASMA_strsm_Tile(*side, *uplo, *transA, *diag, *alpha, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B)); }

void PLASMA_SGEMM_TILE(PLASMA_enum *transA, PLASMA_enum *transB, int *alpha, intptr_t *A, intptr_t *B, int *beta, intptr_t *C, int *INFO)
{   *INFO = PLASMA_sgemm_Tile(*transA, *transB, *alpha, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B), *beta, (PLASMA_desc *)(*C)); }

  void PLASMA_SSYMM_TILE(PLASMA_enum *side, PLASMA_enum *uplo, float *alpha, intptr_t *A, intptr_t *B, float *beta, intptr_t *C, int *INFO)
  {   *INFO = PLASMA_ssymm_Tile(*side, *uplo, *alpha, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B), *beta, (PLASMA_desc *)(*C)); }

void PLASMA_SSYRK_TILE(PLASMA_enum *uplo, PLASMA_enum *trans, float *alpha, intptr_t *A, float *beta, intptr_t *C, int *INFO)
{   *INFO = PLASMA_ssyrk_Tile(*uplo, *trans, *alpha, (PLASMA_desc *)(*A), *beta, (PLASMA_desc *)(*C)); }

#ifdef COMPLEX
void PLASMA_SSYMM_TILE(PLASMA_enum *side, PLASMA_enum *uplo, float *alpha, intptr_t *A, intptr_t *B, float *beta, intptr_t *C, int *INFO)
{   *INFO = PLASMA_ssymm_Tile(*side, *uplo, *alpha, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B), *beta, (PLASMA_desc *)(*C)); }

void PLASMA_SSYRK_TILE(PLASMA_enum *uplo, PLASMA_enum *trans, float *alpha, intptr_t *A, float *beta, intptr_t *C, int *INFO)
{   *INFO = PLASMA_ssyrk_Tile(*uplo, *trans, *alpha, (PLASMA_desc *)(*A), *beta, (PLASMA_desc *)(*C)); }
#endif

/***************************************************************************//**
 *  FORTRAN API - workspace allocation
 **/
void PLASMA_ALLOC_WORKSPACE_SGEBRD(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_sgebrd(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_SGELQF(int *M, int *N, float **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_sgelqf(*M, *N, T); }

void PLASMA_ALLOC_WORKSPACE_SGELS(int *M, int *N, float **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_sgels(*M, *N, T); }

void PLASMA_ALLOC_WORKSPACE_SGEQRF(int *M, int *N, float **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_sgeqrf(*M, *N, T); }

void PLASMA_ALLOC_WORKSPACE_SGESV_INCPIV(int *N, float **L, int **IPIV, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_sgesv_incpiv(*N, L, IPIV); }

void PLASMA_ALLOC_WORKSPACE_SGESVD(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_sgesvd(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_SGETRF_INCPIV(int *M, int *N, float **L, int **IPIV, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_sgetrf_incpiv(*M, *N, L, IPIV); }

void PLASMA_ALLOC_WORKSPACE_SSYEV(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_ssyev(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_SSYGV(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_ssygv(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_SSYTRD(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_ssytrd(*M, *N, (PLASMA_desc **)T); }


/***************************************************************************//**
 *  FORTRAN API - tiled workspace allocation
 **/
void PLASMA_ALLOC_WORKSPACE_SGELQF_TILE(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_sgelqf_Tile(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_SGELS_TILE(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_sgels_Tile(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_SGEQRF_TILE(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_sgeqrf_Tile(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_SGESV_INCPIV_TILE(int *N, intptr_t **L, int **IPIV, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_sgesv_incpiv_Tile(*N, (PLASMA_desc **)L, IPIV); }

void PLASMA_ALLOC_WORKSPACE_SGETRF_INCPIV_TILE(int *N, intptr_t **L, int **IPIV, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_sgetrf_incpiv_Tile(*N, (PLASMA_desc **)L, IPIV); }

/***************************************************************************//**
 *  FORTRAN API - conversion from LAPACK F77 matrix layout to tile layout
 **/
void PLASMA_SLAPACK_TO_TILE(float **Af77, int *LDA, intptr_t *A, int *INFO)
{   *INFO = PLASMA_sLapack_to_Tile( *Af77, *LDA, (PLASMA_desc *)(*A) ); }

void PLASMA_STILE_TO_LAPACK(intptr_t *A, float **Af77, int *LDA, int *INFO)
{   *INFO = PLASMA_sTile_to_Lapack( (PLASMA_desc *)(*A), *Af77, *LDA ); }

#ifdef __cplusplus
}
#endif
