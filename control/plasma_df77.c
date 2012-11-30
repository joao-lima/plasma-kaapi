/**
 *
 * @file plasma_df77.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Bilel Hadri
 * @date 2010-11-15
 * @generated d Thu Sep 15 12:09:27 2011
 *
 **/
#include <stdlib.h>
#include "common.h"
#undef COMPLEX
#define REAL

/*
 * Lapack interface
 */
#define PLASMA_DGEBRD        PLASMA_FNAME(dgebrd,  DGEBRD )
#define PLASMA_DGEEV         PLASMA_FNAME(dgeev,   DGEEV  )
#define PLASMA_DGEHRD        PLASMA_FNAME(dgehrd,  DGEHRD )
#define PLASMA_DGELQF        PLASMA_FNAME(dgelqf,  DGELQF )
#define PLASMA_DGELQS        PLASMA_FNAME(dgelqs,  DGELQS )
#define PLASMA_DGELS         PLASMA_FNAME(dgels,   DGELS  )
#define PLASMA_DGEQRF        PLASMA_FNAME(dgeqrf,  DGEQRF )
#define PLASMA_DGEQRS        PLASMA_FNAME(dgeqrs,  DGEQRS )
#define PLASMA_DGESV         PLASMA_FNAME(dgesv,   DGESV  )
#define PLASMA_DGESVD        PLASMA_FNAME(dgesvd,  DGESVD )
#define PLASMA_DGETRF        PLASMA_FNAME(dgetrf,  DGETRF )
#define PLASMA_DGETRS        PLASMA_FNAME(dgetrs,  DGETRS )
#define PLASMA_DGESV_INCPIV  PLASMA_FNAME(dgesv_incpiv,  DGESV_INCPIV  )
#define PLASMA_DGETRF_INCPIV PLASMA_FNAME(dgetrf_incpiv, DGETRF_INCPIV )
#define PLASMA_DGETRS_INCPIV PLASMA_FNAME(dgetrs_incpiv, DGETRS_INCPIV )
#define PLASMA_DSYEV         PLASMA_FNAME(dsyev,   DSYEV  )
#define PLASMA_DSYGV         PLASMA_FNAME(dsygv,   DSYGV  )
#define PLASMA_DSYGST        PLASMA_FNAME(dsygst,  DSYGST )
#define PLASMA_DSYTRD        PLASMA_FNAME(dsytrd,  DSYTRD )
#define PLASMA_DPOSV         PLASMA_FNAME(dposv,   DPOSV  )
#define PLASMA_DPOTRF        PLASMA_FNAME(dpotrf,  DPOTRF )
#define PLASMA_DPOTRI        PLASMA_FNAME(dpotri,  DPOTRI )
#define PLASMA_DPOTRS        PLASMA_FNAME(dpotrs,  DPOTRS )
#define PLASMA_DTRSMPL       PLASMA_FNAME(dtrsmpl, DTRSMPL)
#define PLASMA_DORGBR        PLASMA_FNAME(dorgbr,  DORGBR )
#define PLASMA_DORGHR        PLASMA_FNAME(dorghr,  DORGHR )
#define PLASMA_DORGLQ        PLASMA_FNAME(dorglq,  DORGLQ )
#define PLASMA_DORGQR        PLASMA_FNAME(dorgqr,  DORGQR )
#define PLASMA_DORGTR        PLASMA_FNAME(dorgtr,  DORGTR )
#define PLASMA_DORMLQ        PLASMA_FNAME(dormlq,  DORMLQ )
#define PLASMA_DORMQR        PLASMA_FNAME(dormqr,  DORMQR )
#define PLASMA_DTRSM         PLASMA_FNAME(dtrsm,   DTRSM  )
#define PLASMA_DGEMM         PLASMA_FNAME(dgemm,   DGEMM  )
#define PLASMA_DSYMM         PLASMA_FNAME(dsymm,   DSYMM  )
#define PLASMA_DSYRK         PLASMA_FNAME(dsyrk,   DSYRK  )
#ifdef COMPLEX
#define PLASMA_DSYMM         PLASMA_FNAME(dsymm,   DSYMM  )
#define PLASMA_DSYRK         PLASMA_FNAME(dsyrk,   DSYRK  )
#endif

/*
 * Tile interface
 */
#define PLASMA_DGEBRD_TILE        PLASMA_TILE_FNAME(dgebrd,  DGEBRD )
#define PLASMA_DGEEV_TILE         PLASMA_TILE_FNAME(dgeev,   DGEEV  )
#define PLASMA_DGEHRD_TILE        PLASMA_TILE_FNAME(dgehrd,  DGEHRD )
#define PLASMA_DGELQF_TILE        PLASMA_TILE_FNAME(dgelqf,  DGELQF )
#define PLASMA_DGELQS_TILE        PLASMA_TILE_FNAME(dgelqs,  DGELQS )
#define PLASMA_DGELS_TILE         PLASMA_TILE_FNAME(dgels,   DGELS  )
#define PLASMA_DGEQRF_TILE        PLASMA_TILE_FNAME(dgeqrf,  DGEQRF )
#define PLASMA_DGEQRS_TILE        PLASMA_TILE_FNAME(dgeqrs,  DGEQRS )
#define PLASMA_DGESV_TILE         PLASMA_TILE_FNAME(dgesv,   DGESV  )
#define PLASMA_DGESVD_TILE        PLASMA_TILE_FNAME(dgesvd,  DGESVD )
#define PLASMA_DGETRF_TILE        PLASMA_TILE_FNAME(dgetrf,  DGETRF )
#define PLASMA_DGETRS_TILE        PLASMA_TILE_FNAME(dgetrs,  DGETRS ) 
#define PLASMA_DGESV_INCPIV_TILE  PLASMA_TILE_FNAME(dgesv_incpiv,  DGESV_INCPIV  )
#define PLASMA_DGETRF_INCPIV_TILE PLASMA_TILE_FNAME(dgetrf_incpiv, DGETRF_INCPIV )
#define PLASMA_DGETRS_INCPIV_TILE PLASMA_TILE_FNAME(dgetrs_incpiv, DGETRS_INCPIV ) 
#define PLASMA_DSYEV_TILE         PLASMA_TILE_FNAME(dsyev,   DSYEV  )
#define PLASMA_DSYGV_TILE         PLASMA_TILE_FNAME(dsygv,   DSYGV  )
#define PLASMA_DSYGST_TILE        PLASMA_TILE_FNAME(dsygst,  DSYGST )
#define PLASMA_DSYTRD_TILE        PLASMA_TILE_FNAME(dsytrd,  DSYTRD )
#define PLASMA_DPOSV_TILE         PLASMA_TILE_FNAME(dposv,   DPOSV  )
#define PLASMA_DPOTRF_TILE        PLASMA_TILE_FNAME(dpotrf,  DPOTRF )
#define PLASMA_DPOTRI_TILE        PLASMA_TILE_FNAME(dpotri,  DPOTRI )
#define PLASMA_DPOTRS_TILE        PLASMA_TILE_FNAME(dpotrs,  DPOTRS )
#define PLASMA_DTRSM_TILE         PLASMA_TILE_FNAME(dtrsm,   DTRSM  )
#define PLASMA_DTRSMPL_TILE       PLASMA_TILE_FNAME(dtrsmpl, DTRSMPL)
#define PLASMA_DORGBR_TILE        PLASMA_TILE_FNAME(dorgbr,  DORGBR )
#define PLASMA_DORGHR_TILE        PLASMA_TILE_FNAME(dorghr,  DORGHR )
#define PLASMA_DORGLQ_TILE        PLASMA_TILE_FNAME(dorglq,  DORGLQ )
#define PLASMA_DORGQR_TILE        PLASMA_TILE_FNAME(dorgqr,  DORGQR )
#define PLASMA_DORGTR_TILE        PLASMA_TILE_FNAME(dorgtr,  DORGTR )
#define PLASMA_DORMLQ_TILE        PLASMA_TILE_FNAME(dormlq,  DORMLQ )
#define PLASMA_DORMQR_TILE        PLASMA_TILE_FNAME(dormqr,  DORMQR )
#define PLASMA_DGEMM_TILE         PLASMA_TILE_FNAME(dgemm,   DGEMM  )
#define PLASMA_DSYMM_TILE         PLASMA_TILE_FNAME(dsymm,   DSYMM  )
#define PLASMA_DSYRK_TILE         PLASMA_TILE_FNAME(dsyrk,   DSYRK  )
#ifdef COMPLEX                    
#define PLASMA_DSYMM_TILE         PLASMA_TILE_FNAME(dsymm,   DSYMM  )
#define PLASMA_DSYRK_TILE         PLASMA_TILE_FNAME(dsyrk,   DSYRK  )
#endif

/*
 * Workspaces
 */
#define PLASMA_ALLOC_WORKSPACE_DGEBRD        PLASMA_WS_FNAME(dgehrd, DGEBRD)
#define PLASMA_ALLOC_WORKSPACE_DGEEV         PLASMA_WS_FNAME(dgeev,  DGEEV )
#define PLASMA_ALLOC_WORKSPACE_DGEHRD        PLASMA_WS_FNAME(dgehrd, DGEHRD)
#define PLASMA_ALLOC_WORKSPACE_DGELQF        PLASMA_WS_FNAME(dgelqf, DGELQF) 
#define PLASMA_ALLOC_WORKSPACE_DGELS         PLASMA_WS_FNAME(dgels,  DGELS )  
#define PLASMA_ALLOC_WORKSPACE_DGEQRF        PLASMA_WS_FNAME(dgeqrf, DGEQRF) 
#define PLASMA_ALLOC_WORKSPACE_DGESV_INCPIV  PLASMA_WS_FNAME(dgesv_incpiv,  DGESV_INCPIV ) 
#define PLASMA_ALLOC_WORKSPACE_DGETRF_INCPIV PLASMA_WS_FNAME(dgetrf_incpiv, DGETRF_INCPIV) 
#define PLASMA_ALLOC_WORKSPACE_DGESVD        PLASMA_WS_FNAME(dgesvd, DGESVD)
#define PLASMA_ALLOC_WORKSPACE_DSYEV         PLASMA_WS_FNAME(dsyev,  DSYEV )
#define PLASMA_ALLOC_WORKSPACE_DSYGV         PLASMA_WS_FNAME(dsygv,  DSYGV )
#define PLASMA_ALLOC_WORKSPACE_DSYTRD        PLASMA_WS_FNAME(dsytrd, DSYTRD)

/*
 * Tile Workspaces
 */
#define PLASMA_ALLOC_WORKSPACE_DGELQF_TILE        PLASMA_WST_FNAME(dgelqf,        DGELQF       )
#define PLASMA_ALLOC_WORKSPACE_DGELS_TILE         PLASMA_WST_FNAME(dgels,         DGELS        )
#define PLASMA_ALLOC_WORKSPACE_DGEQRF_TILE        PLASMA_WST_FNAME(dgeqrf,        DGEQRF       )
#define PLASMA_ALLOC_WORKSPACE_DGESV_INCPIV_TILE  PLASMA_WST_FNAME(dgesv_incpiv,  DGESV_INCPIV )
#define PLASMA_ALLOC_WORKSPACE_DGETRF_INCPIV_TILE PLASMA_WST_FNAME(dgetrf_incpiv, DGETRF_INCPIV)

#define PLASMA_DLAPACK_TO_TILE   PLASMA_FNAME(dlapack_to_tile, DLAPACK_TO_TILE)
#define PLASMA_DTILE_TO_LAPACK   PLASMA_FNAME(dtile_to_lapack, DTILE_TO_LAPACK)

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  FORTRAN API - math functions (simple interface)
 **/
void PLASMA_DGEBRD(PLASMA_enum *jobu, PLASMA_enum *jobvt, int *M, int *N, double *A, int *LDA, double *D, double *E, double *U, int *LDU, double *VT, int *LDVT, intptr_t *descT, int *INFO)
{   *INFO = PLASMA_dgebrd(*jobu, *jobvt, *M, *N, A, *LDA, D, E, U, *LDU, VT, *LDVT, (PLASMA_desc *)(*descT)); }

void PLASMA_DGELQF(int *M, int *N, double *A, int *LDA, double **T, int *INFO)
{   *INFO = PLASMA_dgelqf(*M, *N, A, *LDA, *T); }

void PLASMA_DGELQS(int *M, int *N, int *NRHS, double *A, int *LDA, double **T, double *B, int *LDB, int *INFO)
{   *INFO = PLASMA_dgelqs(*M, *N, *NRHS, A, *LDA, *T, B, *LDB); }

void PLASMA_DGELS(PLASMA_enum *trans, int *M, int *N, int *NRHS, double *A, int *LDA, double **T, double *B, int *LDB, int *INFO)
{   *INFO = PLASMA_dgels(*trans, *M, *N, *NRHS, A, *LDA, *T, B, *LDB); }

void PLASMA_DGEQRF(int *M, int *N, double *A, int *LDA, double **T, int *INFO)
{   *INFO = PLASMA_dgeqrf(*M, *N, A, *LDA, *T); }

void PLASMA_DGEQRS(int *M, int *N, int *NRHS, double *A, int *LDA, double **T, double *B, int *LDB, int *INFO)
{   *INFO = PLASMA_dgeqrs(*M, *N, *NRHS, A, *LDA, *T, B, *LDB); }

void PLASMA_DGESV(int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO)
{   *INFO = PLASMA_dgesv(*N, *NRHS, A, *LDA, IPIV, B, *LDB); }

void PLASMA_DGESVD(PLASMA_enum *jobu, PLASMA_enum *jobvt, int *M, int *N, double *A, int *LDA, double *S, double *U, int *LDU, double *VT, int *LDVT, intptr_t *T, int *INFO)
{   *INFO = PLASMA_dgesvd(*jobu, *jobvt, *M, *N, A, *LDA, S, U, *LDU, VT, *LDVT, (PLASMA_desc *)(*T)); }

void PLASMA_DGETRF(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO)
{   *INFO = PLASMA_dgetrf(*M, *N, A, *LDA, IPIV); }

void PLASMA_DGETRS(PLASMA_enum *trans, int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO)
{   *INFO = PLASMA_dgetrs(*trans, *N, *NRHS, A, *LDA, IPIV, B, *LDB); }

void PLASMA_DGESV_INCPIV(int *N, int *NRHS, double *A, int *LDA, double **LH, int **IPIVH, double *B, int *LDB, int *INFO)
{   *INFO = PLASMA_dgesv_incpiv(*N, *NRHS, A, *LDA, *LH, *IPIVH, B, *LDB); }

void PLASMA_DGETRF_INCPIV(int *M, int *N, double *A, int *LDA, double **LH, int **IPIVH, int *INFO)
{   *INFO = PLASMA_dgetrf_incpiv(*M, *N, A, *LDA, *LH, *IPIVH); }

void PLASMA_DGETRS_INCPIV(PLASMA_enum *uplo, int *N, int *NRHS, double *A, int *LDA, double **LH, int **IPIVH, double *B, int *LDB, int *INFO)
{   *INFO = PLASMA_dgetrs_incpiv(*uplo, *N, *NRHS, A, *LDA, *LH, *IPIVH, B, *LDB); }

void PLASMA_DSYEV(PLASMA_enum *jobz, PLASMA_enum *uplo, int *N, double *A, int *LDA, double *W, intptr_t *T, double *Q, int *LDQ, int *INFO)
{   *INFO = PLASMA_dsyev(*jobz, *uplo, *N, A, *LDA, W, (PLASMA_desc *)(*T), Q, *LDQ); }

void PLASMA_DSYGV(PLASMA_enum *itype, PLASMA_enum *jobz, PLASMA_enum *uplo, int *N, double *A, int *LDA, double *B, int *LDB, double *W, intptr_t *T, double *Q, int *LDQ, int *INFO)
{   *INFO = PLASMA_dsygv(*itype, *jobz, *uplo, *N, A, *LDA, B, *LDB, W, (PLASMA_desc*)(*T), Q, *LDQ); }

void PLASMA_DSYGST(PLASMA_enum *itype, PLASMA_enum *uplo, int *N, double *A, int *LDA, double *B, int *LDB, int *INFO)
{   *INFO = PLASMA_dsygst(*itype, *uplo, *N, A, *LDA, B, *LDB); }

void PLASMA_DSYTRD(PLASMA_enum *jobz, PLASMA_enum *uplo, int *N, double *A, int *LDA, double *D, double *E, intptr_t *T, double *Q, int *LDQ, int *INFO)
{   *INFO = PLASMA_dsytrd(*jobz, *uplo, *N, A, *LDA, D, E, (PLASMA_desc *)(*T), Q, *LDQ); }

void PLASMA_DPOSV(PLASMA_enum *uplo, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, int *INFO)
{   *INFO = PLASMA_dposv(*uplo, *N, *NRHS, A, *LDA, B, *LDB); }

void PLASMA_DPOTRF(PLASMA_enum *uplo, int *N, double *A, int *LDA, int *INFO)
{   *INFO = PLASMA_dpotrf(*uplo, *N, A, *LDA); }

void PLASMA_DPOTRI(PLASMA_enum *uplo, int *N, double *A, int *LDA, int *INFO)
{   *INFO = PLASMA_dpotri(*uplo, *N, A, *LDA); }

void PLASMA_DPOTRS(PLASMA_enum *uplo, int *N, int *NRHS, double *A, int *LDA, double *B, int* LDB, int * INFO)
{   *INFO = PLASMA_dpotrs(*uplo, *N, *NRHS, A, *LDA, B, *LDB); }

void PLASMA_DTRSMPL(int *N, int *NRHS, double *A, int *LDA, double **LH, int **IPIVH, double *B, int *LDB, int *INFO)
{   *INFO = PLASMA_dtrsmpl(*N, *NRHS, A, *LDA, *LH, *IPIVH, B, *LDB); }

void PLASMA_DORGLQ(int *M, int *N, int *K, double *A, int *LDA, double **T, double *B, int *LDB, int *INFO)
{   *INFO = PLASMA_dorglq(*M, *N, *K, A, *LDA, *T, B, *LDB); }

void PLASMA_DORGQR(int *M, int *N, int *K, double *A, int *LDA, double **T, double *B, int *LDB, int *INFO)
{   *INFO = PLASMA_dorgqr(*M, *N, *K, A, *LDA, *T, B, *LDB); }

void PLASMA_DORMLQ(PLASMA_enum *side, PLASMA_enum *trans, int *M, int *N, int *K, double *A, int *LDA, double **T, double *B, int *LDB, int *INFO)
{   *INFO = PLASMA_dormlq(*side, *trans, *M, *N, *K, A, *LDA, *T, B, *LDB); }

void PLASMA_DORMQR(PLASMA_enum *side, PLASMA_enum *trans, int *M, int *N, int *K, double *A, int *LDA, double **T, double *B, int *LDB, int *INFO)
{   *INFO = PLASMA_dormqr(*side, *trans, *M, *N, *K, A, *LDA, *T, B, *LDB); }

void PLASMA_DTRSM(PLASMA_enum *side, PLASMA_enum *uplo, PLASMA_enum *transA, PLASMA_enum *diag, int *N, int *NRHS, double *alpha, double *A, int *LDA, double *B, int *LDB, int *INFO)
{   *INFO = PLASMA_dtrsm(*side, *uplo, *transA, *diag, *N, *NRHS, *alpha, A, *LDA, B, *LDB); }

void PLASMA_DGEMM(PLASMA_enum *transA, PLASMA_enum *transB, int *M, int *N, int *K, double *alpha, double *A, int *LDA, double *B, int *LDB, double *beta, double *C, int *LDC, int *INFO)
{   *INFO = PLASMA_dgemm(*transA, *transB, *M, *N, *K, *alpha, A, *LDA, B, *LDB, *beta, C, *LDC); }

void PLASMA_DSYMM(PLASMA_enum *side, PLASMA_enum *uplo, int *M, int *N, double *alpha, double *A, int *LDA, double *B, int *LDB, double *beta, double *C, int *LDC, int *INFO)
{   *INFO = PLASMA_dsymm(*side, *uplo, *M, *N, *alpha, A, *LDA, B, *LDB, *beta, C, *LDC); }

void PLASMA_DSYRK(PLASMA_enum *uplo, PLASMA_enum *trans, int *N, int *K, double *alpha, double *A, int *LDA, double *beta, double *C, int *LDC, int *INFO)
{   *INFO = PLASMA_dsyrk(*uplo, *trans, *N, *K, *alpha, A, *LDA, *beta, C, *LDC); }

#ifdef COMPLEX
void PLASMA_DSYMM(PLASMA_enum *side, PLASMA_enum *uplo, int *M, int *N, double *alpha, double *A, int *LDA, double *B, int *LDB, double *beta, double *C, int *LDC, int *INFO)
{   *INFO = PLASMA_dsymm(*side, *uplo, *M, *N, *alpha, A, *LDA, B, *LDB, *beta, C, *LDC); }

void PLASMA_DSYRK(PLASMA_enum *uplo, PLASMA_enum *trans, int *N, int *K, double *alpha, double *A, int *LDA, double *beta, double *C, int *LDC, int *INFO)
{   *INFO = PLASMA_dsyrk(*uplo, *trans, *N, *K, *alpha, A, *LDA, *beta, C, *LDC); }
#endif

/***************************************************************************//**
 *  FORTRAN API - math functions (native interface)
 **/
void PLASMA_DGEBRD_TILE(PLASMA_enum *jobu, PLASMA_enum *jobvt, intptr_t *A, double *D, double *E, intptr_t *U, intptr_t *VT, intptr_t *T, int *INFO)
{   *INFO = PLASMA_dgebrd_Tile(*jobu, *jobvt, (PLASMA_desc *)(*A), D, E, (PLASMA_desc *)(*U), (PLASMA_desc *)(*VT), (PLASMA_desc *)(*T)); }

void PLASMA_DGELQF_TILE(intptr_t *A, intptr_t *T, int *INFO)
{   *INFO = PLASMA_dgelqf_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*T)); }

void PLASMA_DGELQS_TILE(intptr_t *A, intptr_t *B, intptr_t *T, int *INFO)
{   *INFO = PLASMA_dgelqs_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*B), (PLASMA_desc *)(*T)); }

void PLASMA_DGELS_TILE(PLASMA_enum *trans, intptr_t *A, intptr_t *B, intptr_t *T, int *INFO)
{   *INFO = PLASMA_dgels_Tile(*trans, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B), (PLASMA_desc *)(*T)); }

void PLASMA_DGEQRF_TILE(intptr_t *A, intptr_t *T, int *INFO)
{   *INFO = PLASMA_dgeqrf_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*T)); }

void PLASMA_DGEQRS_TILE(intptr_t *A, intptr_t *B, intptr_t *T, int *INFO)
{   *INFO = PLASMA_dgeqrs_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*B), (PLASMA_desc *)(*T)); }

void PLASMA_DGESV_TILE(intptr_t *A, int *IPIV, intptr_t *B, int *INFO)
{   *INFO = PLASMA_dgesv_Tile((PLASMA_desc *)(*A), IPIV, (PLASMA_desc *)(*B)); }

void PLASMA_DGESVD_TILE(PLASMA_enum *jobu, PLASMA_enum *jobvt, intptr_t *A, double *S, intptr_t *U, intptr_t *VT, intptr_t *T, int *INFO)
{   *INFO = PLASMA_dgesvd_Tile(*jobu, *jobvt, (PLASMA_desc *)(*A), S, (PLASMA_desc *)(*U), (PLASMA_desc *)(*VT), (PLASMA_desc *)(*T)); }

void PLASMA_DGETRF_TILE(intptr_t *A, int *IPIV, int *INFO)
{   *INFO = PLASMA_dgetrf_Tile((PLASMA_desc *)(*A), IPIV); }

void PLASMA_DGETRS_TILE(PLASMA_enum *trans, intptr_t *A, int *IPIV, intptr_t *B, int *INFO)
{   *INFO = PLASMA_dgetrs_Tile(*trans, (PLASMA_desc *)(*A), IPIV, (PLASMA_desc *)(*B)); }

void PLASMA_DGESV_INCPIV_TILE(intptr_t *A, intptr_t *L, int **IPIVH, intptr_t *B, int *INFO)
{   *INFO = PLASMA_dgesv_incpiv_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*L), *IPIVH, (PLASMA_desc *)(*B)); }

void PLASMA_DGETRF_INCPIV_TILE(intptr_t *A, intptr_t *L, int **IPIVH, int *INFO)
{   *INFO = PLASMA_dgetrf_incpiv_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*L), *IPIVH); }

void PLASMA_DGETRS_INCPIV_TILE(intptr_t *A, intptr_t *L, int **IPIVH, intptr_t *B, int *INFO)
{   *INFO = PLASMA_dgetrs_incpiv_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*L), *IPIVH, (PLASMA_desc *)(*B)); }

void PLASMA_DSYEV_TILE(PLASMA_enum *jobz, PLASMA_enum *uplo, intptr_t *A, double *W, intptr_t *T, intptr_t *Q, int *INFO)
{   *INFO = PLASMA_dsyev_Tile(*jobz, *uplo, (PLASMA_desc *)(*A), W, (PLASMA_desc *)(*T), (PLASMA_desc *)(*Q)); }

void PLASMA_DSYGV_TILE(PLASMA_enum *itype, PLASMA_enum *jobz, PLASMA_enum *uplo, intptr_t *A, intptr_t *B, double *W, intptr_t *T, intptr_t *Q, int *INFO)
{   *INFO = PLASMA_dsygv_Tile(*itype, *jobz, *uplo, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B), W, (PLASMA_desc *)(*T), (PLASMA_desc *)(*Q)); }

void PLASMA_DSYGST_TILE(PLASMA_enum *itype, PLASMA_enum *uplo, intptr_t *A, intptr_t *B, int *INFO)
{   *INFO = PLASMA_dsygst_Tile(*itype, *uplo, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B)); }

void PLASMA_DSYTRD_TILE(PLASMA_enum *jobz, PLASMA_enum *uplo, intptr_t *A, double *D, double *E, intptr_t *T, intptr_t *Q, int *INFO)
{   *INFO = PLASMA_dsytrd_Tile(*jobz, *uplo, (PLASMA_desc *)(*A), D, E, (PLASMA_desc *)(*T), (PLASMA_desc *)(*Q)); }

void PLASMA_DPOSV_TILE(PLASMA_enum *uplo, intptr_t *A, intptr_t *B, int *INFO)
{   *INFO = PLASMA_dposv_Tile(*uplo, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B)); }

void PLASMA_DPOTRF_TILE(PLASMA_enum *uplo, intptr_t *A, int *INFO)
{   *INFO = PLASMA_dpotrf_Tile(*uplo, (PLASMA_desc *)(*A)); }

void PLASMA_DPOTRS_TILE(PLASMA_enum *uplo, intptr_t *A, intptr_t *B, int *INFO)
{   *INFO = PLASMA_dpotrs_Tile(*uplo, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B)); }

void PLASMA_DTRSMPL_TILE(intptr_t *A, intptr_t *L, int **IPIVH, intptr_t *B, int *INFO)
{   *INFO = PLASMA_dtrsmpl_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*L), *IPIVH, (PLASMA_desc *)(*B)); }

void PLASMA_DORGLQ_TILE(intptr_t *A, intptr_t *T, intptr_t *B, int *INFO)
{   *INFO = PLASMA_dorglq_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*T), (PLASMA_desc *)(*B)); }

void PLASMA_DORGQR_TILE(intptr_t *A, intptr_t *T, intptr_t *B, int *INFO)
{   *INFO = PLASMA_dorgqr_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*T), (PLASMA_desc *)(*B)); }

void PLASMA_DORMLQ_TILE(PLASMA_enum *side, PLASMA_enum *trans, intptr_t *A, intptr_t *T, intptr_t *B, int *INFO)
{   *INFO = PLASMA_dormlq_Tile(*side, *trans, (PLASMA_desc *)(*A), (PLASMA_desc *)(*T), (PLASMA_desc *)(*B)); }

void PLASMA_DORMQR_TILE(PLASMA_enum *side, PLASMA_enum *trans, intptr_t *A, intptr_t *T, intptr_t *B, int *INFO)
{   *INFO = PLASMA_dormqr_Tile(*side, *trans, (PLASMA_desc *)(*A), (PLASMA_desc *)(*T), (PLASMA_desc *)(*B)); }

void PLASMA_DTRSM_TILE(PLASMA_enum *side, PLASMA_enum *uplo, PLASMA_enum *transA, PLASMA_enum *diag, double *alpha, intptr_t *A, intptr_t *B, int *INFO)
{   *INFO = PLASMA_dtrsm_Tile(*side, *uplo, *transA, *diag, *alpha, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B)); }

void PLASMA_DGEMM_TILE(PLASMA_enum *transA, PLASMA_enum *transB, int *alpha, intptr_t *A, intptr_t *B, int *beta, intptr_t *C, int *INFO)
{   *INFO = PLASMA_dgemm_Tile(*transA, *transB, *alpha, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B), *beta, (PLASMA_desc *)(*C)); }

  void PLASMA_DSYMM_TILE(PLASMA_enum *side, PLASMA_enum *uplo, double *alpha, intptr_t *A, intptr_t *B, double *beta, intptr_t *C, int *INFO)
  {   *INFO = PLASMA_dsymm_Tile(*side, *uplo, *alpha, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B), *beta, (PLASMA_desc *)(*C)); }

void PLASMA_DSYRK_TILE(PLASMA_enum *uplo, PLASMA_enum *trans, double *alpha, intptr_t *A, double *beta, intptr_t *C, int *INFO)
{   *INFO = PLASMA_dsyrk_Tile(*uplo, *trans, *alpha, (PLASMA_desc *)(*A), *beta, (PLASMA_desc *)(*C)); }

#ifdef COMPLEX
void PLASMA_DSYMM_TILE(PLASMA_enum *side, PLASMA_enum *uplo, double *alpha, intptr_t *A, intptr_t *B, double *beta, intptr_t *C, int *INFO)
{   *INFO = PLASMA_dsymm_Tile(*side, *uplo, *alpha, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B), *beta, (PLASMA_desc *)(*C)); }

void PLASMA_DSYRK_TILE(PLASMA_enum *uplo, PLASMA_enum *trans, double *alpha, intptr_t *A, double *beta, intptr_t *C, int *INFO)
{   *INFO = PLASMA_dsyrk_Tile(*uplo, *trans, *alpha, (PLASMA_desc *)(*A), *beta, (PLASMA_desc *)(*C)); }
#endif

/***************************************************************************//**
 *  FORTRAN API - workspace allocation
 **/
void PLASMA_ALLOC_WORKSPACE_DGEBRD(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_dgebrd(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_DGELQF(int *M, int *N, double **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_dgelqf(*M, *N, T); }

void PLASMA_ALLOC_WORKSPACE_DGELS(int *M, int *N, double **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_dgels(*M, *N, T); }

void PLASMA_ALLOC_WORKSPACE_DGEQRF(int *M, int *N, double **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_dgeqrf(*M, *N, T); }

void PLASMA_ALLOC_WORKSPACE_DGESV_INCPIV(int *N, double **L, int **IPIV, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_dgesv_incpiv(*N, L, IPIV); }

void PLASMA_ALLOC_WORKSPACE_DGESVD(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_dgesvd(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_DGETRF_INCPIV(int *M, int *N, double **L, int **IPIV, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_dgetrf_incpiv(*M, *N, L, IPIV); }

void PLASMA_ALLOC_WORKSPACE_DSYEV(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_dsyev(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_DSYGV(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_dsygv(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_DSYTRD(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_dsytrd(*M, *N, (PLASMA_desc **)T); }


/***************************************************************************//**
 *  FORTRAN API - tiled workspace allocation
 **/
void PLASMA_ALLOC_WORKSPACE_DGELQF_TILE(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_dgelqf_Tile(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_DGELS_TILE(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_dgels_Tile(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_DGEQRF_TILE(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_dgeqrf_Tile(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_DGESV_INCPIV_TILE(int *N, intptr_t **L, int **IPIV, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_dgesv_incpiv_Tile(*N, (PLASMA_desc **)L, IPIV); }

void PLASMA_ALLOC_WORKSPACE_DGETRF_INCPIV_TILE(int *N, intptr_t **L, int **IPIV, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_dgetrf_incpiv_Tile(*N, (PLASMA_desc **)L, IPIV); }

/***************************************************************************//**
 *  FORTRAN API - conversion from LAPACK F77 matrix layout to tile layout
 **/
void PLASMA_DLAPACK_TO_TILE(double **Af77, int *LDA, intptr_t *A, int *INFO)
{   *INFO = PLASMA_dLapack_to_Tile( *Af77, *LDA, (PLASMA_desc *)(*A) ); }

void PLASMA_DTILE_TO_LAPACK(intptr_t *A, double **Af77, int *LDA, int *INFO)
{   *INFO = PLASMA_dTile_to_Lapack( (PLASMA_desc *)(*A), *Af77, *LDA ); }

#ifdef __cplusplus
}
#endif
