/**
 *
 * @file plasma_zf77.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Bilel Hadri
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include <stdlib.h>
#include "common.h"
#undef REAL
#define COMPLEX

/*
 * Lapack interface
 */
#define PLASMA_ZGEBRD        PLASMA_FNAME(zgebrd,  ZGEBRD )
#define PLASMA_ZGEEV         PLASMA_FNAME(zgeev,   ZGEEV  )
#define PLASMA_ZGEHRD        PLASMA_FNAME(zgehrd,  ZGEHRD )
#define PLASMA_ZGELQF        PLASMA_FNAME(zgelqf,  ZGELQF )
#define PLASMA_ZGELQS        PLASMA_FNAME(zgelqs,  ZGELQS )
#define PLASMA_ZGELS         PLASMA_FNAME(zgels,   ZGELS  )
#define PLASMA_ZGEQRF        PLASMA_FNAME(zgeqrf,  ZGEQRF )
#define PLASMA_ZGEQRS        PLASMA_FNAME(zgeqrs,  ZGEQRS )
#define PLASMA_ZGESV         PLASMA_FNAME(zgesv,   ZGESV  )
#define PLASMA_ZGESVD        PLASMA_FNAME(zgesvd,  ZGESVD )
#define PLASMA_ZGETRF        PLASMA_FNAME(zgetrf,  ZGETRF )
#define PLASMA_ZGETRS        PLASMA_FNAME(zgetrs,  ZGETRS )
#define PLASMA_ZGESV_INCPIV  PLASMA_FNAME(zgesv_incpiv,  ZGESV_INCPIV  )
#define PLASMA_ZGETRF_INCPIV PLASMA_FNAME(zgetrf_incpiv, ZGETRF_INCPIV )
#define PLASMA_ZGETRS_INCPIV PLASMA_FNAME(zgetrs_incpiv, ZGETRS_INCPIV )
#define PLASMA_ZHEEV         PLASMA_FNAME(zheev,   ZHEEV  )
#define PLASMA_ZHEGV         PLASMA_FNAME(zhegv,   ZHEGV  )
#define PLASMA_ZHEGST        PLASMA_FNAME(zhegst,  ZHEGST )
#define PLASMA_ZHETRD        PLASMA_FNAME(zhetrd,  ZHETRD )
#define PLASMA_ZPOSV         PLASMA_FNAME(zposv,   ZPOSV  )
#define PLASMA_ZPOTRF        PLASMA_FNAME(zpotrf,  ZPOTRF )
#define PLASMA_ZPOTRI        PLASMA_FNAME(zpotri,  ZPOTRI )
#define PLASMA_ZPOTRS        PLASMA_FNAME(zpotrs,  ZPOTRS )
#define PLASMA_ZTRSMPL       PLASMA_FNAME(ztrsmpl, ZTRSMPL)
#define PLASMA_ZUNGBR        PLASMA_FNAME(zungbr,  ZUNGBR )
#define PLASMA_ZUNGHR        PLASMA_FNAME(zunghr,  ZUNGHR )
#define PLASMA_ZUNGLQ        PLASMA_FNAME(zunglq,  ZUNGLQ )
#define PLASMA_ZUNGQR        PLASMA_FNAME(zungqr,  ZUNGQR )
#define PLASMA_ZUNGTR        PLASMA_FNAME(zungtr,  ZUNGTR )
#define PLASMA_ZUNMLQ        PLASMA_FNAME(zunmlq,  ZUNMLQ )
#define PLASMA_ZUNMQR        PLASMA_FNAME(zunmqr,  ZUNMQR )
#define PLASMA_ZTRSM         PLASMA_FNAME(ztrsm,   ZTRSM  )
#define PLASMA_ZGEMM         PLASMA_FNAME(zgemm,   ZGEMM  )
#define PLASMA_ZSYMM         PLASMA_FNAME(zsymm,   ZSYMM  )
#define PLASMA_ZSYRK         PLASMA_FNAME(zsyrk,   ZSYRK  )
#ifdef COMPLEX
#define PLASMA_ZHEMM         PLASMA_FNAME(zhemm,   ZHEMM  )
#define PLASMA_ZHERK         PLASMA_FNAME(zherk,   ZHERK  )
#endif

/*
 * Tile interface
 */
#define PLASMA_ZGEBRD_TILE        PLASMA_TILE_FNAME(zgebrd,  ZGEBRD )
#define PLASMA_ZGEEV_TILE         PLASMA_TILE_FNAME(zgeev,   ZGEEV  )
#define PLASMA_ZGEHRD_TILE        PLASMA_TILE_FNAME(zgehrd,  ZGEHRD )
#define PLASMA_ZGELQF_TILE        PLASMA_TILE_FNAME(zgelqf,  ZGELQF )
#define PLASMA_ZGELQS_TILE        PLASMA_TILE_FNAME(zgelqs,  ZGELQS )
#define PLASMA_ZGELS_TILE         PLASMA_TILE_FNAME(zgels,   ZGELS  )
#define PLASMA_ZGEQRF_TILE        PLASMA_TILE_FNAME(zgeqrf,  ZGEQRF )
#define PLASMA_ZGEQRS_TILE        PLASMA_TILE_FNAME(zgeqrs,  ZGEQRS )
#define PLASMA_ZGESV_TILE         PLASMA_TILE_FNAME(zgesv,   ZGESV  )
#define PLASMA_ZGESVD_TILE        PLASMA_TILE_FNAME(zgesvd,  ZGESVD )
#define PLASMA_ZGETRF_TILE        PLASMA_TILE_FNAME(zgetrf,  ZGETRF )
#define PLASMA_ZGETRS_TILE        PLASMA_TILE_FNAME(zgetrs,  ZGETRS ) 
#define PLASMA_ZGESV_INCPIV_TILE  PLASMA_TILE_FNAME(zgesv_incpiv,  ZGESV_INCPIV  )
#define PLASMA_ZGETRF_INCPIV_TILE PLASMA_TILE_FNAME(zgetrf_incpiv, ZGETRF_INCPIV )
#define PLASMA_ZGETRS_INCPIV_TILE PLASMA_TILE_FNAME(zgetrs_incpiv, ZGETRS_INCPIV ) 
#define PLASMA_ZHEEV_TILE         PLASMA_TILE_FNAME(zheev,   ZHEEV  )
#define PLASMA_ZHEGV_TILE         PLASMA_TILE_FNAME(zhegv,   ZHEGV  )
#define PLASMA_ZHEGST_TILE        PLASMA_TILE_FNAME(zhegst,  ZHEGST )
#define PLASMA_ZHETRD_TILE        PLASMA_TILE_FNAME(zhetrd,  ZHETRD )
#define PLASMA_ZPOSV_TILE         PLASMA_TILE_FNAME(zposv,   ZPOSV  )
#define PLASMA_ZPOTRF_TILE        PLASMA_TILE_FNAME(zpotrf,  ZPOTRF )
#define PLASMA_ZPOTRI_TILE        PLASMA_TILE_FNAME(zpotri,  ZPOTRI )
#define PLASMA_ZPOTRS_TILE        PLASMA_TILE_FNAME(zpotrs,  ZPOTRS )
#define PLASMA_ZTRSM_TILE         PLASMA_TILE_FNAME(ztrsm,   ZTRSM  )
#define PLASMA_ZTRSMPL_TILE       PLASMA_TILE_FNAME(ztrsmpl, ZTRSMPL)
#define PLASMA_ZUNGBR_TILE        PLASMA_TILE_FNAME(zungbr,  ZUNGBR )
#define PLASMA_ZUNGHR_TILE        PLASMA_TILE_FNAME(zunghr,  ZUNGHR )
#define PLASMA_ZUNGLQ_TILE        PLASMA_TILE_FNAME(zunglq,  ZUNGLQ )
#define PLASMA_ZUNGQR_TILE        PLASMA_TILE_FNAME(zungqr,  ZUNGQR )
#define PLASMA_ZUNGTR_TILE        PLASMA_TILE_FNAME(zungtr,  ZUNGTR )
#define PLASMA_ZUNMLQ_TILE        PLASMA_TILE_FNAME(zunmlq,  ZUNMLQ )
#define PLASMA_ZUNMQR_TILE        PLASMA_TILE_FNAME(zunmqr,  ZUNMQR )
#define PLASMA_ZGEMM_TILE         PLASMA_TILE_FNAME(zgemm,   ZGEMM  )
#define PLASMA_ZSYMM_TILE         PLASMA_TILE_FNAME(zsymm,   ZSYMM  )
#define PLASMA_ZSYRK_TILE         PLASMA_TILE_FNAME(zsyrk,   ZSYRK  )
#ifdef COMPLEX                    
#define PLASMA_ZHEMM_TILE         PLASMA_TILE_FNAME(zhemm,   ZHEMM  )
#define PLASMA_ZHERK_TILE         PLASMA_TILE_FNAME(zherk,   ZHERK  )
#endif

/*
 * Workspaces
 */
#define PLASMA_ALLOC_WORKSPACE_ZGEBRD        PLASMA_WS_FNAME(zgehrd, ZGEBRD)
#define PLASMA_ALLOC_WORKSPACE_ZGEEV         PLASMA_WS_FNAME(zgeev,  ZGEEV )
#define PLASMA_ALLOC_WORKSPACE_ZGEHRD        PLASMA_WS_FNAME(zgehrd, ZGEHRD)
#define PLASMA_ALLOC_WORKSPACE_ZGELQF        PLASMA_WS_FNAME(zgelqf, ZGELQF) 
#define PLASMA_ALLOC_WORKSPACE_ZGELS         PLASMA_WS_FNAME(zgels,  ZGELS )  
#define PLASMA_ALLOC_WORKSPACE_ZGEQRF        PLASMA_WS_FNAME(zgeqrf, ZGEQRF) 
#define PLASMA_ALLOC_WORKSPACE_ZGESV_INCPIV  PLASMA_WS_FNAME(zgesv_incpiv,  ZGESV_INCPIV ) 
#define PLASMA_ALLOC_WORKSPACE_ZGETRF_INCPIV PLASMA_WS_FNAME(zgetrf_incpiv, ZGETRF_INCPIV) 
#define PLASMA_ALLOC_WORKSPACE_ZGESVD        PLASMA_WS_FNAME(zgesvd, ZGESVD)
#define PLASMA_ALLOC_WORKSPACE_ZHEEV         PLASMA_WS_FNAME(zheev,  ZHEEV )
#define PLASMA_ALLOC_WORKSPACE_ZHEGV         PLASMA_WS_FNAME(zhegv,  ZHEGV )
#define PLASMA_ALLOC_WORKSPACE_ZHETRD        PLASMA_WS_FNAME(zhetrd, ZHETRD)

/*
 * Tile Workspaces
 */
#define PLASMA_ALLOC_WORKSPACE_ZGELQF_TILE        PLASMA_WST_FNAME(zgelqf,        ZGELQF       )
#define PLASMA_ALLOC_WORKSPACE_ZGELS_TILE         PLASMA_WST_FNAME(zgels,         ZGELS        )
#define PLASMA_ALLOC_WORKSPACE_ZGEQRF_TILE        PLASMA_WST_FNAME(zgeqrf,        ZGEQRF       )
#define PLASMA_ALLOC_WORKSPACE_ZGESV_INCPIV_TILE  PLASMA_WST_FNAME(zgesv_incpiv,  ZGESV_INCPIV )
#define PLASMA_ALLOC_WORKSPACE_ZGETRF_INCPIV_TILE PLASMA_WST_FNAME(zgetrf_incpiv, ZGETRF_INCPIV)

#define PLASMA_ZLAPACK_TO_TILE   PLASMA_FNAME(zlapack_to_tile, ZLAPACK_TO_TILE)
#define PLASMA_ZTILE_TO_LAPACK   PLASMA_FNAME(ztile_to_lapack, ZTILE_TO_LAPACK)

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  FORTRAN API - math functions (simple interface)
 **/
void PLASMA_ZGEBRD(PLASMA_enum *jobu, PLASMA_enum *jobvt, int *M, int *N, PLASMA_Complex64_t *A, int *LDA, double *D, double *E, PLASMA_Complex64_t *U, int *LDU, PLASMA_Complex64_t *VT, int *LDVT, intptr_t *descT, int *INFO)
{   *INFO = PLASMA_zgebrd(*jobu, *jobvt, *M, *N, A, *LDA, D, E, U, *LDU, VT, *LDVT, (PLASMA_desc *)(*descT)); }

void PLASMA_ZGELQF(int *M, int *N, PLASMA_Complex64_t *A, int *LDA, PLASMA_Complex64_t **T, int *INFO)
{   *INFO = PLASMA_zgelqf(*M, *N, A, *LDA, *T); }

void PLASMA_ZGELQS(int *M, int *N, int *NRHS, PLASMA_Complex64_t *A, int *LDA, PLASMA_Complex64_t **T, PLASMA_Complex64_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_zgelqs(*M, *N, *NRHS, A, *LDA, *T, B, *LDB); }

void PLASMA_ZGELS(PLASMA_enum *trans, int *M, int *N, int *NRHS, PLASMA_Complex64_t *A, int *LDA, PLASMA_Complex64_t **T, PLASMA_Complex64_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_zgels(*trans, *M, *N, *NRHS, A, *LDA, *T, B, *LDB); }

void PLASMA_ZGEQRF(int *M, int *N, PLASMA_Complex64_t *A, int *LDA, PLASMA_Complex64_t **T, int *INFO)
{   *INFO = PLASMA_zgeqrf(*M, *N, A, *LDA, *T); }

void PLASMA_ZGEQRS(int *M, int *N, int *NRHS, PLASMA_Complex64_t *A, int *LDA, PLASMA_Complex64_t **T, PLASMA_Complex64_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_zgeqrs(*M, *N, *NRHS, A, *LDA, *T, B, *LDB); }

void PLASMA_ZGESV(int *N, int *NRHS, PLASMA_Complex64_t *A, int *LDA, int *IPIV, PLASMA_Complex64_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_zgesv(*N, *NRHS, A, *LDA, IPIV, B, *LDB); }

void PLASMA_ZGESVD(PLASMA_enum *jobu, PLASMA_enum *jobvt, int *M, int *N, PLASMA_Complex64_t *A, int *LDA, double *S, PLASMA_Complex64_t *U, int *LDU, PLASMA_Complex64_t *VT, int *LDVT, intptr_t *T, int *INFO)
{   *INFO = PLASMA_zgesvd(*jobu, *jobvt, *M, *N, A, *LDA, S, U, *LDU, VT, *LDVT, (PLASMA_desc *)(*T)); }

void PLASMA_ZGETRF(int *M, int *N, PLASMA_Complex64_t *A, int *LDA, int *IPIV, int *INFO)
{   *INFO = PLASMA_zgetrf(*M, *N, A, *LDA, IPIV); }

void PLASMA_ZGETRS(PLASMA_enum *trans, int *N, int *NRHS, PLASMA_Complex64_t *A, int *LDA, int *IPIV, PLASMA_Complex64_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_zgetrs(*trans, *N, *NRHS, A, *LDA, IPIV, B, *LDB); }

void PLASMA_ZGESV_INCPIV(int *N, int *NRHS, PLASMA_Complex64_t *A, int *LDA, PLASMA_Complex64_t **LH, int **IPIVH, PLASMA_Complex64_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_zgesv_incpiv(*N, *NRHS, A, *LDA, *LH, *IPIVH, B, *LDB); }

void PLASMA_ZGETRF_INCPIV(int *M, int *N, PLASMA_Complex64_t *A, int *LDA, PLASMA_Complex64_t **LH, int **IPIVH, int *INFO)
{   *INFO = PLASMA_zgetrf_incpiv(*M, *N, A, *LDA, *LH, *IPIVH); }

void PLASMA_ZGETRS_INCPIV(PLASMA_enum *uplo, int *N, int *NRHS, PLASMA_Complex64_t *A, int *LDA, PLASMA_Complex64_t **LH, int **IPIVH, PLASMA_Complex64_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_zgetrs_incpiv(*uplo, *N, *NRHS, A, *LDA, *LH, *IPIVH, B, *LDB); }

void PLASMA_ZHEEV(PLASMA_enum *jobz, PLASMA_enum *uplo, int *N, PLASMA_Complex64_t *A, int *LDA, double *W, intptr_t *T, PLASMA_Complex64_t *Q, int *LDQ, int *INFO)
{   *INFO = PLASMA_zheev(*jobz, *uplo, *N, A, *LDA, W, (PLASMA_desc *)(*T), Q, *LDQ); }

void PLASMA_ZHEGV(PLASMA_enum *itype, PLASMA_enum *jobz, PLASMA_enum *uplo, int *N, PLASMA_Complex64_t *A, int *LDA, PLASMA_Complex64_t *B, int *LDB, double *W, intptr_t *T, PLASMA_Complex64_t *Q, int *LDQ, int *INFO)
{   *INFO = PLASMA_zhegv(*itype, *jobz, *uplo, *N, A, *LDA, B, *LDB, W, (PLASMA_desc*)(*T), Q, *LDQ); }

void PLASMA_ZHEGST(PLASMA_enum *itype, PLASMA_enum *uplo, int *N, PLASMA_Complex64_t *A, int *LDA, PLASMA_Complex64_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_zhegst(*itype, *uplo, *N, A, *LDA, B, *LDB); }

void PLASMA_ZHETRD(PLASMA_enum *jobz, PLASMA_enum *uplo, int *N, PLASMA_Complex64_t *A, int *LDA, double *D, double *E, intptr_t *T, PLASMA_Complex64_t *Q, int *LDQ, int *INFO)
{   *INFO = PLASMA_zhetrd(*jobz, *uplo, *N, A, *LDA, D, E, (PLASMA_desc *)(*T), Q, *LDQ); }

void PLASMA_ZPOSV(PLASMA_enum *uplo, int *N, int *NRHS, PLASMA_Complex64_t *A, int *LDA, PLASMA_Complex64_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_zposv(*uplo, *N, *NRHS, A, *LDA, B, *LDB); }

void PLASMA_ZPOTRF(PLASMA_enum *uplo, int *N, PLASMA_Complex64_t *A, int *LDA, int *INFO)
{   *INFO = PLASMA_zpotrf(*uplo, *N, A, *LDA); }

void PLASMA_ZPOTRI(PLASMA_enum *uplo, int *N, PLASMA_Complex64_t *A, int *LDA, int *INFO)
{   *INFO = PLASMA_zpotri(*uplo, *N, A, *LDA); }

void PLASMA_ZPOTRS(PLASMA_enum *uplo, int *N, int *NRHS, PLASMA_Complex64_t *A, int *LDA, PLASMA_Complex64_t *B, int* LDB, int * INFO)
{   *INFO = PLASMA_zpotrs(*uplo, *N, *NRHS, A, *LDA, B, *LDB); }

void PLASMA_ZTRSMPL(int *N, int *NRHS, PLASMA_Complex64_t *A, int *LDA, PLASMA_Complex64_t **LH, int **IPIVH, PLASMA_Complex64_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_ztrsmpl(*N, *NRHS, A, *LDA, *LH, *IPIVH, B, *LDB); }

void PLASMA_ZUNGLQ(int *M, int *N, int *K, PLASMA_Complex64_t *A, int *LDA, PLASMA_Complex64_t **T, PLASMA_Complex64_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_zunglq(*M, *N, *K, A, *LDA, *T, B, *LDB); }

void PLASMA_ZUNGQR(int *M, int *N, int *K, PLASMA_Complex64_t *A, int *LDA, PLASMA_Complex64_t **T, PLASMA_Complex64_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_zungqr(*M, *N, *K, A, *LDA, *T, B, *LDB); }

void PLASMA_ZUNMLQ(PLASMA_enum *side, PLASMA_enum *trans, int *M, int *N, int *K, PLASMA_Complex64_t *A, int *LDA, PLASMA_Complex64_t **T, PLASMA_Complex64_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_zunmlq(*side, *trans, *M, *N, *K, A, *LDA, *T, B, *LDB); }

void PLASMA_ZUNMQR(PLASMA_enum *side, PLASMA_enum *trans, int *M, int *N, int *K, PLASMA_Complex64_t *A, int *LDA, PLASMA_Complex64_t **T, PLASMA_Complex64_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_zunmqr(*side, *trans, *M, *N, *K, A, *LDA, *T, B, *LDB); }

void PLASMA_ZTRSM(PLASMA_enum *side, PLASMA_enum *uplo, PLASMA_enum *transA, PLASMA_enum *diag, int *N, int *NRHS, PLASMA_Complex64_t *alpha, PLASMA_Complex64_t *A, int *LDA, PLASMA_Complex64_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_ztrsm(*side, *uplo, *transA, *diag, *N, *NRHS, *alpha, A, *LDA, B, *LDB); }

void PLASMA_ZGEMM(PLASMA_enum *transA, PLASMA_enum *transB, int *M, int *N, int *K, PLASMA_Complex64_t *alpha, PLASMA_Complex64_t *A, int *LDA, PLASMA_Complex64_t *B, int *LDB, PLASMA_Complex64_t *beta, PLASMA_Complex64_t *C, int *LDC, int *INFO)
{   *INFO = PLASMA_zgemm(*transA, *transB, *M, *N, *K, *alpha, A, *LDA, B, *LDB, *beta, C, *LDC); }

void PLASMA_ZSYMM(PLASMA_enum *side, PLASMA_enum *uplo, int *M, int *N, PLASMA_Complex64_t *alpha, PLASMA_Complex64_t *A, int *LDA, PLASMA_Complex64_t *B, int *LDB, PLASMA_Complex64_t *beta, PLASMA_Complex64_t *C, int *LDC, int *INFO)
{   *INFO = PLASMA_zsymm(*side, *uplo, *M, *N, *alpha, A, *LDA, B, *LDB, *beta, C, *LDC); }

void PLASMA_ZSYRK(PLASMA_enum *uplo, PLASMA_enum *trans, int *N, int *K, PLASMA_Complex64_t *alpha, PLASMA_Complex64_t *A, int *LDA, PLASMA_Complex64_t *beta, PLASMA_Complex64_t *C, int *LDC, int *INFO)
{   *INFO = PLASMA_zsyrk(*uplo, *trans, *N, *K, *alpha, A, *LDA, *beta, C, *LDC); }

#ifdef COMPLEX
void PLASMA_ZHEMM(PLASMA_enum *side, PLASMA_enum *uplo, int *M, int *N, PLASMA_Complex64_t *alpha, PLASMA_Complex64_t *A, int *LDA, PLASMA_Complex64_t *B, int *LDB, PLASMA_Complex64_t *beta, PLASMA_Complex64_t *C, int *LDC, int *INFO)
{   *INFO = PLASMA_zhemm(*side, *uplo, *M, *N, *alpha, A, *LDA, B, *LDB, *beta, C, *LDC); }

void PLASMA_ZHERK(PLASMA_enum *uplo, PLASMA_enum *trans, int *N, int *K, PLASMA_Complex64_t *alpha, PLASMA_Complex64_t *A, int *LDA, double *beta, PLASMA_Complex64_t *C, int *LDC, int *INFO)
{   *INFO = PLASMA_zherk(*uplo, *trans, *N, *K, *alpha, A, *LDA, *beta, C, *LDC); }
#endif

/***************************************************************************//**
 *  FORTRAN API - math functions (native interface)
 **/
void PLASMA_ZGEBRD_TILE(PLASMA_enum *jobu, PLASMA_enum *jobvt, intptr_t *A, double *D, double *E, intptr_t *U, intptr_t *VT, intptr_t *T, int *INFO)
{   *INFO = PLASMA_zgebrd_Tile(*jobu, *jobvt, (PLASMA_desc *)(*A), D, E, (PLASMA_desc *)(*U), (PLASMA_desc *)(*VT), (PLASMA_desc *)(*T)); }

void PLASMA_ZGELQF_TILE(intptr_t *A, intptr_t *T, int *INFO)
{   *INFO = PLASMA_zgelqf_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*T)); }

void PLASMA_ZGELQS_TILE(intptr_t *A, intptr_t *B, intptr_t *T, int *INFO)
{   *INFO = PLASMA_zgelqs_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*B), (PLASMA_desc *)(*T)); }

void PLASMA_ZGELS_TILE(PLASMA_enum *trans, intptr_t *A, intptr_t *B, intptr_t *T, int *INFO)
{   *INFO = PLASMA_zgels_Tile(*trans, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B), (PLASMA_desc *)(*T)); }

void PLASMA_ZGEQRF_TILE(intptr_t *A, intptr_t *T, int *INFO)
{   *INFO = PLASMA_zgeqrf_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*T)); }

void PLASMA_ZGEQRS_TILE(intptr_t *A, intptr_t *B, intptr_t *T, int *INFO)
{   *INFO = PLASMA_zgeqrs_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*B), (PLASMA_desc *)(*T)); }

void PLASMA_ZGESV_TILE(intptr_t *A, int *IPIV, intptr_t *B, int *INFO)
{   *INFO = PLASMA_zgesv_Tile((PLASMA_desc *)(*A), IPIV, (PLASMA_desc *)(*B)); }

void PLASMA_ZGESVD_TILE(PLASMA_enum *jobu, PLASMA_enum *jobvt, intptr_t *A, double *S, intptr_t *U, intptr_t *VT, intptr_t *T, int *INFO)
{   *INFO = PLASMA_zgesvd_Tile(*jobu, *jobvt, (PLASMA_desc *)(*A), S, (PLASMA_desc *)(*U), (PLASMA_desc *)(*VT), (PLASMA_desc *)(*T)); }

void PLASMA_ZGETRF_TILE(intptr_t *A, int *IPIV, int *INFO)
{   *INFO = PLASMA_zgetrf_Tile((PLASMA_desc *)(*A), IPIV); }

void PLASMA_ZGETRS_TILE(PLASMA_enum *trans, intptr_t *A, int *IPIV, intptr_t *B, int *INFO)
{   *INFO = PLASMA_zgetrs_Tile(*trans, (PLASMA_desc *)(*A), IPIV, (PLASMA_desc *)(*B)); }

void PLASMA_ZGESV_INCPIV_TILE(intptr_t *A, intptr_t *L, int **IPIVH, intptr_t *B, int *INFO)
{   *INFO = PLASMA_zgesv_incpiv_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*L), *IPIVH, (PLASMA_desc *)(*B)); }

void PLASMA_ZGETRF_INCPIV_TILE(intptr_t *A, intptr_t *L, int **IPIVH, int *INFO)
{   *INFO = PLASMA_zgetrf_incpiv_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*L), *IPIVH); }

void PLASMA_ZGETRS_INCPIV_TILE(intptr_t *A, intptr_t *L, int **IPIVH, intptr_t *B, int *INFO)
{   *INFO = PLASMA_zgetrs_incpiv_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*L), *IPIVH, (PLASMA_desc *)(*B)); }

void PLASMA_ZHEEV_TILE(PLASMA_enum *jobz, PLASMA_enum *uplo, intptr_t *A, double *W, intptr_t *T, intptr_t *Q, int *INFO)
{   *INFO = PLASMA_zheev_Tile(*jobz, *uplo, (PLASMA_desc *)(*A), W, (PLASMA_desc *)(*T), (PLASMA_desc *)(*Q)); }

void PLASMA_ZHEGV_TILE(PLASMA_enum *itype, PLASMA_enum *jobz, PLASMA_enum *uplo, intptr_t *A, intptr_t *B, double *W, intptr_t *T, intptr_t *Q, int *INFO)
{   *INFO = PLASMA_zhegv_Tile(*itype, *jobz, *uplo, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B), W, (PLASMA_desc *)(*T), (PLASMA_desc *)(*Q)); }

void PLASMA_ZHEGST_TILE(PLASMA_enum *itype, PLASMA_enum *uplo, intptr_t *A, intptr_t *B, int *INFO)
{   *INFO = PLASMA_zhegst_Tile(*itype, *uplo, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B)); }

void PLASMA_ZHETRD_TILE(PLASMA_enum *jobz, PLASMA_enum *uplo, intptr_t *A, double *D, double *E, intptr_t *T, intptr_t *Q, int *INFO)
{   *INFO = PLASMA_zhetrd_Tile(*jobz, *uplo, (PLASMA_desc *)(*A), D, E, (PLASMA_desc *)(*T), (PLASMA_desc *)(*Q)); }

void PLASMA_ZPOSV_TILE(PLASMA_enum *uplo, intptr_t *A, intptr_t *B, int *INFO)
{   *INFO = PLASMA_zposv_Tile(*uplo, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B)); }

void PLASMA_ZPOTRF_TILE(PLASMA_enum *uplo, intptr_t *A, int *INFO)
{   *INFO = PLASMA_zpotrf_Tile(*uplo, (PLASMA_desc *)(*A)); }

void PLASMA_ZPOTRS_TILE(PLASMA_enum *uplo, intptr_t *A, intptr_t *B, int *INFO)
{   *INFO = PLASMA_zpotrs_Tile(*uplo, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B)); }

void PLASMA_ZTRSMPL_TILE(intptr_t *A, intptr_t *L, int **IPIVH, intptr_t *B, int *INFO)
{   *INFO = PLASMA_ztrsmpl_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*L), *IPIVH, (PLASMA_desc *)(*B)); }

void PLASMA_ZUNGLQ_TILE(intptr_t *A, intptr_t *T, intptr_t *B, int *INFO)
{   *INFO = PLASMA_zunglq_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*T), (PLASMA_desc *)(*B)); }

void PLASMA_ZUNGQR_TILE(intptr_t *A, intptr_t *T, intptr_t *B, int *INFO)
{   *INFO = PLASMA_zungqr_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*T), (PLASMA_desc *)(*B)); }

void PLASMA_ZUNMLQ_TILE(PLASMA_enum *side, PLASMA_enum *trans, intptr_t *A, intptr_t *T, intptr_t *B, int *INFO)
{   *INFO = PLASMA_zunmlq_Tile(*side, *trans, (PLASMA_desc *)(*A), (PLASMA_desc *)(*T), (PLASMA_desc *)(*B)); }

void PLASMA_ZUNMQR_TILE(PLASMA_enum *side, PLASMA_enum *trans, intptr_t *A, intptr_t *T, intptr_t *B, int *INFO)
{   *INFO = PLASMA_zunmqr_Tile(*side, *trans, (PLASMA_desc *)(*A), (PLASMA_desc *)(*T), (PLASMA_desc *)(*B)); }

void PLASMA_ZTRSM_TILE(PLASMA_enum *side, PLASMA_enum *uplo, PLASMA_enum *transA, PLASMA_enum *diag, PLASMA_Complex64_t *alpha, intptr_t *A, intptr_t *B, int *INFO)
{   *INFO = PLASMA_ztrsm_Tile(*side, *uplo, *transA, *diag, *alpha, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B)); }

void PLASMA_ZGEMM_TILE(PLASMA_enum *transA, PLASMA_enum *transB, int *alpha, intptr_t *A, intptr_t *B, int *beta, intptr_t *C, int *INFO)
{   *INFO = PLASMA_zgemm_Tile(*transA, *transB, *alpha, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B), *beta, (PLASMA_desc *)(*C)); }

  void PLASMA_ZSYMM_TILE(PLASMA_enum *side, PLASMA_enum *uplo, PLASMA_Complex64_t *alpha, intptr_t *A, intptr_t *B, PLASMA_Complex64_t *beta, intptr_t *C, int *INFO)
  {   *INFO = PLASMA_zsymm_Tile(*side, *uplo, *alpha, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B), *beta, (PLASMA_desc *)(*C)); }

void PLASMA_ZSYRK_TILE(PLASMA_enum *uplo, PLASMA_enum *trans, PLASMA_Complex64_t *alpha, intptr_t *A, PLASMA_Complex64_t *beta, intptr_t *C, int *INFO)
{   *INFO = PLASMA_zsyrk_Tile(*uplo, *trans, *alpha, (PLASMA_desc *)(*A), *beta, (PLASMA_desc *)(*C)); }

#ifdef COMPLEX
void PLASMA_ZHEMM_TILE(PLASMA_enum *side, PLASMA_enum *uplo, PLASMA_Complex64_t *alpha, intptr_t *A, intptr_t *B, PLASMA_Complex64_t *beta, intptr_t *C, int *INFO)
{   *INFO = PLASMA_zhemm_Tile(*side, *uplo, *alpha, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B), *beta, (PLASMA_desc *)(*C)); }

void PLASMA_ZHERK_TILE(PLASMA_enum *uplo, PLASMA_enum *trans, PLASMA_Complex64_t *alpha, intptr_t *A, double *beta, intptr_t *C, int *INFO)
{   *INFO = PLASMA_zherk_Tile(*uplo, *trans, *alpha, (PLASMA_desc *)(*A), *beta, (PLASMA_desc *)(*C)); }
#endif

/***************************************************************************//**
 *  FORTRAN API - workspace allocation
 **/
void PLASMA_ALLOC_WORKSPACE_ZGEBRD(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_zgebrd(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_ZGELQF(int *M, int *N, PLASMA_Complex64_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_zgelqf(*M, *N, T); }

void PLASMA_ALLOC_WORKSPACE_ZGELS(int *M, int *N, PLASMA_Complex64_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_zgels(*M, *N, T); }

void PLASMA_ALLOC_WORKSPACE_ZGEQRF(int *M, int *N, PLASMA_Complex64_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_zgeqrf(*M, *N, T); }

void PLASMA_ALLOC_WORKSPACE_ZGESV_INCPIV(int *N, PLASMA_Complex64_t **L, int **IPIV, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_zgesv_incpiv(*N, L, IPIV); }

void PLASMA_ALLOC_WORKSPACE_ZGESVD(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_zgesvd(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_ZGETRF_INCPIV(int *M, int *N, PLASMA_Complex64_t **L, int **IPIV, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_zgetrf_incpiv(*M, *N, L, IPIV); }

void PLASMA_ALLOC_WORKSPACE_ZHEEV(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_zheev(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_ZHEGV(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_zhegv(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_ZHETRD(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_zhetrd(*M, *N, (PLASMA_desc **)T); }


/***************************************************************************//**
 *  FORTRAN API - tiled workspace allocation
 **/
void PLASMA_ALLOC_WORKSPACE_ZGELQF_TILE(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_zgelqf_Tile(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_ZGELS_TILE(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_zgels_Tile(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_ZGEQRF_TILE(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_zgeqrf_Tile(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_ZGESV_INCPIV_TILE(int *N, intptr_t **L, int **IPIV, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_zgesv_incpiv_Tile(*N, (PLASMA_desc **)L, IPIV); }

void PLASMA_ALLOC_WORKSPACE_ZGETRF_INCPIV_TILE(int *N, intptr_t **L, int **IPIV, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_zgetrf_incpiv_Tile(*N, (PLASMA_desc **)L, IPIV); }

/***************************************************************************//**
 *  FORTRAN API - conversion from LAPACK F77 matrix layout to tile layout
 **/
void PLASMA_ZLAPACK_TO_TILE(PLASMA_Complex64_t **Af77, int *LDA, intptr_t *A, int *INFO)
{   *INFO = PLASMA_zLapack_to_Tile( *Af77, *LDA, (PLASMA_desc *)(*A) ); }

void PLASMA_ZTILE_TO_LAPACK(intptr_t *A, PLASMA_Complex64_t **Af77, int *LDA, int *INFO)
{   *INFO = PLASMA_zTile_to_Lapack( (PLASMA_desc *)(*A), *Af77, *LDA ); }

#ifdef __cplusplus
}
#endif
