/**
 *
 * @file plasma_cf77.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Bilel Hadri
 * @date 2010-11-15
 * @generated c Thu Sep 15 12:09:27 2011
 *
 **/
#include <stdlib.h>
#include "common.h"
#undef REAL
#define COMPLEX

/*
 * Lapack interface
 */
#define PLASMA_CGEBRD        PLASMA_FNAME(cgebrd,  CGEBRD )
#define PLASMA_CGEEV         PLASMA_FNAME(cgeev,   CGEEV  )
#define PLASMA_CGEHRD        PLASMA_FNAME(cgehrd,  CGEHRD )
#define PLASMA_CGELQF        PLASMA_FNAME(cgelqf,  CGELQF )
#define PLASMA_CGELQS        PLASMA_FNAME(cgelqs,  CGELQS )
#define PLASMA_CGELS         PLASMA_FNAME(cgels,   CGELS  )
#define PLASMA_CGEQRF        PLASMA_FNAME(cgeqrf,  CGEQRF )
#define PLASMA_CGEQRS        PLASMA_FNAME(cgeqrs,  CGEQRS )
#define PLASMA_CGESV         PLASMA_FNAME(cgesv,   CGESV  )
#define PLASMA_CGESVD        PLASMA_FNAME(cgesvd,  CGESVD )
#define PLASMA_CGETRF        PLASMA_FNAME(cgetrf,  CGETRF )
#define PLASMA_CGETRS        PLASMA_FNAME(cgetrs,  CGETRS )
#define PLASMA_CGESV_INCPIV  PLASMA_FNAME(cgesv_incpiv,  CGESV_INCPIV  )
#define PLASMA_CGETRF_INCPIV PLASMA_FNAME(cgetrf_incpiv, CGETRF_INCPIV )
#define PLASMA_CGETRS_INCPIV PLASMA_FNAME(cgetrs_incpiv, CGETRS_INCPIV )
#define PLASMA_CHEEV         PLASMA_FNAME(cheev,   CHEEV  )
#define PLASMA_CHEGV         PLASMA_FNAME(chegv,   CHEGV  )
#define PLASMA_CHEGST        PLASMA_FNAME(chegst,  CHEGST )
#define PLASMA_CHETRD        PLASMA_FNAME(chetrd,  CHETRD )
#define PLASMA_CPOSV         PLASMA_FNAME(cposv,   CPOSV  )
#define PLASMA_CPOTRF        PLASMA_FNAME(cpotrf,  CPOTRF )
#define PLASMA_CPOTRI        PLASMA_FNAME(cpotri,  CPOTRI )
#define PLASMA_CPOTRS        PLASMA_FNAME(cpotrs,  CPOTRS )
#define PLASMA_CTRSMPL       PLASMA_FNAME(ctrsmpl, CTRSMPL)
#define PLASMA_CUNGBR        PLASMA_FNAME(cungbr,  CUNGBR )
#define PLASMA_CUNGHR        PLASMA_FNAME(cunghr,  CUNGHR )
#define PLASMA_CUNGLQ        PLASMA_FNAME(cunglq,  CUNGLQ )
#define PLASMA_CUNGQR        PLASMA_FNAME(cungqr,  CUNGQR )
#define PLASMA_CUNGTR        PLASMA_FNAME(cungtr,  CUNGTR )
#define PLASMA_CUNMLQ        PLASMA_FNAME(cunmlq,  CUNMLQ )
#define PLASMA_CUNMQR        PLASMA_FNAME(cunmqr,  CUNMQR )
#define PLASMA_CTRSM         PLASMA_FNAME(ctrsm,   CTRSM  )
#define PLASMA_CGEMM         PLASMA_FNAME(cgemm,   CGEMM  )
#define PLASMA_CSYMM         PLASMA_FNAME(csymm,   CSYMM  )
#define PLASMA_CSYRK         PLASMA_FNAME(csyrk,   CSYRK  )
#ifdef COMPLEX
#define PLASMA_CHEMM         PLASMA_FNAME(chemm,   CHEMM  )
#define PLASMA_CHERK         PLASMA_FNAME(cherk,   CHERK  )
#endif

/*
 * Tile interface
 */
#define PLASMA_CGEBRD_TILE        PLASMA_TILE_FNAME(cgebrd,  CGEBRD )
#define PLASMA_CGEEV_TILE         PLASMA_TILE_FNAME(cgeev,   CGEEV  )
#define PLASMA_CGEHRD_TILE        PLASMA_TILE_FNAME(cgehrd,  CGEHRD )
#define PLASMA_CGELQF_TILE        PLASMA_TILE_FNAME(cgelqf,  CGELQF )
#define PLASMA_CGELQS_TILE        PLASMA_TILE_FNAME(cgelqs,  CGELQS )
#define PLASMA_CGELS_TILE         PLASMA_TILE_FNAME(cgels,   CGELS  )
#define PLASMA_CGEQRF_TILE        PLASMA_TILE_FNAME(cgeqrf,  CGEQRF )
#define PLASMA_CGEQRS_TILE        PLASMA_TILE_FNAME(cgeqrs,  CGEQRS )
#define PLASMA_CGESV_TILE         PLASMA_TILE_FNAME(cgesv,   CGESV  )
#define PLASMA_CGESVD_TILE        PLASMA_TILE_FNAME(cgesvd,  CGESVD )
#define PLASMA_CGETRF_TILE        PLASMA_TILE_FNAME(cgetrf,  CGETRF )
#define PLASMA_CGETRS_TILE        PLASMA_TILE_FNAME(cgetrs,  CGETRS ) 
#define PLASMA_CGESV_INCPIV_TILE  PLASMA_TILE_FNAME(cgesv_incpiv,  CGESV_INCPIV  )
#define PLASMA_CGETRF_INCPIV_TILE PLASMA_TILE_FNAME(cgetrf_incpiv, CGETRF_INCPIV )
#define PLASMA_CGETRS_INCPIV_TILE PLASMA_TILE_FNAME(cgetrs_incpiv, CGETRS_INCPIV ) 
#define PLASMA_CHEEV_TILE         PLASMA_TILE_FNAME(cheev,   CHEEV  )
#define PLASMA_CHEGV_TILE         PLASMA_TILE_FNAME(chegv,   CHEGV  )
#define PLASMA_CHEGST_TILE        PLASMA_TILE_FNAME(chegst,  CHEGST )
#define PLASMA_CHETRD_TILE        PLASMA_TILE_FNAME(chetrd,  CHETRD )
#define PLASMA_CPOSV_TILE         PLASMA_TILE_FNAME(cposv,   CPOSV  )
#define PLASMA_CPOTRF_TILE        PLASMA_TILE_FNAME(cpotrf,  CPOTRF )
#define PLASMA_CPOTRI_TILE        PLASMA_TILE_FNAME(cpotri,  CPOTRI )
#define PLASMA_CPOTRS_TILE        PLASMA_TILE_FNAME(cpotrs,  CPOTRS )
#define PLASMA_CTRSM_TILE         PLASMA_TILE_FNAME(ctrsm,   CTRSM  )
#define PLASMA_CTRSMPL_TILE       PLASMA_TILE_FNAME(ctrsmpl, CTRSMPL)
#define PLASMA_CUNGBR_TILE        PLASMA_TILE_FNAME(cungbr,  CUNGBR )
#define PLASMA_CUNGHR_TILE        PLASMA_TILE_FNAME(cunghr,  CUNGHR )
#define PLASMA_CUNGLQ_TILE        PLASMA_TILE_FNAME(cunglq,  CUNGLQ )
#define PLASMA_CUNGQR_TILE        PLASMA_TILE_FNAME(cungqr,  CUNGQR )
#define PLASMA_CUNGTR_TILE        PLASMA_TILE_FNAME(cungtr,  CUNGTR )
#define PLASMA_CUNMLQ_TILE        PLASMA_TILE_FNAME(cunmlq,  CUNMLQ )
#define PLASMA_CUNMQR_TILE        PLASMA_TILE_FNAME(cunmqr,  CUNMQR )
#define PLASMA_CGEMM_TILE         PLASMA_TILE_FNAME(cgemm,   CGEMM  )
#define PLASMA_CSYMM_TILE         PLASMA_TILE_FNAME(csymm,   CSYMM  )
#define PLASMA_CSYRK_TILE         PLASMA_TILE_FNAME(csyrk,   CSYRK  )
#ifdef COMPLEX                    
#define PLASMA_CHEMM_TILE         PLASMA_TILE_FNAME(chemm,   CHEMM  )
#define PLASMA_CHERK_TILE         PLASMA_TILE_FNAME(cherk,   CHERK  )
#endif

/*
 * Workspaces
 */
#define PLASMA_ALLOC_WORKSPACE_CGEBRD        PLASMA_WS_FNAME(cgehrd, CGEBRD)
#define PLASMA_ALLOC_WORKSPACE_CGEEV         PLASMA_WS_FNAME(cgeev,  CGEEV )
#define PLASMA_ALLOC_WORKSPACE_CGEHRD        PLASMA_WS_FNAME(cgehrd, CGEHRD)
#define PLASMA_ALLOC_WORKSPACE_CGELQF        PLASMA_WS_FNAME(cgelqf, CGELQF) 
#define PLASMA_ALLOC_WORKSPACE_CGELS         PLASMA_WS_FNAME(cgels,  CGELS )  
#define PLASMA_ALLOC_WORKSPACE_CGEQRF        PLASMA_WS_FNAME(cgeqrf, CGEQRF) 
#define PLASMA_ALLOC_WORKSPACE_CGESV_INCPIV  PLASMA_WS_FNAME(cgesv_incpiv,  CGESV_INCPIV ) 
#define PLASMA_ALLOC_WORKSPACE_CGETRF_INCPIV PLASMA_WS_FNAME(cgetrf_incpiv, CGETRF_INCPIV) 
#define PLASMA_ALLOC_WORKSPACE_CGESVD        PLASMA_WS_FNAME(cgesvd, CGESVD)
#define PLASMA_ALLOC_WORKSPACE_CHEEV         PLASMA_WS_FNAME(cheev,  CHEEV )
#define PLASMA_ALLOC_WORKSPACE_CHEGV         PLASMA_WS_FNAME(chegv,  CHEGV )
#define PLASMA_ALLOC_WORKSPACE_CHETRD        PLASMA_WS_FNAME(chetrd, CHETRD)

/*
 * Tile Workspaces
 */
#define PLASMA_ALLOC_WORKSPACE_CGELQF_TILE        PLASMA_WST_FNAME(cgelqf,        CGELQF       )
#define PLASMA_ALLOC_WORKSPACE_CGELS_TILE         PLASMA_WST_FNAME(cgels,         CGELS        )
#define PLASMA_ALLOC_WORKSPACE_CGEQRF_TILE        PLASMA_WST_FNAME(cgeqrf,        CGEQRF       )
#define PLASMA_ALLOC_WORKSPACE_CGESV_INCPIV_TILE  PLASMA_WST_FNAME(cgesv_incpiv,  CGESV_INCPIV )
#define PLASMA_ALLOC_WORKSPACE_CGETRF_INCPIV_TILE PLASMA_WST_FNAME(cgetrf_incpiv, CGETRF_INCPIV)

#define PLASMA_CLAPACK_TO_TILE   PLASMA_FNAME(clapack_to_tile, CLAPACK_TO_TILE)
#define PLASMA_CTILE_TO_LAPACK   PLASMA_FNAME(ctile_to_lapack, CTILE_TO_LAPACK)

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  FORTRAN API - math functions (simple interface)
 **/
void PLASMA_CGEBRD(PLASMA_enum *jobu, PLASMA_enum *jobvt, int *M, int *N, PLASMA_Complex32_t *A, int *LDA, float *D, float *E, PLASMA_Complex32_t *U, int *LDU, PLASMA_Complex32_t *VT, int *LDVT, intptr_t *descT, int *INFO)
{   *INFO = PLASMA_cgebrd(*jobu, *jobvt, *M, *N, A, *LDA, D, E, U, *LDU, VT, *LDVT, (PLASMA_desc *)(*descT)); }

void PLASMA_CGELQF(int *M, int *N, PLASMA_Complex32_t *A, int *LDA, PLASMA_Complex32_t **T, int *INFO)
{   *INFO = PLASMA_cgelqf(*M, *N, A, *LDA, *T); }

void PLASMA_CGELQS(int *M, int *N, int *NRHS, PLASMA_Complex32_t *A, int *LDA, PLASMA_Complex32_t **T, PLASMA_Complex32_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_cgelqs(*M, *N, *NRHS, A, *LDA, *T, B, *LDB); }

void PLASMA_CGELS(PLASMA_enum *trans, int *M, int *N, int *NRHS, PLASMA_Complex32_t *A, int *LDA, PLASMA_Complex32_t **T, PLASMA_Complex32_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_cgels(*trans, *M, *N, *NRHS, A, *LDA, *T, B, *LDB); }

void PLASMA_CGEQRF(int *M, int *N, PLASMA_Complex32_t *A, int *LDA, PLASMA_Complex32_t **T, int *INFO)
{   *INFO = PLASMA_cgeqrf(*M, *N, A, *LDA, *T); }

void PLASMA_CGEQRS(int *M, int *N, int *NRHS, PLASMA_Complex32_t *A, int *LDA, PLASMA_Complex32_t **T, PLASMA_Complex32_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_cgeqrs(*M, *N, *NRHS, A, *LDA, *T, B, *LDB); }

void PLASMA_CGESV(int *N, int *NRHS, PLASMA_Complex32_t *A, int *LDA, int *IPIV, PLASMA_Complex32_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_cgesv(*N, *NRHS, A, *LDA, IPIV, B, *LDB); }

void PLASMA_CGESVD(PLASMA_enum *jobu, PLASMA_enum *jobvt, int *M, int *N, PLASMA_Complex32_t *A, int *LDA, float *S, PLASMA_Complex32_t *U, int *LDU, PLASMA_Complex32_t *VT, int *LDVT, intptr_t *T, int *INFO)
{   *INFO = PLASMA_cgesvd(*jobu, *jobvt, *M, *N, A, *LDA, S, U, *LDU, VT, *LDVT, (PLASMA_desc *)(*T)); }

void PLASMA_CGETRF(int *M, int *N, PLASMA_Complex32_t *A, int *LDA, int *IPIV, int *INFO)
{   *INFO = PLASMA_cgetrf(*M, *N, A, *LDA, IPIV); }

void PLASMA_CGETRS(PLASMA_enum *trans, int *N, int *NRHS, PLASMA_Complex32_t *A, int *LDA, int *IPIV, PLASMA_Complex32_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_cgetrs(*trans, *N, *NRHS, A, *LDA, IPIV, B, *LDB); }

void PLASMA_CGESV_INCPIV(int *N, int *NRHS, PLASMA_Complex32_t *A, int *LDA, PLASMA_Complex32_t **LH, int **IPIVH, PLASMA_Complex32_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_cgesv_incpiv(*N, *NRHS, A, *LDA, *LH, *IPIVH, B, *LDB); }

void PLASMA_CGETRF_INCPIV(int *M, int *N, PLASMA_Complex32_t *A, int *LDA, PLASMA_Complex32_t **LH, int **IPIVH, int *INFO)
{   *INFO = PLASMA_cgetrf_incpiv(*M, *N, A, *LDA, *LH, *IPIVH); }

void PLASMA_CGETRS_INCPIV(PLASMA_enum *uplo, int *N, int *NRHS, PLASMA_Complex32_t *A, int *LDA, PLASMA_Complex32_t **LH, int **IPIVH, PLASMA_Complex32_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_cgetrs_incpiv(*uplo, *N, *NRHS, A, *LDA, *LH, *IPIVH, B, *LDB); }

void PLASMA_CHEEV(PLASMA_enum *jobz, PLASMA_enum *uplo, int *N, PLASMA_Complex32_t *A, int *LDA, float *W, intptr_t *T, PLASMA_Complex32_t *Q, int *LDQ, int *INFO)
{   *INFO = PLASMA_cheev(*jobz, *uplo, *N, A, *LDA, W, (PLASMA_desc *)(*T), Q, *LDQ); }

void PLASMA_CHEGV(PLASMA_enum *itype, PLASMA_enum *jobz, PLASMA_enum *uplo, int *N, PLASMA_Complex32_t *A, int *LDA, PLASMA_Complex32_t *B, int *LDB, float *W, intptr_t *T, PLASMA_Complex32_t *Q, int *LDQ, int *INFO)
{   *INFO = PLASMA_chegv(*itype, *jobz, *uplo, *N, A, *LDA, B, *LDB, W, (PLASMA_desc*)(*T), Q, *LDQ); }

void PLASMA_CHEGST(PLASMA_enum *itype, PLASMA_enum *uplo, int *N, PLASMA_Complex32_t *A, int *LDA, PLASMA_Complex32_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_chegst(*itype, *uplo, *N, A, *LDA, B, *LDB); }

void PLASMA_CHETRD(PLASMA_enum *jobz, PLASMA_enum *uplo, int *N, PLASMA_Complex32_t *A, int *LDA, float *D, float *E, intptr_t *T, PLASMA_Complex32_t *Q, int *LDQ, int *INFO)
{   *INFO = PLASMA_chetrd(*jobz, *uplo, *N, A, *LDA, D, E, (PLASMA_desc *)(*T), Q, *LDQ); }

void PLASMA_CPOSV(PLASMA_enum *uplo, int *N, int *NRHS, PLASMA_Complex32_t *A, int *LDA, PLASMA_Complex32_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_cposv(*uplo, *N, *NRHS, A, *LDA, B, *LDB); }

void PLASMA_CPOTRF(PLASMA_enum *uplo, int *N, PLASMA_Complex32_t *A, int *LDA, int *INFO)
{   *INFO = PLASMA_cpotrf(*uplo, *N, A, *LDA); }

void PLASMA_CPOTRI(PLASMA_enum *uplo, int *N, PLASMA_Complex32_t *A, int *LDA, int *INFO)
{   *INFO = PLASMA_cpotri(*uplo, *N, A, *LDA); }

void PLASMA_CPOTRS(PLASMA_enum *uplo, int *N, int *NRHS, PLASMA_Complex32_t *A, int *LDA, PLASMA_Complex32_t *B, int* LDB, int * INFO)
{   *INFO = PLASMA_cpotrs(*uplo, *N, *NRHS, A, *LDA, B, *LDB); }

void PLASMA_CTRSMPL(int *N, int *NRHS, PLASMA_Complex32_t *A, int *LDA, PLASMA_Complex32_t **LH, int **IPIVH, PLASMA_Complex32_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_ctrsmpl(*N, *NRHS, A, *LDA, *LH, *IPIVH, B, *LDB); }

void PLASMA_CUNGLQ(int *M, int *N, int *K, PLASMA_Complex32_t *A, int *LDA, PLASMA_Complex32_t **T, PLASMA_Complex32_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_cunglq(*M, *N, *K, A, *LDA, *T, B, *LDB); }

void PLASMA_CUNGQR(int *M, int *N, int *K, PLASMA_Complex32_t *A, int *LDA, PLASMA_Complex32_t **T, PLASMA_Complex32_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_cungqr(*M, *N, *K, A, *LDA, *T, B, *LDB); }

void PLASMA_CUNMLQ(PLASMA_enum *side, PLASMA_enum *trans, int *M, int *N, int *K, PLASMA_Complex32_t *A, int *LDA, PLASMA_Complex32_t **T, PLASMA_Complex32_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_cunmlq(*side, *trans, *M, *N, *K, A, *LDA, *T, B, *LDB); }

void PLASMA_CUNMQR(PLASMA_enum *side, PLASMA_enum *trans, int *M, int *N, int *K, PLASMA_Complex32_t *A, int *LDA, PLASMA_Complex32_t **T, PLASMA_Complex32_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_cunmqr(*side, *trans, *M, *N, *K, A, *LDA, *T, B, *LDB); }

void PLASMA_CTRSM(PLASMA_enum *side, PLASMA_enum *uplo, PLASMA_enum *transA, PLASMA_enum *diag, int *N, int *NRHS, PLASMA_Complex32_t *alpha, PLASMA_Complex32_t *A, int *LDA, PLASMA_Complex32_t *B, int *LDB, int *INFO)
{   *INFO = PLASMA_ctrsm(*side, *uplo, *transA, *diag, *N, *NRHS, *alpha, A, *LDA, B, *LDB); }

void PLASMA_CGEMM(PLASMA_enum *transA, PLASMA_enum *transB, int *M, int *N, int *K, PLASMA_Complex32_t *alpha, PLASMA_Complex32_t *A, int *LDA, PLASMA_Complex32_t *B, int *LDB, PLASMA_Complex32_t *beta, PLASMA_Complex32_t *C, int *LDC, int *INFO)
{   *INFO = PLASMA_cgemm(*transA, *transB, *M, *N, *K, *alpha, A, *LDA, B, *LDB, *beta, C, *LDC); }

void PLASMA_CSYMM(PLASMA_enum *side, PLASMA_enum *uplo, int *M, int *N, PLASMA_Complex32_t *alpha, PLASMA_Complex32_t *A, int *LDA, PLASMA_Complex32_t *B, int *LDB, PLASMA_Complex32_t *beta, PLASMA_Complex32_t *C, int *LDC, int *INFO)
{   *INFO = PLASMA_csymm(*side, *uplo, *M, *N, *alpha, A, *LDA, B, *LDB, *beta, C, *LDC); }

void PLASMA_CSYRK(PLASMA_enum *uplo, PLASMA_enum *trans, int *N, int *K, PLASMA_Complex32_t *alpha, PLASMA_Complex32_t *A, int *LDA, PLASMA_Complex32_t *beta, PLASMA_Complex32_t *C, int *LDC, int *INFO)
{   *INFO = PLASMA_csyrk(*uplo, *trans, *N, *K, *alpha, A, *LDA, *beta, C, *LDC); }

#ifdef COMPLEX
void PLASMA_CHEMM(PLASMA_enum *side, PLASMA_enum *uplo, int *M, int *N, PLASMA_Complex32_t *alpha, PLASMA_Complex32_t *A, int *LDA, PLASMA_Complex32_t *B, int *LDB, PLASMA_Complex32_t *beta, PLASMA_Complex32_t *C, int *LDC, int *INFO)
{   *INFO = PLASMA_chemm(*side, *uplo, *M, *N, *alpha, A, *LDA, B, *LDB, *beta, C, *LDC); }

void PLASMA_CHERK(PLASMA_enum *uplo, PLASMA_enum *trans, int *N, int *K, PLASMA_Complex32_t *alpha, PLASMA_Complex32_t *A, int *LDA, float *beta, PLASMA_Complex32_t *C, int *LDC, int *INFO)
{   *INFO = PLASMA_cherk(*uplo, *trans, *N, *K, *alpha, A, *LDA, *beta, C, *LDC); }
#endif

/***************************************************************************//**
 *  FORTRAN API - math functions (native interface)
 **/
void PLASMA_CGEBRD_TILE(PLASMA_enum *jobu, PLASMA_enum *jobvt, intptr_t *A, float *D, float *E, intptr_t *U, intptr_t *VT, intptr_t *T, int *INFO)
{   *INFO = PLASMA_cgebrd_Tile(*jobu, *jobvt, (PLASMA_desc *)(*A), D, E, (PLASMA_desc *)(*U), (PLASMA_desc *)(*VT), (PLASMA_desc *)(*T)); }

void PLASMA_CGELQF_TILE(intptr_t *A, intptr_t *T, int *INFO)
{   *INFO = PLASMA_cgelqf_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*T)); }

void PLASMA_CGELQS_TILE(intptr_t *A, intptr_t *B, intptr_t *T, int *INFO)
{   *INFO = PLASMA_cgelqs_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*B), (PLASMA_desc *)(*T)); }

void PLASMA_CGELS_TILE(PLASMA_enum *trans, intptr_t *A, intptr_t *B, intptr_t *T, int *INFO)
{   *INFO = PLASMA_cgels_Tile(*trans, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B), (PLASMA_desc *)(*T)); }

void PLASMA_CGEQRF_TILE(intptr_t *A, intptr_t *T, int *INFO)
{   *INFO = PLASMA_cgeqrf_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*T)); }

void PLASMA_CGEQRS_TILE(intptr_t *A, intptr_t *B, intptr_t *T, int *INFO)
{   *INFO = PLASMA_cgeqrs_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*B), (PLASMA_desc *)(*T)); }

void PLASMA_CGESV_TILE(intptr_t *A, int *IPIV, intptr_t *B, int *INFO)
{   *INFO = PLASMA_cgesv_Tile((PLASMA_desc *)(*A), IPIV, (PLASMA_desc *)(*B)); }

void PLASMA_CGESVD_TILE(PLASMA_enum *jobu, PLASMA_enum *jobvt, intptr_t *A, float *S, intptr_t *U, intptr_t *VT, intptr_t *T, int *INFO)
{   *INFO = PLASMA_cgesvd_Tile(*jobu, *jobvt, (PLASMA_desc *)(*A), S, (PLASMA_desc *)(*U), (PLASMA_desc *)(*VT), (PLASMA_desc *)(*T)); }

void PLASMA_CGETRF_TILE(intptr_t *A, int *IPIV, int *INFO)
{   *INFO = PLASMA_cgetrf_Tile((PLASMA_desc *)(*A), IPIV); }

void PLASMA_CGETRS_TILE(PLASMA_enum *trans, intptr_t *A, int *IPIV, intptr_t *B, int *INFO)
{   *INFO = PLASMA_cgetrs_Tile(*trans, (PLASMA_desc *)(*A), IPIV, (PLASMA_desc *)(*B)); }

void PLASMA_CGESV_INCPIV_TILE(intptr_t *A, intptr_t *L, int **IPIVH, intptr_t *B, int *INFO)
{   *INFO = PLASMA_cgesv_incpiv_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*L), *IPIVH, (PLASMA_desc *)(*B)); }

void PLASMA_CGETRF_INCPIV_TILE(intptr_t *A, intptr_t *L, int **IPIVH, int *INFO)
{   *INFO = PLASMA_cgetrf_incpiv_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*L), *IPIVH); }

void PLASMA_CGETRS_INCPIV_TILE(intptr_t *A, intptr_t *L, int **IPIVH, intptr_t *B, int *INFO)
{   *INFO = PLASMA_cgetrs_incpiv_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*L), *IPIVH, (PLASMA_desc *)(*B)); }

void PLASMA_CHEEV_TILE(PLASMA_enum *jobz, PLASMA_enum *uplo, intptr_t *A, float *W, intptr_t *T, intptr_t *Q, int *INFO)
{   *INFO = PLASMA_cheev_Tile(*jobz, *uplo, (PLASMA_desc *)(*A), W, (PLASMA_desc *)(*T), (PLASMA_desc *)(*Q)); }

void PLASMA_CHEGV_TILE(PLASMA_enum *itype, PLASMA_enum *jobz, PLASMA_enum *uplo, intptr_t *A, intptr_t *B, float *W, intptr_t *T, intptr_t *Q, int *INFO)
{   *INFO = PLASMA_chegv_Tile(*itype, *jobz, *uplo, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B), W, (PLASMA_desc *)(*T), (PLASMA_desc *)(*Q)); }

void PLASMA_CHEGST_TILE(PLASMA_enum *itype, PLASMA_enum *uplo, intptr_t *A, intptr_t *B, int *INFO)
{   *INFO = PLASMA_chegst_Tile(*itype, *uplo, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B)); }

void PLASMA_CHETRD_TILE(PLASMA_enum *jobz, PLASMA_enum *uplo, intptr_t *A, float *D, float *E, intptr_t *T, intptr_t *Q, int *INFO)
{   *INFO = PLASMA_chetrd_Tile(*jobz, *uplo, (PLASMA_desc *)(*A), D, E, (PLASMA_desc *)(*T), (PLASMA_desc *)(*Q)); }

void PLASMA_CPOSV_TILE(PLASMA_enum *uplo, intptr_t *A, intptr_t *B, int *INFO)
{   *INFO = PLASMA_cposv_Tile(*uplo, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B)); }

void PLASMA_CPOTRF_TILE(PLASMA_enum *uplo, intptr_t *A, int *INFO)
{   *INFO = PLASMA_cpotrf_Tile(*uplo, (PLASMA_desc *)(*A)); }

void PLASMA_CPOTRS_TILE(PLASMA_enum *uplo, intptr_t *A, intptr_t *B, int *INFO)
{   *INFO = PLASMA_cpotrs_Tile(*uplo, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B)); }

void PLASMA_CTRSMPL_TILE(intptr_t *A, intptr_t *L, int **IPIVH, intptr_t *B, int *INFO)
{   *INFO = PLASMA_ctrsmpl_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*L), *IPIVH, (PLASMA_desc *)(*B)); }

void PLASMA_CUNGLQ_TILE(intptr_t *A, intptr_t *T, intptr_t *B, int *INFO)
{   *INFO = PLASMA_cunglq_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*T), (PLASMA_desc *)(*B)); }

void PLASMA_CUNGQR_TILE(intptr_t *A, intptr_t *T, intptr_t *B, int *INFO)
{   *INFO = PLASMA_cungqr_Tile((PLASMA_desc *)(*A), (PLASMA_desc *)(*T), (PLASMA_desc *)(*B)); }

void PLASMA_CUNMLQ_TILE(PLASMA_enum *side, PLASMA_enum *trans, intptr_t *A, intptr_t *T, intptr_t *B, int *INFO)
{   *INFO = PLASMA_cunmlq_Tile(*side, *trans, (PLASMA_desc *)(*A), (PLASMA_desc *)(*T), (PLASMA_desc *)(*B)); }

void PLASMA_CUNMQR_TILE(PLASMA_enum *side, PLASMA_enum *trans, intptr_t *A, intptr_t *T, intptr_t *B, int *INFO)
{   *INFO = PLASMA_cunmqr_Tile(*side, *trans, (PLASMA_desc *)(*A), (PLASMA_desc *)(*T), (PLASMA_desc *)(*B)); }

void PLASMA_CTRSM_TILE(PLASMA_enum *side, PLASMA_enum *uplo, PLASMA_enum *transA, PLASMA_enum *diag, PLASMA_Complex32_t *alpha, intptr_t *A, intptr_t *B, int *INFO)
{   *INFO = PLASMA_ctrsm_Tile(*side, *uplo, *transA, *diag, *alpha, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B)); }

void PLASMA_CGEMM_TILE(PLASMA_enum *transA, PLASMA_enum *transB, int *alpha, intptr_t *A, intptr_t *B, int *beta, intptr_t *C, int *INFO)
{   *INFO = PLASMA_cgemm_Tile(*transA, *transB, *alpha, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B), *beta, (PLASMA_desc *)(*C)); }

  void PLASMA_CSYMM_TILE(PLASMA_enum *side, PLASMA_enum *uplo, PLASMA_Complex32_t *alpha, intptr_t *A, intptr_t *B, PLASMA_Complex32_t *beta, intptr_t *C, int *INFO)
  {   *INFO = PLASMA_csymm_Tile(*side, *uplo, *alpha, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B), *beta, (PLASMA_desc *)(*C)); }

void PLASMA_CSYRK_TILE(PLASMA_enum *uplo, PLASMA_enum *trans, PLASMA_Complex32_t *alpha, intptr_t *A, PLASMA_Complex32_t *beta, intptr_t *C, int *INFO)
{   *INFO = PLASMA_csyrk_Tile(*uplo, *trans, *alpha, (PLASMA_desc *)(*A), *beta, (PLASMA_desc *)(*C)); }

#ifdef COMPLEX
void PLASMA_CHEMM_TILE(PLASMA_enum *side, PLASMA_enum *uplo, PLASMA_Complex32_t *alpha, intptr_t *A, intptr_t *B, PLASMA_Complex32_t *beta, intptr_t *C, int *INFO)
{   *INFO = PLASMA_chemm_Tile(*side, *uplo, *alpha, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B), *beta, (PLASMA_desc *)(*C)); }

void PLASMA_CHERK_TILE(PLASMA_enum *uplo, PLASMA_enum *trans, PLASMA_Complex32_t *alpha, intptr_t *A, float *beta, intptr_t *C, int *INFO)
{   *INFO = PLASMA_cherk_Tile(*uplo, *trans, *alpha, (PLASMA_desc *)(*A), *beta, (PLASMA_desc *)(*C)); }
#endif

/***************************************************************************//**
 *  FORTRAN API - workspace allocation
 **/
void PLASMA_ALLOC_WORKSPACE_CGEBRD(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_cgebrd(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_CGELQF(int *M, int *N, PLASMA_Complex32_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_cgelqf(*M, *N, T); }

void PLASMA_ALLOC_WORKSPACE_CGELS(int *M, int *N, PLASMA_Complex32_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_cgels(*M, *N, T); }

void PLASMA_ALLOC_WORKSPACE_CGEQRF(int *M, int *N, PLASMA_Complex32_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_cgeqrf(*M, *N, T); }

void PLASMA_ALLOC_WORKSPACE_CGESV_INCPIV(int *N, PLASMA_Complex32_t **L, int **IPIV, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_cgesv_incpiv(*N, L, IPIV); }

void PLASMA_ALLOC_WORKSPACE_CGESVD(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_cgesvd(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_CGETRF_INCPIV(int *M, int *N, PLASMA_Complex32_t **L, int **IPIV, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_cgetrf_incpiv(*M, *N, L, IPIV); }

void PLASMA_ALLOC_WORKSPACE_CHEEV(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_cheev(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_CHEGV(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_chegv(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_CHETRD(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_chetrd(*M, *N, (PLASMA_desc **)T); }


/***************************************************************************//**
 *  FORTRAN API - tiled workspace allocation
 **/
void PLASMA_ALLOC_WORKSPACE_CGELQF_TILE(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_cgelqf_Tile(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_CGELS_TILE(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_cgels_Tile(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_CGEQRF_TILE(int *M, int *N, intptr_t **T, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_cgeqrf_Tile(*M, *N, (PLASMA_desc **)T); }

void PLASMA_ALLOC_WORKSPACE_CGESV_INCPIV_TILE(int *N, intptr_t **L, int **IPIV, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_cgesv_incpiv_Tile(*N, (PLASMA_desc **)L, IPIV); }

void PLASMA_ALLOC_WORKSPACE_CGETRF_INCPIV_TILE(int *N, intptr_t **L, int **IPIV, int *INFO)
{   *INFO = PLASMA_Alloc_Workspace_cgetrf_incpiv_Tile(*N, (PLASMA_desc **)L, IPIV); }

/***************************************************************************//**
 *  FORTRAN API - conversion from LAPACK F77 matrix layout to tile layout
 **/
void PLASMA_CLAPACK_TO_TILE(PLASMA_Complex32_t **Af77, int *LDA, intptr_t *A, int *INFO)
{   *INFO = PLASMA_cLapack_to_Tile( *Af77, *LDA, (PLASMA_desc *)(*A) ); }

void PLASMA_CTILE_TO_LAPACK(intptr_t *A, PLASMA_Complex32_t **Af77, int *LDA, int *INFO)
{   *INFO = PLASMA_cTile_to_Lapack( (PLASMA_desc *)(*A), *Af77, *LDA ); }

#ifdef __cplusplus
}
#endif
