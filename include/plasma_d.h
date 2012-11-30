/**
 *
 * @file plasma_d.h
 *
 *  PLASMA header file for double routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Azzam Haidar
 * @date 2010-11-15
 * @generated d Thu Sep 15 12:08:48 2011
 *
 **/
#ifndef _PLASMA_D_H_
#define _PLASMA_D_H_
#undef COMPLEX
#define REAL

#ifdef __cplusplus
extern "C" {
#endif

/** ****************************************************************************
 *  Declarations of math functions (LAPACK layout) - alphabetical order
 **/
int PLASMA_dgebrd(PLASMA_enum jobu, PLASMA_enum jobvt, int M, int N, double *A, int LDA, double *D, double *E, double *U, int LDU, double *VT, int LDVT, PLASMA_desc *T);
int PLASMA_dgeev(PLASMA_enum jobvl, PLASMA_enum jobvr, int N, double *A, int LDA, double *W, double *VL, int LDVL, double *VR, int LDVR, double *T);
int PLASMA_dgehrd(int N, int ILO, int IHI, double *A, int LDA, double *T);
int PLASMA_dgelqf(int M, int N, double *A, int LDA, double *T);
int PLASMA_dgelqs(int M, int N, int NRHS, double *A, int LDA, double *T, double *B, int LDB);
int PLASMA_dgels(PLASMA_enum trans, int M, int N, int NRHS, double *A, int LDA, double *T, double *B, int LDB);
int PLASMA_dgemm(PLASMA_enum transA, PLASMA_enum transB, int M, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC);
int PLASMA_dgeqrf(int M, int N, double *A, int LDA, double *T);
int PLASMA_dgeqrs(int M, int N, int NRHS, double *A, int LDA, double *T, double *B, int LDB);
int PLASMA_dgesv(int N, int NRHS, double *A, int LDA, int *IPIV, double *B, int LDB);
int PLASMA_dgesv_incpiv(int N, int NRHS, double *A, int LDA, double *L, int *IPIV, double *B, int LDB);
int PLASMA_dgesvd(PLASMA_enum jobu, PLASMA_enum jobvt, int M, int N, double *A, int LDA, double *S, double *U, int LDU, double *VT, int LDVT, PLASMA_desc *T);
int PLASMA_dgetrf(int M, int N, double *A, int LDA, int *IPIV);
int PLASMA_dgetrf_incpiv(int M, int N, double *A, int LDA, double *L, int *IPIV);
int PLASMA_dgetrs(PLASMA_enum trans, int N, int NRHS, double *A, int LDA, int *IPIV, double *B, int LDB);
int PLASMA_dgetrs_incpiv(PLASMA_enum trans, int N, int NRHS, double *A, int LDA, double *L, int *IPIV, double *B, int LDB);
#ifdef COMPLEX
int PLASMA_dsymm(PLASMA_enum side, PLASMA_enum uplo, int M, int N, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC);
int PLASMA_dsyrk(PLASMA_enum uplo, PLASMA_enum trans, int N, int K, double alpha, double *A, int LDA, double beta, double *C, int LDC);
int PLASMA_dsyr2k(PLASMA_enum uplo, PLASMA_enum trans, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC);
#endif
int PLASMA_dsyev(PLASMA_enum jobz, PLASMA_enum uplo, int N, double *A, int LDA, double *W, PLASMA_desc *T, double *Q, int LDQ);
int PLASMA_dsygv(PLASMA_enum itype, PLASMA_enum jobz, PLASMA_enum uplo, int N, double *A, int LDA, double *B, int LDB, double *W, PLASMA_desc *T, double *Q, int LDQ);
int PLASMA_dsygst(PLASMA_enum itype, PLASMA_enum uplo, int N, double *A, int LDA, double *B, int LDB);
int PLASMA_dsytrd(PLASMA_enum jobz, PLASMA_enum uplo, int N, double *A, int LDA, double *D, double *E, PLASMA_desc *T, double *Q, int LDQ);
double PLASMA_dlange(PLASMA_enum norm, int M, int N, double *A, int LDA, double *work);
#ifdef COMPLEX
double PLASMA_dlansy(PLASMA_enum norm, PLASMA_enum uplo, int N, double *A, int LDA, double *work);
#endif
double PLASMA_dlansy(PLASMA_enum norm, PLASMA_enum uplo, int N, double *A, int LDA, double *work);
int PLASMA_dlaswp(int N, double *A, int LDA, int K1, int K2, int *IPIV, int INCX);
int PLASMA_dlauum(PLASMA_enum uplo, int N, double *A, int LDA);
#ifdef COMPLEX
int PLASMA_dplgsy( double bump, int N, double *A, int LDA, unsigned long long int seed );
#endif
int PLASMA_dplgsy( double bump, int N, double *A, int LDA, unsigned long long int seed );
int PLASMA_dplrnt( int M, int N, double *A, int LDA, unsigned long long int seed );
int PLASMA_dposv(PLASMA_enum uplo, int N, int NRHS, double *A, int LDA, double *B, int LDB);
int PLASMA_dpotrf(PLASMA_enum uplo, int N, double *A, int LDA);
int PLASMA_dpotri(PLASMA_enum uplo, int N, double *A, int LDA);
int PLASMA_dpotrs(PLASMA_enum uplo, int N, int NRHS, double *A, int LDA, double *B, int LDB);
int PLASMA_dsymm(PLASMA_enum side, PLASMA_enum uplo, int M, int N, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC);
int PLASMA_dsyrk(PLASMA_enum uplo, PLASMA_enum trans, int N, int K, double alpha, double *A, int LDA, double beta, double *C, int LDC);
int PLASMA_dsyr2k(PLASMA_enum uplo, PLASMA_enum trans, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC);
int PLASMA_dtrmm(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, int N, int NRHS, double alpha, double *A, int LDA, double *B, int LDB);
int PLASMA_dtrsm(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, int N, int NRHS, double alpha, double *A, int LDA, double *B, int LDB);
int PLASMA_dtrsmpl(int N, int NRHS, double *A, int LDA, double *L, int *IPIV, double *B, int LDB);
int PLASMA_dtrtri(PLASMA_enum uplo, PLASMA_enum diag, int N, double *A, int LDA);
int PLASMA_dorgbr(PLASMA_enum side, int M, int N, int K, double *A, int LDA, double *T, double *Q, int LDQ);
int PLASMA_dorghr(int N, int ILO, int IHI, double *A, int LDA, double *T, double *Q, int LDQ);
int PLASMA_dorglq(int M, int N, int K, double *A, int LDA, double *T, double *B, int LDB);
int PLASMA_dorgqr(int M, int N, int K, double *A, int LDA, double *T, double *B, int LDB);
int PLASMA_dorgtr(PLASMA_enum uplo, int N, double *A, int LDA, double *T, double *B, int LDB);
int PLASMA_dormlq(PLASMA_enum side, PLASMA_enum trans, int M, int N, int K, double *A, int LDA, double *T, double *B, int LDB);
int PLASMA_dormqr(PLASMA_enum side, PLASMA_enum trans, int M, int N, int K, double *A, int LDA, double *T, double *B, int LDB);

int PLASMA_dgecfi(int m, int n, double *A, PLASMA_enum fin, int imb, int inb, PLASMA_enum fout, int omb, int onb);
int PLASMA_dgetmi(int m, int n, double *A, PLASMA_enum fin, int mb,  int nb);

/** ****************************************************************************
 *  Declarations of math functions (tile layout) - alphabetical order
 **/
int PLASMA_dgebrd_Tile(PLASMA_enum jobu, PLASMA_enum jobvt, PLASMA_desc *A, double *D, double *E, PLASMA_desc *U, PLASMA_desc *VT, PLASMA_desc *T);
int PLASMA_dgeev_Tile(PLASMA_enum jobvl, PLASMA_enum jobvr, PLASMA_desc *A, double *W, PLASMA_desc *VL, PLASMA_desc *VR, PLASMA_desc *T);
int PLASMA_dgehrd_Tile(PLASMA_desc *A, PLASMA_desc *T);
int PLASMA_dgelqf_Tile(PLASMA_desc *A, PLASMA_desc *T);
int PLASMA_dgelqs_Tile(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B);
int PLASMA_dgels_Tile(PLASMA_enum trans, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B);
int PLASMA_dgemm_Tile(PLASMA_enum transA, PLASMA_enum transB, double alpha, PLASMA_desc *A, PLASMA_desc *B, double beta, PLASMA_desc *C);
int PLASMA_dgeqrf_Tile(PLASMA_desc *A, PLASMA_desc *T);
int PLASMA_dgeqrs_Tile(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B);
int PLASMA_dgesv_Tile(PLASMA_desc *A, int *IPIV, PLASMA_desc *B);
int PLASMA_dgesv_incpiv_Tile(PLASMA_desc *A, PLASMA_desc *L, int *IPIV, PLASMA_desc *B);
int PLASMA_dgesvd_Tile(PLASMA_enum jobu, PLASMA_enum jobvt, PLASMA_desc *A, double *S, PLASMA_desc *U, PLASMA_desc *VT, PLASMA_desc *T);
int PLASMA_dgetrf_Tile(PLASMA_desc *A, int *IPIV);
int PLASMA_dgetrf_incpiv_Tile(PLASMA_desc *A, PLASMA_desc *L, int *IPIV);
int PLASMA_dgetrs_Tile(PLASMA_enum trans, PLASMA_desc *A, int *IPIV, PLASMA_desc *B);
int PLASMA_dgetrs_incpiv_Tile(PLASMA_desc *A, PLASMA_desc *L, int *IPIV, PLASMA_desc *B);
#ifdef COMPLEX
int PLASMA_dsymm_Tile(PLASMA_enum side, PLASMA_enum uplo, double alpha, PLASMA_desc *A, PLASMA_desc *B, double beta, PLASMA_desc *C);
int PLASMA_dsyrk_Tile(PLASMA_enum uplo, PLASMA_enum trans, double alpha, PLASMA_desc *A, double beta, PLASMA_desc *C);
int PLASMA_dsyr2k_Tile(PLASMA_enum uplo, PLASMA_enum trans, double alpha, PLASMA_desc *A, PLASMA_desc *B, double beta, PLASMA_desc *C);
#endif
int PLASMA_dsyev_Tile(PLASMA_enum jobz, PLASMA_enum uplo, PLASMA_desc *A, double *W, PLASMA_desc *T, PLASMA_desc *Q);
int PLASMA_dsygv_Tile(PLASMA_enum itype, PLASMA_enum jobz, PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B, double *W, PLASMA_desc *T, PLASMA_desc *Q);
int PLASMA_dsygst_Tile(PLASMA_enum itype, PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B);
int PLASMA_dsytrd_Tile(PLASMA_enum jobz, PLASMA_enum uplo, PLASMA_desc *A, double *D, double *E, PLASMA_desc *T, PLASMA_desc *Q);
double PLASMA_dlange_Tile(PLASMA_enum norm, PLASMA_desc *A, double *work);
#ifdef COMPLEX
double PLASMA_dlansy_Tile(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc *A, double *work);
#endif
double PLASMA_dlansy_Tile(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc *A, double *work);
int PLASMA_dlaswp_Tile(PLASMA_desc *A, int K1, int K2, int *IPIV, int INCX);
int PLASMA_dlauum_Tile(PLASMA_enum uplo, PLASMA_desc *A);
#ifdef COMPLEX
int PLASMA_dplgsy_Tile(double bump, PLASMA_desc *A, unsigned long long int seed );
#endif
int PLASMA_dplgsy_Tile(double bump, PLASMA_desc *A, unsigned long long int seed );
int PLASMA_dplrnt_Tile(PLASMA_desc *A, unsigned long long int seed );
int PLASMA_dposv_Tile(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B);
int PLASMA_dpotrf_Tile(PLASMA_enum uplo, PLASMA_desc *A);
int PLASMA_dpotri_Tile(PLASMA_enum uplo, PLASMA_desc *A);
int PLASMA_dpotrs_Tile(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B);
int PLASMA_dsymm_Tile(PLASMA_enum side, PLASMA_enum uplo, double alpha, PLASMA_desc *A, PLASMA_desc *B, double beta, PLASMA_desc *C);
int PLASMA_dsyrk_Tile(PLASMA_enum uplo, PLASMA_enum trans, double alpha, PLASMA_desc *A, double beta, PLASMA_desc *C);
int PLASMA_dsyr2k_Tile(PLASMA_enum uplo, PLASMA_enum trans, double alpha, PLASMA_desc *A, PLASMA_desc *B, double beta, PLASMA_desc *C);
int PLASMA_dtrmm_Tile(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, double alpha, PLASMA_desc *A, PLASMA_desc *B);
int PLASMA_dtrsm_Tile(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, double alpha, PLASMA_desc *A, PLASMA_desc *B);
int PLASMA_dtrsmpl_Tile(PLASMA_desc *A, PLASMA_desc *L, int *IPIV, PLASMA_desc *B);
int PLASMA_dtrtri_Tile(PLASMA_enum uplo, PLASMA_enum diag, PLASMA_desc *A);
int PLASMA_dorgbr_Tile(PLASMA_enum size, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *Q);
int PLASMA_dorghr_Tile(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *Q);
int PLASMA_dorglq_Tile(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B);
int PLASMA_dorgqr_Tile(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B);
int PLASMA_dorgtr_Tile(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B);
int PLASMA_dormlq_Tile(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B);
int PLASMA_dormqr_Tile(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B);

/** ****************************************************************************
 *  Declarations of math functions (tile layout, asynchronous execution) - alphabetical order
 **/
int PLASMA_dgebrd_Tile_Async(PLASMA_enum jobu, PLASMA_enum jobvt, PLASMA_desc *A, double *D, double *E, PLASMA_desc *U, PLASMA_desc *VT, PLASMA_desc *T, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dgeev_Tile_Async(PLASMA_enum jobvl, PLASMA_enum jobvr, PLASMA_desc *A, double *W, PLASMA_desc *VL, PLASMA_desc *VR, PLASMA_desc *T, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dgehrd_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dgelqf_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dgelqs_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dgels_Tile_Async(PLASMA_enum trans, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dgemm_Tile_Async(PLASMA_enum transA, PLASMA_enum transB, double alpha, PLASMA_desc *A, PLASMA_desc *B, double beta, PLASMA_desc *C, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dgeqrf_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dgeqrs_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dgesv_Tile_Async(PLASMA_desc *A, int *IPIV, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dgesv_incpiv_Tile_Async(PLASMA_desc *A, PLASMA_desc *L, int *IPIV, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dgesvd_Tile_Async(PLASMA_enum jobu, PLASMA_enum jobvt, PLASMA_desc *A, double *S, PLASMA_desc *U, PLASMA_desc *VT, PLASMA_desc *T, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dgetrf_Tile_Async(PLASMA_desc *A, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dgetrf_incpiv_Tile_Async(PLASMA_desc *A, PLASMA_desc *L, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dgetrs_Tile_Async(PLASMA_enum trans, PLASMA_desc *A, int *IPIV, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dgetrs_incpiv_Tile_Async(PLASMA_desc *A, PLASMA_desc *L, int *IPIV, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
#ifdef COMPLEX
int PLASMA_dsymm_Tile_Async(PLASMA_enum side, PLASMA_enum uplo, double alpha, PLASMA_desc *A, PLASMA_desc *B, double beta, PLASMA_desc *C, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dsyrk_Tile_Async(PLASMA_enum uplo, PLASMA_enum trans, double alpha, PLASMA_desc *A, double beta, PLASMA_desc *C, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dsyr2k_Tile_Async(PLASMA_enum uplo, PLASMA_enum trans, double alpha, PLASMA_desc *A, PLASMA_desc *B, double beta, PLASMA_desc *C, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
int PLASMA_dsyev_Tile_Async(PLASMA_enum jobz, PLASMA_enum uplo, PLASMA_desc *A, double *W, PLASMA_desc *T, PLASMA_desc *Q, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dsygv_Tile_Async(PLASMA_enum itype, PLASMA_enum jobz, PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B, double *W, PLASMA_desc *T, PLASMA_desc *Q, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dsygst_Tile_Async(PLASMA_enum itype, PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dsytrd_Tile_Async(PLASMA_enum jobz, PLASMA_enum uplo, PLASMA_desc *A, double *D, double *E, PLASMA_desc *T, PLASMA_desc *Q, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dlange_Tile_Async(PLASMA_enum norm, PLASMA_desc *A, double *work, double *value, PLASMA_sequence *sequence, PLASMA_request *request);
#ifdef COMPLEX
int PLASMA_dlansy_Tile_Async(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc *A, double *work, double *value, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
int PLASMA_dlansy_Tile_Async(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc *A, double *work, double *value, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dlaswp_Tile_Async(PLASMA_desc *A, int K1, int K2, int *IPIV, int INCX, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dlauum_Tile_Async(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_sequence *sequence, PLASMA_request *request);
#ifdef COMPLEX
int PLASMA_dplgsy_Tile_Async(double bump, PLASMA_desc *A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
#endif
int PLASMA_dplgsy_Tile_Async(double bump, PLASMA_desc *A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
int PLASMA_dplrnt_Tile_Async(PLASMA_desc *A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
int PLASMA_dposv_Tile_Async(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dpotrf_Tile_Async(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dpotri_Tile_Async(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dpotrs_Tile_Async(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dsymm_Tile_Async(PLASMA_enum side, PLASMA_enum uplo, double alpha, PLASMA_desc *A, PLASMA_desc *B, double beta, PLASMA_desc *C, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dsyrk_Tile_Async(PLASMA_enum uplo, PLASMA_enum trans, double alpha, PLASMA_desc *A, double beta, PLASMA_desc *C, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dsyr2k_Tile_Async(PLASMA_enum uplo, PLASMA_enum trans, double alpha, PLASMA_desc *A, PLASMA_desc *B, double beta, PLASMA_desc *C, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dtrmm_Tile_Async(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, double alpha, PLASMA_desc *A, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dtrsm_Tile_Async(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, double alpha, PLASMA_desc *A, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dtrsmpl_Tile_Async(PLASMA_desc *A, PLASMA_desc *L, int *IPIV, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dtrtri_Tile_Async(PLASMA_enum uplo, PLASMA_enum diag, PLASMA_desc *A, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dorgbr_Tile_Async(PLASMA_enum side, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *Q, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dorghr_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *Q, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dorglq_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dorgqr_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dorgtr_Tile_Async(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dormlq_Tile_Async(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dormqr_Tile_Async(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);

int PLASMA_dgecfi_Async(int m, int n, double *A, PLASMA_enum f_in, int imb, int inb, PLASMA_enum f_out, int omb, int onb, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dgetmi_Async(int m, int n, double *A, PLASMA_enum f_in, int mb,  int inb, PLASMA_sequence *sequence, PLASMA_request *request);

/** ****************************************************************************
 *  Declarations of workspace allocation functions (tile layout) - alphabetical order
 **/
int PLASMA_Alloc_Workspace_dgelqf(int M, int N, double **T);
int PLASMA_Alloc_Workspace_dgels( int M, int N, double **T);
int PLASMA_Alloc_Workspace_dgeqrf(int M, int N, double **T);
int PLASMA_Alloc_Workspace_dgesv_incpiv(        int N, double **L, int **IPIV);
int PLASMA_Alloc_Workspace_dgetrf_incpiv(int M, int N, double **L, int **IPIV);

int PLASMA_Alloc_Workspace_dgebrd(int M, int N, PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_dgeev( int N,        PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_dgehrd(int N,        PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_dgesvd(int M, int N, PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_dsyev( int M, int N, PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_dsygv( int M, int N, PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_dsytrd(int M, int N, PLASMA_desc **descT);

/** ****************************************************************************
 *  Declarations of workspace allocation functions (tile layout, asynchronous execution) - alphabetical order
 **/
int PLASMA_Alloc_Workspace_dgelqf_Tile(int M, int N, PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_dgels_Tile( int M, int N, PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_dgeqrf_Tile(int M, int N, PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_dgesv_incpiv_Tile (int N, PLASMA_desc **descL, int **IPIV);
int PLASMA_Alloc_Workspace_dgetrf_incpiv_Tile(int N, PLASMA_desc **descL, int **IPIV);

/** ****************************************************************************
 *  Auxiliary function prototypes
 **/
int PLASMA_dLapack_to_Tile(double *Af77, int LDA, PLASMA_desc *A);
int PLASMA_dTile_to_Lapack(PLASMA_desc *A, double *Af77, int LDA);
int PLASMA_dLapack_to_Tile_Async(double *Af77, int LDA, PLASMA_desc *A, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_dTile_to_Lapack_Async(PLASMA_desc *A, double *Af77, int LDA, PLASMA_sequence *sequence, PLASMA_request *request);

#ifdef __cplusplus
}
#endif

#endif
