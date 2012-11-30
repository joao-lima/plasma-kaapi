/**
 *
 * @file plasma_s.h
 *
 *  PLASMA header file for float routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Azzam Haidar
 * @date 2010-11-15
 * @generated s Thu Sep 15 12:08:48 2011
 *
 **/
#ifndef _PLASMA_S_H_
#define _PLASMA_S_H_
#undef COMPLEX
#define REAL

#ifdef __cplusplus
extern "C" {
#endif

/** ****************************************************************************
 *  Declarations of math functions (LAPACK layout) - alphabetical order
 **/
int PLASMA_sgebrd(PLASMA_enum jobu, PLASMA_enum jobvt, int M, int N, float *A, int LDA, float *D, float *E, float *U, int LDU, float *VT, int LDVT, PLASMA_desc *T);
int PLASMA_sgeev(PLASMA_enum jobvl, PLASMA_enum jobvr, int N, float *A, int LDA, float *W, float *VL, int LDVL, float *VR, int LDVR, float *T);
int PLASMA_sgehrd(int N, int ILO, int IHI, float *A, int LDA, float *T);
int PLASMA_sgelqf(int M, int N, float *A, int LDA, float *T);
int PLASMA_sgelqs(int M, int N, int NRHS, float *A, int LDA, float *T, float *B, int LDB);
int PLASMA_sgels(PLASMA_enum trans, int M, int N, int NRHS, float *A, int LDA, float *T, float *B, int LDB);
int PLASMA_sgemm(PLASMA_enum transA, PLASMA_enum transB, int M, int N, int K, float alpha, float *A, int LDA, float *B, int LDB, float beta, float *C, int LDC);
int PLASMA_sgeqrf(int M, int N, float *A, int LDA, float *T);
int PLASMA_sgeqrs(int M, int N, int NRHS, float *A, int LDA, float *T, float *B, int LDB);
int PLASMA_sgesv(int N, int NRHS, float *A, int LDA, int *IPIV, float *B, int LDB);
int PLASMA_sgesv_incpiv(int N, int NRHS, float *A, int LDA, float *L, int *IPIV, float *B, int LDB);
int PLASMA_sgesvd(PLASMA_enum jobu, PLASMA_enum jobvt, int M, int N, float *A, int LDA, float *S, float *U, int LDU, float *VT, int LDVT, PLASMA_desc *T);
int PLASMA_sgetrf(int M, int N, float *A, int LDA, int *IPIV);
int PLASMA_sgetrf_incpiv(int M, int N, float *A, int LDA, float *L, int *IPIV);
int PLASMA_sgetrs(PLASMA_enum trans, int N, int NRHS, float *A, int LDA, int *IPIV, float *B, int LDB);
int PLASMA_sgetrs_incpiv(PLASMA_enum trans, int N, int NRHS, float *A, int LDA, float *L, int *IPIV, float *B, int LDB);
#ifdef COMPLEX
int PLASMA_ssymm(PLASMA_enum side, PLASMA_enum uplo, int M, int N, float alpha, float *A, int LDA, float *B, int LDB, float beta, float *C, int LDC);
int PLASMA_ssyrk(PLASMA_enum uplo, PLASMA_enum trans, int N, int K, float alpha, float *A, int LDA, float beta, float *C, int LDC);
int PLASMA_ssyr2k(PLASMA_enum uplo, PLASMA_enum trans, int N, int K, float alpha, float *A, int LDA, float *B, int LDB, float beta, float *C, int LDC);
#endif
int PLASMA_ssyev(PLASMA_enum jobz, PLASMA_enum uplo, int N, float *A, int LDA, float *W, PLASMA_desc *T, float *Q, int LDQ);
int PLASMA_ssygv(PLASMA_enum itype, PLASMA_enum jobz, PLASMA_enum uplo, int N, float *A, int LDA, float *B, int LDB, float *W, PLASMA_desc *T, float *Q, int LDQ);
int PLASMA_ssygst(PLASMA_enum itype, PLASMA_enum uplo, int N, float *A, int LDA, float *B, int LDB);
int PLASMA_ssytrd(PLASMA_enum jobz, PLASMA_enum uplo, int N, float *A, int LDA, float *D, float *E, PLASMA_desc *T, float *Q, int LDQ);
float PLASMA_slange(PLASMA_enum norm, int M, int N, float *A, int LDA, float *work);
#ifdef COMPLEX
float PLASMA_slansy(PLASMA_enum norm, PLASMA_enum uplo, int N, float *A, int LDA, float *work);
#endif
float PLASMA_slansy(PLASMA_enum norm, PLASMA_enum uplo, int N, float *A, int LDA, float *work);
int PLASMA_slaswp(int N, float *A, int LDA, int K1, int K2, int *IPIV, int INCX);
int PLASMA_slauum(PLASMA_enum uplo, int N, float *A, int LDA);
#ifdef COMPLEX
int PLASMA_splgsy( float bump, int N, float *A, int LDA, unsigned long long int seed );
#endif
int PLASMA_splgsy( float bump, int N, float *A, int LDA, unsigned long long int seed );
int PLASMA_splrnt( int M, int N, float *A, int LDA, unsigned long long int seed );
int PLASMA_sposv(PLASMA_enum uplo, int N, int NRHS, float *A, int LDA, float *B, int LDB);
int PLASMA_spotrf(PLASMA_enum uplo, int N, float *A, int LDA);
int PLASMA_spotri(PLASMA_enum uplo, int N, float *A, int LDA);
int PLASMA_spotrs(PLASMA_enum uplo, int N, int NRHS, float *A, int LDA, float *B, int LDB);
int PLASMA_ssymm(PLASMA_enum side, PLASMA_enum uplo, int M, int N, float alpha, float *A, int LDA, float *B, int LDB, float beta, float *C, int LDC);
int PLASMA_ssyrk(PLASMA_enum uplo, PLASMA_enum trans, int N, int K, float alpha, float *A, int LDA, float beta, float *C, int LDC);
int PLASMA_ssyr2k(PLASMA_enum uplo, PLASMA_enum trans, int N, int K, float alpha, float *A, int LDA, float *B, int LDB, float beta, float *C, int LDC);
int PLASMA_strmm(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, int N, int NRHS, float alpha, float *A, int LDA, float *B, int LDB);
int PLASMA_strsm(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, int N, int NRHS, float alpha, float *A, int LDA, float *B, int LDB);
int PLASMA_strsmpl(int N, int NRHS, float *A, int LDA, float *L, int *IPIV, float *B, int LDB);
int PLASMA_strtri(PLASMA_enum uplo, PLASMA_enum diag, int N, float *A, int LDA);
int PLASMA_sorgbr(PLASMA_enum side, int M, int N, int K, float *A, int LDA, float *T, float *Q, int LDQ);
int PLASMA_sorghr(int N, int ILO, int IHI, float *A, int LDA, float *T, float *Q, int LDQ);
int PLASMA_sorglq(int M, int N, int K, float *A, int LDA, float *T, float *B, int LDB);
int PLASMA_sorgqr(int M, int N, int K, float *A, int LDA, float *T, float *B, int LDB);
int PLASMA_sorgtr(PLASMA_enum uplo, int N, float *A, int LDA, float *T, float *B, int LDB);
int PLASMA_sormlq(PLASMA_enum side, PLASMA_enum trans, int M, int N, int K, float *A, int LDA, float *T, float *B, int LDB);
int PLASMA_sormqr(PLASMA_enum side, PLASMA_enum trans, int M, int N, int K, float *A, int LDA, float *T, float *B, int LDB);

int PLASMA_sgecfi(int m, int n, float *A, PLASMA_enum fin, int imb, int inb, PLASMA_enum fout, int omb, int onb);
int PLASMA_sgetmi(int m, int n, float *A, PLASMA_enum fin, int mb,  int nb);

/** ****************************************************************************
 *  Declarations of math functions (tile layout) - alphabetical order
 **/
int PLASMA_sgebrd_Tile(PLASMA_enum jobu, PLASMA_enum jobvt, PLASMA_desc *A, float *D, float *E, PLASMA_desc *U, PLASMA_desc *VT, PLASMA_desc *T);
int PLASMA_sgeev_Tile(PLASMA_enum jobvl, PLASMA_enum jobvr, PLASMA_desc *A, float *W, PLASMA_desc *VL, PLASMA_desc *VR, PLASMA_desc *T);
int PLASMA_sgehrd_Tile(PLASMA_desc *A, PLASMA_desc *T);
int PLASMA_sgelqf_Tile(PLASMA_desc *A, PLASMA_desc *T);
int PLASMA_sgelqs_Tile(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B);
int PLASMA_sgels_Tile(PLASMA_enum trans, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B);
int PLASMA_sgemm_Tile(PLASMA_enum transA, PLASMA_enum transB, float alpha, PLASMA_desc *A, PLASMA_desc *B, float beta, PLASMA_desc *C);
int PLASMA_sgeqrf_Tile(PLASMA_desc *A, PLASMA_desc *T);
int PLASMA_sgeqrs_Tile(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B);
int PLASMA_sgesv_Tile(PLASMA_desc *A, int *IPIV, PLASMA_desc *B);
int PLASMA_sgesv_incpiv_Tile(PLASMA_desc *A, PLASMA_desc *L, int *IPIV, PLASMA_desc *B);
int PLASMA_sgesvd_Tile(PLASMA_enum jobu, PLASMA_enum jobvt, PLASMA_desc *A, float *S, PLASMA_desc *U, PLASMA_desc *VT, PLASMA_desc *T);
int PLASMA_sgetrf_Tile(PLASMA_desc *A, int *IPIV);
int PLASMA_sgetrf_incpiv_Tile(PLASMA_desc *A, PLASMA_desc *L, int *IPIV);
int PLASMA_sgetrs_Tile(PLASMA_enum trans, PLASMA_desc *A, int *IPIV, PLASMA_desc *B);
int PLASMA_sgetrs_incpiv_Tile(PLASMA_desc *A, PLASMA_desc *L, int *IPIV, PLASMA_desc *B);
#ifdef COMPLEX
int PLASMA_ssymm_Tile(PLASMA_enum side, PLASMA_enum uplo, float alpha, PLASMA_desc *A, PLASMA_desc *B, float beta, PLASMA_desc *C);
int PLASMA_ssyrk_Tile(PLASMA_enum uplo, PLASMA_enum trans, float alpha, PLASMA_desc *A, float beta, PLASMA_desc *C);
int PLASMA_ssyr2k_Tile(PLASMA_enum uplo, PLASMA_enum trans, float alpha, PLASMA_desc *A, PLASMA_desc *B, float beta, PLASMA_desc *C);
#endif
int PLASMA_ssyev_Tile(PLASMA_enum jobz, PLASMA_enum uplo, PLASMA_desc *A, float *W, PLASMA_desc *T, PLASMA_desc *Q);
int PLASMA_ssygv_Tile(PLASMA_enum itype, PLASMA_enum jobz, PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B, float *W, PLASMA_desc *T, PLASMA_desc *Q);
int PLASMA_ssygst_Tile(PLASMA_enum itype, PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B);
int PLASMA_ssytrd_Tile(PLASMA_enum jobz, PLASMA_enum uplo, PLASMA_desc *A, float *D, float *E, PLASMA_desc *T, PLASMA_desc *Q);
float PLASMA_slange_Tile(PLASMA_enum norm, PLASMA_desc *A, float *work);
#ifdef COMPLEX
float PLASMA_slansy_Tile(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc *A, float *work);
#endif
float PLASMA_slansy_Tile(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc *A, float *work);
int PLASMA_slaswp_Tile(PLASMA_desc *A, int K1, int K2, int *IPIV, int INCX);
int PLASMA_slauum_Tile(PLASMA_enum uplo, PLASMA_desc *A);
#ifdef COMPLEX
int PLASMA_splgsy_Tile(float bump, PLASMA_desc *A, unsigned long long int seed );
#endif
int PLASMA_splgsy_Tile(float bump, PLASMA_desc *A, unsigned long long int seed );
int PLASMA_splrnt_Tile(PLASMA_desc *A, unsigned long long int seed );
int PLASMA_sposv_Tile(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B);
int PLASMA_spotrf_Tile(PLASMA_enum uplo, PLASMA_desc *A);
int PLASMA_spotri_Tile(PLASMA_enum uplo, PLASMA_desc *A);
int PLASMA_spotrs_Tile(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B);
int PLASMA_ssymm_Tile(PLASMA_enum side, PLASMA_enum uplo, float alpha, PLASMA_desc *A, PLASMA_desc *B, float beta, PLASMA_desc *C);
int PLASMA_ssyrk_Tile(PLASMA_enum uplo, PLASMA_enum trans, float alpha, PLASMA_desc *A, float beta, PLASMA_desc *C);
int PLASMA_ssyr2k_Tile(PLASMA_enum uplo, PLASMA_enum trans, float alpha, PLASMA_desc *A, PLASMA_desc *B, float beta, PLASMA_desc *C);
int PLASMA_strmm_Tile(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, float alpha, PLASMA_desc *A, PLASMA_desc *B);
int PLASMA_strsm_Tile(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, float alpha, PLASMA_desc *A, PLASMA_desc *B);
int PLASMA_strsmpl_Tile(PLASMA_desc *A, PLASMA_desc *L, int *IPIV, PLASMA_desc *B);
int PLASMA_strtri_Tile(PLASMA_enum uplo, PLASMA_enum diag, PLASMA_desc *A);
int PLASMA_sorgbr_Tile(PLASMA_enum size, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *Q);
int PLASMA_sorghr_Tile(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *Q);
int PLASMA_sorglq_Tile(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B);
int PLASMA_sorgqr_Tile(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B);
int PLASMA_sorgtr_Tile(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B);
int PLASMA_sormlq_Tile(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B);
int PLASMA_sormqr_Tile(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B);

/** ****************************************************************************
 *  Declarations of math functions (tile layout, asynchronous execution) - alphabetical order
 **/
int PLASMA_sgebrd_Tile_Async(PLASMA_enum jobu, PLASMA_enum jobvt, PLASMA_desc *A, float *D, float *E, PLASMA_desc *U, PLASMA_desc *VT, PLASMA_desc *T, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_sgeev_Tile_Async(PLASMA_enum jobvl, PLASMA_enum jobvr, PLASMA_desc *A, float *W, PLASMA_desc *VL, PLASMA_desc *VR, PLASMA_desc *T, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_sgehrd_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_sgelqf_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_sgelqs_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_sgels_Tile_Async(PLASMA_enum trans, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_sgemm_Tile_Async(PLASMA_enum transA, PLASMA_enum transB, float alpha, PLASMA_desc *A, PLASMA_desc *B, float beta, PLASMA_desc *C, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_sgeqrf_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_sgeqrs_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_sgesv_Tile_Async(PLASMA_desc *A, int *IPIV, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_sgesv_incpiv_Tile_Async(PLASMA_desc *A, PLASMA_desc *L, int *IPIV, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_sgesvd_Tile_Async(PLASMA_enum jobu, PLASMA_enum jobvt, PLASMA_desc *A, float *S, PLASMA_desc *U, PLASMA_desc *VT, PLASMA_desc *T, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_sgetrf_Tile_Async(PLASMA_desc *A, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_sgetrf_incpiv_Tile_Async(PLASMA_desc *A, PLASMA_desc *L, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_sgetrs_Tile_Async(PLASMA_enum trans, PLASMA_desc *A, int *IPIV, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_sgetrs_incpiv_Tile_Async(PLASMA_desc *A, PLASMA_desc *L, int *IPIV, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
#ifdef COMPLEX
int PLASMA_ssymm_Tile_Async(PLASMA_enum side, PLASMA_enum uplo, float alpha, PLASMA_desc *A, PLASMA_desc *B, float beta, PLASMA_desc *C, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_ssyrk_Tile_Async(PLASMA_enum uplo, PLASMA_enum trans, float alpha, PLASMA_desc *A, float beta, PLASMA_desc *C, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_ssyr2k_Tile_Async(PLASMA_enum uplo, PLASMA_enum trans, float alpha, PLASMA_desc *A, PLASMA_desc *B, float beta, PLASMA_desc *C, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
int PLASMA_ssyev_Tile_Async(PLASMA_enum jobz, PLASMA_enum uplo, PLASMA_desc *A, float *W, PLASMA_desc *T, PLASMA_desc *Q, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_ssygv_Tile_Async(PLASMA_enum itype, PLASMA_enum jobz, PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B, float *W, PLASMA_desc *T, PLASMA_desc *Q, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_ssygst_Tile_Async(PLASMA_enum itype, PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_ssytrd_Tile_Async(PLASMA_enum jobz, PLASMA_enum uplo, PLASMA_desc *A, float *D, float *E, PLASMA_desc *T, PLASMA_desc *Q, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_slange_Tile_Async(PLASMA_enum norm, PLASMA_desc *A, float *work, float *value, PLASMA_sequence *sequence, PLASMA_request *request);
#ifdef COMPLEX
int PLASMA_slansy_Tile_Async(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc *A, float *work, float *value, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
int PLASMA_slansy_Tile_Async(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc *A, float *work, float *value, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_slaswp_Tile_Async(PLASMA_desc *A, int K1, int K2, int *IPIV, int INCX, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_slauum_Tile_Async(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_sequence *sequence, PLASMA_request *request);
#ifdef COMPLEX
int PLASMA_splgsy_Tile_Async(float bump, PLASMA_desc *A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
#endif
int PLASMA_splgsy_Tile_Async(float bump, PLASMA_desc *A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
int PLASMA_splrnt_Tile_Async(PLASMA_desc *A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
int PLASMA_sposv_Tile_Async(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_spotrf_Tile_Async(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_spotri_Tile_Async(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_spotrs_Tile_Async(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_ssymm_Tile_Async(PLASMA_enum side, PLASMA_enum uplo, float alpha, PLASMA_desc *A, PLASMA_desc *B, float beta, PLASMA_desc *C, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_ssyrk_Tile_Async(PLASMA_enum uplo, PLASMA_enum trans, float alpha, PLASMA_desc *A, float beta, PLASMA_desc *C, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_ssyr2k_Tile_Async(PLASMA_enum uplo, PLASMA_enum trans, float alpha, PLASMA_desc *A, PLASMA_desc *B, float beta, PLASMA_desc *C, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_strmm_Tile_Async(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, float alpha, PLASMA_desc *A, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_strsm_Tile_Async(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, float alpha, PLASMA_desc *A, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_strsmpl_Tile_Async(PLASMA_desc *A, PLASMA_desc *L, int *IPIV, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_strtri_Tile_Async(PLASMA_enum uplo, PLASMA_enum diag, PLASMA_desc *A, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_sorgbr_Tile_Async(PLASMA_enum side, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *Q, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_sorghr_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *Q, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_sorglq_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_sorgqr_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_sorgtr_Tile_Async(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_sormlq_Tile_Async(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_sormqr_Tile_Async(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);

int PLASMA_sgecfi_Async(int m, int n, float *A, PLASMA_enum f_in, int imb, int inb, PLASMA_enum f_out, int omb, int onb, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_sgetmi_Async(int m, int n, float *A, PLASMA_enum f_in, int mb,  int inb, PLASMA_sequence *sequence, PLASMA_request *request);

/** ****************************************************************************
 *  Declarations of workspace allocation functions (tile layout) - alphabetical order
 **/
int PLASMA_Alloc_Workspace_sgelqf(int M, int N, float **T);
int PLASMA_Alloc_Workspace_sgels( int M, int N, float **T);
int PLASMA_Alloc_Workspace_sgeqrf(int M, int N, float **T);
int PLASMA_Alloc_Workspace_sgesv_incpiv(        int N, float **L, int **IPIV);
int PLASMA_Alloc_Workspace_sgetrf_incpiv(int M, int N, float **L, int **IPIV);

int PLASMA_Alloc_Workspace_sgebrd(int M, int N, PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_sgeev( int N,        PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_sgehrd(int N,        PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_sgesvd(int M, int N, PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_ssyev( int M, int N, PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_ssygv( int M, int N, PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_ssytrd(int M, int N, PLASMA_desc **descT);

/** ****************************************************************************
 *  Declarations of workspace allocation functions (tile layout, asynchronous execution) - alphabetical order
 **/
int PLASMA_Alloc_Workspace_sgelqf_Tile(int M, int N, PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_sgels_Tile( int M, int N, PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_sgeqrf_Tile(int M, int N, PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_sgesv_incpiv_Tile (int N, PLASMA_desc **descL, int **IPIV);
int PLASMA_Alloc_Workspace_sgetrf_incpiv_Tile(int N, PLASMA_desc **descL, int **IPIV);

/** ****************************************************************************
 *  Auxiliary function prototypes
 **/
int PLASMA_sLapack_to_Tile(float *Af77, int LDA, PLASMA_desc *A);
int PLASMA_sTile_to_Lapack(PLASMA_desc *A, float *Af77, int LDA);
int PLASMA_sLapack_to_Tile_Async(float *Af77, int LDA, PLASMA_desc *A, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_sTile_to_Lapack_Async(PLASMA_desc *A, float *Af77, int LDA, PLASMA_sequence *sequence, PLASMA_request *request);

#ifdef __cplusplus
}
#endif

#endif
