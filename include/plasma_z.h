/**
 *
 * @file plasma_z.h
 *
 *  PLASMA header file for double _Complex routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Azzam Haidar
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#ifndef _PLASMA_Z_H_
#define _PLASMA_Z_H_
#undef REAL
#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

/** ****************************************************************************
 *  Declarations of math functions (LAPACK layout) - alphabetical order
 **/
int PLASMA_zgebrd(PLASMA_enum jobu, PLASMA_enum jobvt, int M, int N, PLASMA_Complex64_t *A, int LDA, double *D, double *E, PLASMA_Complex64_t *U, int LDU, PLASMA_Complex64_t *VT, int LDVT, PLASMA_desc *T);
int PLASMA_zgeev(PLASMA_enum jobvl, PLASMA_enum jobvr, int N, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *W, PLASMA_Complex64_t *VL, int LDVL, PLASMA_Complex64_t *VR, int LDVR, PLASMA_Complex64_t *T);
int PLASMA_zgehrd(int N, int ILO, int IHI, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *T);
int PLASMA_zgelqf(int M, int N, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *T);
int PLASMA_zgelqs(int M, int N, int NRHS, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *T, PLASMA_Complex64_t *B, int LDB);
int PLASMA_zgels(PLASMA_enum trans, int M, int N, int NRHS, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *T, PLASMA_Complex64_t *B, int LDB);
int PLASMA_zgemm(PLASMA_enum transA, PLASMA_enum transB, int M, int N, int K, PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB, PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int LDC);
int PLASMA_zgeqrf(int M, int N, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *T);
int PLASMA_zgeqrs(int M, int N, int NRHS, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *T, PLASMA_Complex64_t *B, int LDB);
int PLASMA_zgesv(int N, int NRHS, PLASMA_Complex64_t *A, int LDA, int *IPIV, PLASMA_Complex64_t *B, int LDB);
int PLASMA_zgesv_incpiv(int N, int NRHS, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *L, int *IPIV, PLASMA_Complex64_t *B, int LDB);
int PLASMA_zgesvd(PLASMA_enum jobu, PLASMA_enum jobvt, int M, int N, PLASMA_Complex64_t *A, int LDA, double *S, PLASMA_Complex64_t *U, int LDU, PLASMA_Complex64_t *VT, int LDVT, PLASMA_desc *T);
int PLASMA_zgetrf(int M, int N, PLASMA_Complex64_t *A, int LDA, int *IPIV);
int PLASMA_zgetrf_incpiv(int M, int N, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *L, int *IPIV);
int PLASMA_zgetrs(PLASMA_enum trans, int N, int NRHS, PLASMA_Complex64_t *A, int LDA, int *IPIV, PLASMA_Complex64_t *B, int LDB);
int PLASMA_zgetrs_incpiv(PLASMA_enum trans, int N, int NRHS, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *L, int *IPIV, PLASMA_Complex64_t *B, int LDB);
#ifdef COMPLEX
int PLASMA_zhemm(PLASMA_enum side, PLASMA_enum uplo, int M, int N, PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB, PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int LDC);
int PLASMA_zherk(PLASMA_enum uplo, PLASMA_enum trans, int N, int K, double alpha, PLASMA_Complex64_t *A, int LDA, double beta, PLASMA_Complex64_t *C, int LDC);
int PLASMA_zher2k(PLASMA_enum uplo, PLASMA_enum trans, int N, int K, PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB, double beta, PLASMA_Complex64_t *C, int LDC);
#endif
int PLASMA_zheev(PLASMA_enum jobz, PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA, double *W, PLASMA_desc *T, PLASMA_Complex64_t *Q, int LDQ);
int PLASMA_zhegv(PLASMA_enum itype, PLASMA_enum jobz, PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB, double *W, PLASMA_desc *T, PLASMA_Complex64_t *Q, int LDQ);
int PLASMA_zhegst(PLASMA_enum itype, PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB);
int PLASMA_zhetrd(PLASMA_enum jobz, PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA, double *D, double *E, PLASMA_desc *T, PLASMA_Complex64_t *Q, int LDQ);
double PLASMA_zlange(PLASMA_enum norm, int M, int N, PLASMA_Complex64_t *A, int LDA, double *work);
#ifdef COMPLEX
double PLASMA_zlanhe(PLASMA_enum norm, PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA, double *work);
#endif
double PLASMA_zlansy(PLASMA_enum norm, PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA, double *work);
int PLASMA_zlaswp(int N, PLASMA_Complex64_t *A, int LDA, int K1, int K2, int *IPIV, int INCX);
int PLASMA_zlauum(PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA);
#ifdef COMPLEX
int PLASMA_zplghe( double bump, int N, PLASMA_Complex64_t *A, int LDA, unsigned long long int seed );
#endif
int PLASMA_zplgsy( PLASMA_Complex64_t bump, int N, PLASMA_Complex64_t *A, int LDA, unsigned long long int seed );
int PLASMA_zplrnt( int M, int N, PLASMA_Complex64_t *A, int LDA, unsigned long long int seed );
int PLASMA_zposv(PLASMA_enum uplo, int N, int NRHS, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB);
int PLASMA_zpotrf(PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA);
int PLASMA_zpotri(PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA);
int PLASMA_zpotrs(PLASMA_enum uplo, int N, int NRHS, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB);
int PLASMA_zsymm(PLASMA_enum side, PLASMA_enum uplo, int M, int N, PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB, PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int LDC);
int PLASMA_zsyrk(PLASMA_enum uplo, PLASMA_enum trans, int N, int K, PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int LDC);
int PLASMA_zsyr2k(PLASMA_enum uplo, PLASMA_enum trans, int N, int K, PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB, PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int LDC);
int PLASMA_ztrmm(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, int N, int NRHS, PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB);
int PLASMA_ztrsm(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, int N, int NRHS, PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB);
int PLASMA_ztrsmpl(int N, int NRHS, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *L, int *IPIV, PLASMA_Complex64_t *B, int LDB);
int PLASMA_ztrtri(PLASMA_enum uplo, PLASMA_enum diag, int N, PLASMA_Complex64_t *A, int LDA);
int PLASMA_zungbr(PLASMA_enum side, int M, int N, int K, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *T, PLASMA_Complex64_t *Q, int LDQ);
int PLASMA_zunghr(int N, int ILO, int IHI, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *T, PLASMA_Complex64_t *Q, int LDQ);
int PLASMA_zunglq(int M, int N, int K, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *T, PLASMA_Complex64_t *B, int LDB);
int PLASMA_zungqr(int M, int N, int K, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *T, PLASMA_Complex64_t *B, int LDB);
int PLASMA_zungtr(PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *T, PLASMA_Complex64_t *B, int LDB);
int PLASMA_zunmlq(PLASMA_enum side, PLASMA_enum trans, int M, int N, int K, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *T, PLASMA_Complex64_t *B, int LDB);
int PLASMA_zunmqr(PLASMA_enum side, PLASMA_enum trans, int M, int N, int K, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *T, PLASMA_Complex64_t *B, int LDB);

int PLASMA_zgecfi(int m, int n, PLASMA_Complex64_t *A, PLASMA_enum fin, int imb, int inb, PLASMA_enum fout, int omb, int onb);
int PLASMA_zgetmi(int m, int n, PLASMA_Complex64_t *A, PLASMA_enum fin, int mb,  int nb);

/** ****************************************************************************
 *  Declarations of math functions (tile layout) - alphabetical order
 **/
int PLASMA_zgebrd_Tile(PLASMA_enum jobu, PLASMA_enum jobvt, PLASMA_desc *A, double *D, double *E, PLASMA_desc *U, PLASMA_desc *VT, PLASMA_desc *T);
int PLASMA_zgeev_Tile(PLASMA_enum jobvl, PLASMA_enum jobvr, PLASMA_desc *A, PLASMA_Complex64_t *W, PLASMA_desc *VL, PLASMA_desc *VR, PLASMA_desc *T);
int PLASMA_zgehrd_Tile(PLASMA_desc *A, PLASMA_desc *T);
int PLASMA_zgelqf_Tile(PLASMA_desc *A, PLASMA_desc *T);
int PLASMA_zgelqs_Tile(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B);
int PLASMA_zgels_Tile(PLASMA_enum trans, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B);
int PLASMA_zgemm_Tile(PLASMA_enum transA, PLASMA_enum transB, PLASMA_Complex64_t alpha, PLASMA_desc *A, PLASMA_desc *B, PLASMA_Complex64_t beta, PLASMA_desc *C);
int PLASMA_zgeqrf_Tile(PLASMA_desc *A, PLASMA_desc *T);
int PLASMA_zgeqrs_Tile(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B);
int PLASMA_zgesv_Tile(PLASMA_desc *A, int *IPIV, PLASMA_desc *B);
int PLASMA_zgesv_incpiv_Tile(PLASMA_desc *A, PLASMA_desc *L, int *IPIV, PLASMA_desc *B);
int PLASMA_zgesvd_Tile(PLASMA_enum jobu, PLASMA_enum jobvt, PLASMA_desc *A, double *S, PLASMA_desc *U, PLASMA_desc *VT, PLASMA_desc *T);
int PLASMA_zgetrf_Tile(PLASMA_desc *A, int *IPIV);
int PLASMA_zgetrf_incpiv_Tile(PLASMA_desc *A, PLASMA_desc *L, int *IPIV);
int PLASMA_zgetrs_Tile(PLASMA_enum trans, PLASMA_desc *A, int *IPIV, PLASMA_desc *B);
int PLASMA_zgetrs_incpiv_Tile(PLASMA_desc *A, PLASMA_desc *L, int *IPIV, PLASMA_desc *B);
#ifdef COMPLEX
int PLASMA_zhemm_Tile(PLASMA_enum side, PLASMA_enum uplo, PLASMA_Complex64_t alpha, PLASMA_desc *A, PLASMA_desc *B, PLASMA_Complex64_t beta, PLASMA_desc *C);
int PLASMA_zherk_Tile(PLASMA_enum uplo, PLASMA_enum trans, double alpha, PLASMA_desc *A, double beta, PLASMA_desc *C);
int PLASMA_zher2k_Tile(PLASMA_enum uplo, PLASMA_enum trans, PLASMA_Complex64_t alpha, PLASMA_desc *A, PLASMA_desc *B, double beta, PLASMA_desc *C);
#endif
int PLASMA_zheev_Tile(PLASMA_enum jobz, PLASMA_enum uplo, PLASMA_desc *A, double *W, PLASMA_desc *T, PLASMA_desc *Q);
int PLASMA_zhegv_Tile(PLASMA_enum itype, PLASMA_enum jobz, PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B, double *W, PLASMA_desc *T, PLASMA_desc *Q);
int PLASMA_zhegst_Tile(PLASMA_enum itype, PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B);
int PLASMA_zhetrd_Tile(PLASMA_enum jobz, PLASMA_enum uplo, PLASMA_desc *A, double *D, double *E, PLASMA_desc *T, PLASMA_desc *Q);
double PLASMA_zlange_Tile(PLASMA_enum norm, PLASMA_desc *A, double *work);
#ifdef COMPLEX
double PLASMA_zlanhe_Tile(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc *A, double *work);
#endif
double PLASMA_zlansy_Tile(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc *A, double *work);
int PLASMA_zlaswp_Tile(PLASMA_desc *A, int K1, int K2, int *IPIV, int INCX);
int PLASMA_zlauum_Tile(PLASMA_enum uplo, PLASMA_desc *A);
#ifdef COMPLEX
int PLASMA_zplghe_Tile(double bump, PLASMA_desc *A, unsigned long long int seed );
#endif
int PLASMA_zplgsy_Tile(PLASMA_Complex64_t bump, PLASMA_desc *A, unsigned long long int seed );
int PLASMA_zplrnt_Tile(PLASMA_desc *A, unsigned long long int seed );
int PLASMA_zposv_Tile(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B);
int PLASMA_zpotrf_Tile(PLASMA_enum uplo, PLASMA_desc *A);
int PLASMA_zpotri_Tile(PLASMA_enum uplo, PLASMA_desc *A);
int PLASMA_zpotrs_Tile(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B);
int PLASMA_zsymm_Tile(PLASMA_enum side, PLASMA_enum uplo, PLASMA_Complex64_t alpha, PLASMA_desc *A, PLASMA_desc *B, PLASMA_Complex64_t beta, PLASMA_desc *C);
int PLASMA_zsyrk_Tile(PLASMA_enum uplo, PLASMA_enum trans, PLASMA_Complex64_t alpha, PLASMA_desc *A, PLASMA_Complex64_t beta, PLASMA_desc *C);
int PLASMA_zsyr2k_Tile(PLASMA_enum uplo, PLASMA_enum trans, PLASMA_Complex64_t alpha, PLASMA_desc *A, PLASMA_desc *B, PLASMA_Complex64_t beta, PLASMA_desc *C);
int PLASMA_ztrmm_Tile(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, PLASMA_Complex64_t alpha, PLASMA_desc *A, PLASMA_desc *B);
int PLASMA_ztrsm_Tile(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, PLASMA_Complex64_t alpha, PLASMA_desc *A, PLASMA_desc *B);
int PLASMA_ztrsmpl_Tile(PLASMA_desc *A, PLASMA_desc *L, int *IPIV, PLASMA_desc *B);
int PLASMA_ztrtri_Tile(PLASMA_enum uplo, PLASMA_enum diag, PLASMA_desc *A);
int PLASMA_zungbr_Tile(PLASMA_enum size, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *Q);
int PLASMA_zunghr_Tile(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *Q);
int PLASMA_zunglq_Tile(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B);
int PLASMA_zungqr_Tile(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B);
int PLASMA_zungtr_Tile(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B);
int PLASMA_zunmlq_Tile(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B);
int PLASMA_zunmqr_Tile(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B);

/** ****************************************************************************
 *  Declarations of math functions (tile layout, asynchronous execution) - alphabetical order
 **/
int PLASMA_zgebrd_Tile_Async(PLASMA_enum jobu, PLASMA_enum jobvt, PLASMA_desc *A, double *D, double *E, PLASMA_desc *U, PLASMA_desc *VT, PLASMA_desc *T, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zgeev_Tile_Async(PLASMA_enum jobvl, PLASMA_enum jobvr, PLASMA_desc *A, PLASMA_Complex64_t *W, PLASMA_desc *VL, PLASMA_desc *VR, PLASMA_desc *T, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zgehrd_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zgelqf_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zgelqs_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zgels_Tile_Async(PLASMA_enum trans, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zgemm_Tile_Async(PLASMA_enum transA, PLASMA_enum transB, PLASMA_Complex64_t alpha, PLASMA_desc *A, PLASMA_desc *B, PLASMA_Complex64_t beta, PLASMA_desc *C, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zgeqrf_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zgeqrs_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zgesv_Tile_Async(PLASMA_desc *A, int *IPIV, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zgesv_incpiv_Tile_Async(PLASMA_desc *A, PLASMA_desc *L, int *IPIV, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zgesvd_Tile_Async(PLASMA_enum jobu, PLASMA_enum jobvt, PLASMA_desc *A, double *S, PLASMA_desc *U, PLASMA_desc *VT, PLASMA_desc *T, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zgetrf_Tile_Async(PLASMA_desc *A, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zgetrf_incpiv_Tile_Async(PLASMA_desc *A, PLASMA_desc *L, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zgetrs_Tile_Async(PLASMA_enum trans, PLASMA_desc *A, int *IPIV, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zgetrs_incpiv_Tile_Async(PLASMA_desc *A, PLASMA_desc *L, int *IPIV, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
#ifdef COMPLEX
int PLASMA_zhemm_Tile_Async(PLASMA_enum side, PLASMA_enum uplo, PLASMA_Complex64_t alpha, PLASMA_desc *A, PLASMA_desc *B, PLASMA_Complex64_t beta, PLASMA_desc *C, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zherk_Tile_Async(PLASMA_enum uplo, PLASMA_enum trans, double alpha, PLASMA_desc *A, double beta, PLASMA_desc *C, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zher2k_Tile_Async(PLASMA_enum uplo, PLASMA_enum trans, PLASMA_Complex64_t alpha, PLASMA_desc *A, PLASMA_desc *B, double beta, PLASMA_desc *C, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
int PLASMA_zheev_Tile_Async(PLASMA_enum jobz, PLASMA_enum uplo, PLASMA_desc *A, double *W, PLASMA_desc *T, PLASMA_desc *Q, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zhegv_Tile_Async(PLASMA_enum itype, PLASMA_enum jobz, PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B, double *W, PLASMA_desc *T, PLASMA_desc *Q, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zhegst_Tile_Async(PLASMA_enum itype, PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zhetrd_Tile_Async(PLASMA_enum jobz, PLASMA_enum uplo, PLASMA_desc *A, double *D, double *E, PLASMA_desc *T, PLASMA_desc *Q, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zlange_Tile_Async(PLASMA_enum norm, PLASMA_desc *A, double *work, double *value, PLASMA_sequence *sequence, PLASMA_request *request);
#ifdef COMPLEX
int PLASMA_zlanhe_Tile_Async(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc *A, double *work, double *value, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
int PLASMA_zlansy_Tile_Async(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc *A, double *work, double *value, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zlaswp_Tile_Async(PLASMA_desc *A, int K1, int K2, int *IPIV, int INCX, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zlauum_Tile_Async(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_sequence *sequence, PLASMA_request *request);
#ifdef COMPLEX
int PLASMA_zplghe_Tile_Async(double bump, PLASMA_desc *A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
#endif
int PLASMA_zplgsy_Tile_Async(PLASMA_Complex64_t bump, PLASMA_desc *A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
int PLASMA_zplrnt_Tile_Async(PLASMA_desc *A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
int PLASMA_zposv_Tile_Async(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zpotrf_Tile_Async(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zpotri_Tile_Async(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zpotrs_Tile_Async(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zsymm_Tile_Async(PLASMA_enum side, PLASMA_enum uplo, PLASMA_Complex64_t alpha, PLASMA_desc *A, PLASMA_desc *B, PLASMA_Complex64_t beta, PLASMA_desc *C, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zsyrk_Tile_Async(PLASMA_enum uplo, PLASMA_enum trans, PLASMA_Complex64_t alpha, PLASMA_desc *A, PLASMA_Complex64_t beta, PLASMA_desc *C, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zsyr2k_Tile_Async(PLASMA_enum uplo, PLASMA_enum trans, PLASMA_Complex64_t alpha, PLASMA_desc *A, PLASMA_desc *B, PLASMA_Complex64_t beta, PLASMA_desc *C, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_ztrmm_Tile_Async(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, PLASMA_Complex64_t alpha, PLASMA_desc *A, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_ztrsm_Tile_Async(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, PLASMA_Complex64_t alpha, PLASMA_desc *A, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_ztrsmpl_Tile_Async(PLASMA_desc *A, PLASMA_desc *L, int *IPIV, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_ztrtri_Tile_Async(PLASMA_enum uplo, PLASMA_enum diag, PLASMA_desc *A, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zungbr_Tile_Async(PLASMA_enum side, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *Q, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zunghr_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *Q, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zunglq_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zungqr_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zungtr_Tile_Async(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zunmlq_Tile_Async(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zunmqr_Tile_Async(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B, PLASMA_sequence *sequence, PLASMA_request *request);

int PLASMA_zgecfi_Async(int m, int n, PLASMA_Complex64_t *A, PLASMA_enum f_in, int imb, int inb, PLASMA_enum f_out, int omb, int onb, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zgetmi_Async(int m, int n, PLASMA_Complex64_t *A, PLASMA_enum f_in, int mb,  int inb, PLASMA_sequence *sequence, PLASMA_request *request);

/** ****************************************************************************
 *  Declarations of workspace allocation functions (tile layout) - alphabetical order
 **/
int PLASMA_Alloc_Workspace_zgelqf(int M, int N, PLASMA_Complex64_t **T);
int PLASMA_Alloc_Workspace_zgels( int M, int N, PLASMA_Complex64_t **T);
int PLASMA_Alloc_Workspace_zgeqrf(int M, int N, PLASMA_Complex64_t **T);
int PLASMA_Alloc_Workspace_zgesv_incpiv(        int N, PLASMA_Complex64_t **L, int **IPIV);
int PLASMA_Alloc_Workspace_zgetrf_incpiv(int M, int N, PLASMA_Complex64_t **L, int **IPIV);

int PLASMA_Alloc_Workspace_zgebrd(int M, int N, PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_zgeev( int N,        PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_zgehrd(int N,        PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_zgesvd(int M, int N, PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_zheev( int M, int N, PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_zhegv( int M, int N, PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_zhetrd(int M, int N, PLASMA_desc **descT);

/** ****************************************************************************
 *  Declarations of workspace allocation functions (tile layout, asynchronous execution) - alphabetical order
 **/
int PLASMA_Alloc_Workspace_zgelqf_Tile(int M, int N, PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_zgels_Tile( int M, int N, PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_zgeqrf_Tile(int M, int N, PLASMA_desc **descT);
int PLASMA_Alloc_Workspace_zgesv_incpiv_Tile (int N, PLASMA_desc **descL, int **IPIV);
int PLASMA_Alloc_Workspace_zgetrf_incpiv_Tile(int N, PLASMA_desc **descL, int **IPIV);

/** ****************************************************************************
 *  Auxiliary function prototypes
 **/
int PLASMA_zLapack_to_Tile(PLASMA_Complex64_t *Af77, int LDA, PLASMA_desc *A);
int PLASMA_zTile_to_Lapack(PLASMA_desc *A, PLASMA_Complex64_t *Af77, int LDA);
int PLASMA_zLapack_to_Tile_Async(PLASMA_Complex64_t *Af77, int LDA, PLASMA_desc *A, PLASMA_sequence *sequence, PLASMA_request *request);
int PLASMA_zTile_to_Lapack_Async(PLASMA_desc *A, PLASMA_Complex64_t *Af77, int LDA, PLASMA_sequence *sequence, PLASMA_request *request);

#ifdef __cplusplus
}
#endif

#endif
