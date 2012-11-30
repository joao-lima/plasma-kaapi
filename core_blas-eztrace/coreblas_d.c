/**                                                               
 *                                                                
 * @file coreblas_d.c                                             
 *                                                                
 *  PLASMA core_blas tracing kernel                               
 *  PLASMA is a software package provided by Univ. of Tennessee,  
 *  Univ. of California Berkeley and Univ. of Colorado Denver     
 *                                                                
 *  This file provides the wrapper for each function of the       
 *  core_blas library which will generate an event before and     
 *  after the execution of the kernel.                            
 *  This file is automatically generated with convert2eztrace.pl  
 *  script.                                                       
 *                                                                
 * @version 2.4.2                                                 
 * @author Mathieu Faverge                                        
 * @date 2010-11-15                                               
 * @generated d Thu Sep 15 12:09:49 2011
 *                                                                
 **/                                                              
#include <eztrace.h>           
#include <ev_codes.h>          
#include "common.h"            
#include "coreblas_ev_codes.h" 
#include "coreblas_macros.h"   
#undef COMPLEX                    
#define REAL                

/*****************************************************************
 *        Core functions                                          
 */

FUNCTION_VOID( CORE_dasum, ASUM, void ,
          (int storev, int uplo, int M, int N, double *A, int lda, double *work),
          (storev, uplo, M, N, A, lda, work) )
FUNCTION_VOID( CORE_daxpy, AXPY, void ,
          (int M, int N, double alpha, double *A, int LDA, double *B, int LDB),
          (M, N, alpha, A, LDA, B, LDB) )
FUNCTION_VOID( CORE_dbrdalg, BRDALG, void ,
          (PLASMA_enum uplo, int N, int NB, PLASMA_desc *pA, double *V, double *TAU, int i, int j, int m, int grsiz),
          (uplo, N, NB, pA, V, TAU, i, j, m, grsiz) )
FUNCTION_TYPE( CORE_dgelqt, GELQT, int ,
          (int M, int N, int IB, double *A, int LDA, double *T, int LDT, double *TAU, double *WORK),
          (M, N, IB, A, LDA, T, LDT, TAU, WORK) )
FUNCTION_VOID( CORE_dgemm, GEMM, void ,
          (int transA, int transB, int M, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC),
          (transA, transB, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC) )
FUNCTION_TYPE( CORE_dgeqrt, GEQRT, int ,
          (int M, int N, int IB, double *A, int LDA, double *T, int LDT, double *TAU, double *WORK),
          (M, N, IB, A, LDA, T, LDT, TAU, WORK) )
FUNCTION_TYPE( CORE_dgessm, GESSM, int ,
          (int M, int N, int K, int IB, int *IPIV, double *L, int LDL, double *A, int LDA),
          (M, N, K, IB, IPIV, L, LDL, A, LDA) )
FUNCTION_TYPE( CORE_dgetrf, GETRF, int ,
          (int m, int n, double *A, int lda, int *IPIV, int *info),
          (m, n, A, lda, IPIV, info) )
FUNCTION_TYPE( CORE_dgetrf_incpiv, GETRF, int ,
          (int M, int N, int IB, double *A, int LDA, int *IPIV, int *INFO),
          (M, N, IB, A, LDA, IPIV, INFO) )
FUNCTION_TYPE( CORE_dgetrf_reclap, GETRF, int ,
          (int M, int N, double *A, int LDA, int *IPIV, int *info),
          (M, N, A, LDA, IPIV, info) )
FUNCTION_TYPE( CORE_dgetrf_rectil, GETRF, int ,
          (const PLASMA_desc A, int *IPIV, int *info),
          (A, IPIV, info) )
FUNCTION_VOID( CORE_dgetrip, GETRIP, void ,
          (int m, int n, double *A, double *W) ,
          (m, n, A, W)  )
FUNCTION_VOID( CORE_dsygst, HEGST, void ,
          (int itype, PLASMA_enum uplo, int N, double *A, int LDA, double *B, int LDB, int *INFO),
          (itype, uplo, N, A, LDA, B, LDB, INFO) )
#ifdef COMPLEX
FUNCTION_VOID( CORE_dsymm, HEMM, void ,
          (int side, int uplo, int M, int N, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC),
          (side, uplo, M, N, alpha, A, LDA, B, LDB, beta, C, LDC) )
#endif
#ifdef COMPLEX
FUNCTION_VOID( CORE_dsyr2k, HER2K, void ,
          (int uplo, int trans, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC),
          (uplo, trans, N, K, alpha, A, LDA, B, LDB, beta, C, LDC) )
#endif
FUNCTION_TYPE( CORE_dsyrfb, HERFB, int ,
          ( PLASMA_enum uplo, int n, int k, int ib, int nb, double *A, int lda, double *T, int ldt, double *C, int ldc, double *WORK, int ldwork ),
          (uplo, n, k, ib, nb, A, lda, T, ldt, C, ldc, WORK, ldwork) )
#ifdef COMPLEX
FUNCTION_VOID( CORE_dsyrk, HERK, void ,
          (int uplo, int trans, int N, int K, double alpha, double *A, int LDA, double beta, double *C, int LDC),
          (uplo, trans, N, K, alpha, A, LDA, beta, C, LDC) )
#endif
FUNCTION_VOID( CORE_dlacpy, LACPY, void ,
          (PLASMA_enum uplo, int M, int N, double *A, int LDA, double *B, int LDB),
          (uplo, M, N, A, LDA, B, LDB) )
FUNCTION_VOID( CORE_dlange, LANGE, void ,
          (int norm, int M, int N, double *A, int LDA, double *work, double *normA),
          (norm, M, N, A, LDA, work, normA) )
#ifdef COMPLEX
FUNCTION_VOID( CORE_dlansy, LANHE, void ,
          (int norm, int uplo, int N, double *A, int LDA, double *work, double *normA),
          (norm, uplo, N, A, LDA, work, normA) )
#endif
FUNCTION_VOID( CORE_dlansy, LANSY, void ,
          (int norm, int uplo, int N, double *A, int LDA, double *work, double *normA),
          (norm, uplo, N, A, LDA, work, normA) )
FUNCTION_VOID( CORE_dlaset2, LASET, void ,
          (PLASMA_enum uplo, int M, int N, double alpha, double *A, int LDA),
          (uplo, M, N, alpha, A, LDA) )
FUNCTION_VOID( CORE_dlaset, LASET, void ,
          (PLASMA_enum uplo, int M, int N, double alpha, double beta, double *A, int LDA),
          (uplo, M, N, alpha, beta, A, LDA) )
FUNCTION_VOID( CORE_dlaswp, LASWP, void ,
          (int N, double *A, int LDA, int I1, int I2, int *IPIV, int INC),
          (N, A, LDA, I1, I2, IPIV, INC) )
FUNCTION_TYPE( CORE_dlaswp_ontile, LASWP, int ,
          (PLASMA_desc descA, int i1, int i2, int *ipiv, int inc),
          (descA, i1, i2, ipiv, inc) )
FUNCTION_TYPE( CORE_dswptr_ontile, TRSM, int ,
          (PLASMA_desc descA, int i1, int i2, int *ipiv, int inc, double *Akk, int ldak),
          (descA, i1, i2, ipiv, inc, Akk, ldak) )
FUNCTION_VOID( CORE_dlauum, LAUUM, void ,
          (int uplo, int N, double *A, int LDA),
          (uplo, N, A, LDA) )
#ifdef COMPLEX
FUNCTION_VOID( CORE_dplgsy, PLGHE, void ,
          ( double bump, int m, int n, double *A, int lda, int bigM, int m0, int n0, unsigned long long int seed ),
          (bump, m, n, A, lda, bigM, m0, n0, seed) )
#endif
FUNCTION_VOID( CORE_dplgsy, PLGSY, void ,
          ( double bump, int m, int n, double *A, int lda, int bigM, int m0, int n0, unsigned long long int seed ),
          (bump, m, n, A, lda, bigM, m0, n0, seed) )
FUNCTION_VOID( CORE_dplrnt, PLRNT, void ,
          ( int m, int n, double *A, int lda, int bigM, int m0, int n0, unsigned long long int seed ),
          (m, n, A, lda, bigM, m0, n0, seed) )
FUNCTION_VOID( CORE_dpotrf, POTRF, void ,
          (int uplo, int N, double *A, int LDA, int *INFO),
          (uplo, N, A, LDA, INFO) )
FUNCTION_VOID( CORE_dshiftw, SHIFTW, void ,
          (int s, int cl, int m, int n, int L, double *A, double *W) ,
          (s, cl, m, n, L, A, W)  )
FUNCTION_VOID( CORE_dshift, SHIFT, void ,
          (int s, int m, int n, int L, double *A) ,
          (s, m, n, L, A)  )
FUNCTION_TYPE( CORE_dssssm, SSSSM, int ,
          (int M1, int N1, int M2, int N2, int K, int IB, double *A1, int LDA1, double *A2, int LDA2, double *L1, int LDL1, double *L2, int LDL2, int *IPIV),
          (M1, N1, M2, N2, K, IB, A1, LDA1, A2, LDA2, L1, LDL1, L2, LDL2, IPIV) )
FUNCTION_VOID( CORE_dswpab, SWPAB, void ,
          (int i, int n1, int n2, double *A, double *work) ,
          (i, n1, n2, A, work)  )
FUNCTION_VOID( CORE_dsymm, SYMM, void ,
          (int side, int uplo, int M, int N, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC),
          (side, uplo, M, N, alpha, A, LDA, B, LDB, beta, C, LDC) )
FUNCTION_VOID( CORE_dsyr2k, SYR2K, void ,
          (int uplo, int trans, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC),
          (uplo, trans, N, K, alpha, A, LDA, B, LDB, beta, C, LDC) )
FUNCTION_VOID( CORE_dsyrk, SYRK, void ,
          (int uplo, int trans, int N, int K, double alpha, double *A, int LDA, double beta, double *C, int LDC),
          (uplo, trans, N, K, alpha, A, LDA, beta, C, LDC) )
FUNCTION_VOID( CORE_dtrdalg, TRDALG, void ,
          (PLASMA_enum uplo, int N, int NB, PLASMA_desc *pA, double *V, double *TAU, int i, int j, int m, int grsiz),
          (uplo, N, NB, pA, V, TAU, i, j, m, grsiz) )
FUNCTION_VOID( CORE_dtrmm, TRMM, void ,
          (int side, int uplo, int transA, int diag, int M, int N, double alpha, double *A, int LDA, double *B, int LDB),
          (side, uplo, transA, diag, M, N, alpha, A, LDA, B, LDB) )
FUNCTION_VOID( CORE_dtrsm, TRSM, void ,
          (int side, int uplo, int transA, int diag, int M, int N, double alpha, double *A, int LDA, double *B, int LDB),
          (side, uplo, transA, diag, M, N, alpha, A, LDA, B, LDB) )
FUNCTION_VOID( CORE_dtrtri, TRTRI, void ,
          (int uplo, int diag, int N, double *A, int LDA, int *info),
          (uplo, diag, N, A, LDA, info) )
FUNCTION_TYPE( CORE_dtslqt, TSLQT, int ,
          (int M, int N, int IB, double *A1, int LDA1, double *A2, int LDA2, double *T, int LDT, double *TAU, double *WORK),
          (M, N, IB, A1, LDA1, A2, LDA2, T, LDT, TAU, WORK) )
FUNCTION_TYPE( CORE_dtsmlq, TSMLQ, int ,
          (int side, int trans, int M1, int N1, int M2, int N2, int K, int IB, double *A1, int LDA1, double *A2, int LDA2, double *V, int LDV, double *T, int LDT, double *WORK, int LDWORK),
          (side, trans, M1, N1, M2, N2, K, IB, A1, LDA1, A2, LDA2, V, LDV, T, LDT, WORK, LDWORK) )
FUNCTION_TYPE( CORE_dtsmlq_corner, TSMLQ, int ,
          ( int m1, int n1, int m2, int n2, int m3, int n3, int k, int ib, int nb, double *A1, int lda1, double *A2, int lda2, double *A3, int lda3, double *V, int ldv, double *T, int ldt, double *WORK, int ldwork),
          (m1, n1, m2, n2, m3, n3, k, ib, nb, A1, lda1, A2, lda2, A3, lda3, V, ldv, T, ldt, WORK, ldwork) )
FUNCTION_TYPE( CORE_dtsmlq_sytra1, TSMLQ, int ,
          ( int side, int trans, int m1, int n1, int m2, int n2, int k, int ib, double *A1, int lda1, double *A2, int lda2, double *V, int ldv, double *T, int ldt, double *WORK, int ldwork),
          (side, trans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork) )
FUNCTION_TYPE( CORE_dtsmqr, TSMQR, int ,
          (int side, int trans, int M1, int N1, int M2, int N2, int K, int IB, double *A1, int LDA1, double *A2, int LDA2, double *V, int LDV, double *T, int LDT, double *WORK, int LDWORK),
          (side, trans, M1, N1, M2, N2, K, IB, A1, LDA1, A2, LDA2, V, LDV, T, LDT, WORK, LDWORK) )
FUNCTION_TYPE( CORE_dtsmqr_corner, TSMQR, int ,
          ( int m1, int n1, int m2, int n2, int m3, int n3, int k, int ib, int nb, double *A1, int lda1, double *A2, int lda2, double *A3, int lda3, double *V, int ldv, double *T, int ldt, double *WORK, int ldwork),
          (m1, n1, m2, n2, m3, n3, k, ib, nb, A1, lda1, A2, lda2, A3, lda3, V, ldv, T, ldt, WORK, ldwork) )
FUNCTION_TYPE( CORE_dtsmqr_sytra1, TSMQR, int ,
          ( int side, int trans, int m1, int n1, int m2, int n2, int k, int ib, double *A1, int lda1, double *A2, int lda2, double *V, int ldv, double *T, int ldt, double *WORK, int ldwork),
          (side, trans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork) )
FUNCTION_TYPE( CORE_dtsqrt, TSQRT, int ,
          (int M, int N, int IB, double *A1, int LDA1, double *A2, int LDA2, double *T, int LDT, double *TAU, double *WORK),
          (M, N, IB, A1, LDA1, A2, LDA2, T, LDT, TAU, WORK) )
FUNCTION_TYPE( CORE_dtsrfb, TSRFB, int ,
          (int side, int trans, int direct, int storev, int M1, int N1, int M2, int N2, int K, double *A1, int LDA1, double *A2, int LDA2, double *V, int LDV, double *T, int LDT, double *WORK, int LDWORK),
          (side, trans, direct, storev, M1, N1, M2, N2, K, A1, LDA1, A2, LDA2, V, LDV, T, LDT, WORK, LDWORK) )
FUNCTION_TYPE( CORE_dtstrf, TSTRF, int ,
          (int M, int N, int IB, int NB, double *U, int LDU, double *A, int LDA, double *L, int LDL, int *IPIV, double *WORK, int LDWORK, int *INFO),
          (M, N, IB, NB, U, LDU, A, LDA, L, LDL, IPIV, WORK, LDWORK, INFO) )
FUNCTION_TYPE( CORE_dttlqt, TTLQT, int ,
          (int M, int N, int IB, double *A1, int LDA1, double *A2, int LDA2, double *T, int LDT, double *TAU, double *WORK),
          (M, N, IB, A1, LDA1, A2, LDA2, T, LDT, TAU, WORK) )
FUNCTION_TYPE( CORE_dttmlq, TTMLQ, int ,
          (int side, int trans, int M1, int N1, int M2, int N2, int K, int IB, double *A1, int LDA1, double *A2, int LDA2, double *V, int LDV, double *T, int LDT, double *WORK, int LDWORK),
          (side, trans, M1, N1, M2, N2, K, IB, A1, LDA1, A2, LDA2, V, LDV, T, LDT, WORK, LDWORK) )
FUNCTION_TYPE( CORE_dttmqr, TTMQR, int ,
          (int side, int trans, int M1, int N1, int M2, int N2, int KK, int IB, double *A1, int LDA1, double *A2, int LDA2, double *V, int LDV, double *T, int LDT, double *WORK, int LDWORK),
          (side, trans, M1, N1, M2, N2, KK, IB, A1, LDA1, A2, LDA2, V, LDV, T, LDT, WORK, LDWORK) )
FUNCTION_TYPE( CORE_dttqrt, TTQRT, int ,
          (int M, int N, int IB, double *A1, int LDA1, double *A2, int LDA2, double *T, int LDT, double *TAU, double *WORK),
          (M, N, IB, A1, LDA1, A2, LDA2, T, LDT, TAU, WORK) )
FUNCTION_TYPE( CORE_dttrfb, TTRFB, int ,
          (int side, int trans, int direct, int storev, int M1, int N1, int M2, int N2, int K, double *A1, int LDA1, double *A2, int LDA2, double *V, int LDV, double *T, int LDT, double *WORK, int LDWORK),
          (side, trans, direct, storev, M1, N1, M2, N2, K, A1, LDA1, A2, LDA2, V, LDV, T, LDT, WORK, LDWORK) )
FUNCTION_TYPE( CORE_dormlq, UNMLQ, int ,
          (int side, int trans, int M, int N, int K, int IB, double *A, int LDA, double *T, int LDT, double *C, int LDC, double *WORK, int LDWORK),
          (side, trans, M, N, K, IB, A, LDA, T, LDT, C, LDC, WORK, LDWORK) )
FUNCTION_TYPE( CORE_dormqr, UNMQR, int ,
          (int side, int trans, int M, int N, int K, int IB, double *A, int LDA, double *T, int LDT, double *C, int LDC, double *WORK, int LDWORK),
          (side, trans, M, N, K, IB, A, LDA, T, LDT, C, LDC, WORK, LDWORK) )

/*****************************************************************
 *        QUARK Wrapper functions                                 
 */

FUNCTION_QUARK( CORE_dasum_quark, ASUM )
FUNCTION_QUARK( CORE_dasum_f1_quark, ASUM )
FUNCTION_QUARK( CORE_daxpy_quark, AXPY )
FUNCTION_QUARK( CORE_dbrdalg_quark, BRDALG )
FUNCTION_QUARK( CORE_dgelqt_quark, GELQT )
FUNCTION_QUARK( CORE_dgemm_quark, GEMM )
FUNCTION_QUARK( CORE_dgemm_f2_quark, GEMM )
FUNCTION_QUARK( CORE_dgemm_p2_quark, GEMM )
FUNCTION_QUARK( CORE_dgemm_p3_quark, GEMM )
FUNCTION_QUARK( CORE_dgemm_p2f1_quark, GEMM )
FUNCTION_QUARK( CORE_dgeqrt_quark, GEQRT )
FUNCTION_QUARK( CORE_dgessm_quark, GESSM )
FUNCTION_QUARK( CORE_dgetrf_quark, GETRF )
FUNCTION_QUARK( CORE_dgetrf_incpiv_quark, GETRF )
FUNCTION_QUARK( CORE_dgetrf_reclap_quark, GETRF )
FUNCTION_QUARK( CORE_dgetrf_rectil_quark, GETRF )
FUNCTION_QUARK( CORE_dgetrip_quark, GETRIP )
FUNCTION_QUARK( CORE_dgetrip_f1_quark, GETRIP )
FUNCTION_QUARK( CORE_dgetrip_f2_quark, GETRIP )
FUNCTION_QUARK( CORE_dsygst_quark, HEGST )
#ifdef COMPLEX
FUNCTION_QUARK( CORE_dsymm_quark, HEMM )
#endif
#ifdef COMPLEX
FUNCTION_QUARK( CORE_dsyr2k_quark, HER2K )
#endif
FUNCTION_QUARK( CORE_dsyrfb_quark, HERFB )
#ifdef COMPLEX
FUNCTION_QUARK( CORE_dsyrk_quark, HERK )
#endif
FUNCTION_QUARK( CORE_dlacpy_quark, LACPY )
FUNCTION_QUARK( CORE_dlange_quark, LANGE )
FUNCTION_QUARK( CORE_dlange_f1_quark, LANGE )
#ifdef COMPLEX
FUNCTION_QUARK( CORE_dlansy_quark, LANHE )
FUNCTION_QUARK( CORE_dlansy_f1_quark, LANHE )
#endif
FUNCTION_QUARK( CORE_dlansy_quark, LANSY )
FUNCTION_QUARK( CORE_dlansy_f1_quark, LANSY )
FUNCTION_QUARK( CORE_dlaset2_quark, LASET )
FUNCTION_QUARK( CORE_dlaset_quark, LASET )
FUNCTION_QUARK( CORE_dlaswp_quark, LASWP )
FUNCTION_QUARK( CORE_dlaswp_f2_quark, LASWP )
FUNCTION_QUARK( CORE_dlaswp_ontile_quark, LASWP )
FUNCTION_QUARK( CORE_dlaswp_ontile_f2_quark, LASWP )
FUNCTION_QUARK( CORE_dswptr_ontile_quark, TRSM )
FUNCTION_QUARK( CORE_dlauum_quark, LAUUM )
#ifdef COMPLEX
FUNCTION_QUARK( CORE_dplgsy_quark, PLGHE )
#endif
FUNCTION_QUARK( CORE_dplgsy_quark, PLGSY )
FUNCTION_QUARK( CORE_dplrnt_quark, PLRNT )
FUNCTION_QUARK( CORE_dpotrf_quark, POTRF )
FUNCTION_QUARK( CORE_dshiftw_quark, SHIFTW )
FUNCTION_QUARK( CORE_dshift_quark, SHIFT )
FUNCTION_QUARK( CORE_dssssm_quark, SSSSM )
FUNCTION_QUARK( CORE_dswpab_quark, SWPAB )
FUNCTION_QUARK( CORE_dsymm_quark, SYMM )
FUNCTION_QUARK( CORE_dsyr2k_quark, SYR2K )
FUNCTION_QUARK( CORE_dsyrk_quark, SYRK )
FUNCTION_QUARK( CORE_dtrdalg_quark, TRDALG )
FUNCTION_QUARK( CORE_dtrmm_quark, TRMM )
FUNCTION_QUARK( CORE_dtrmm_p2_quark, TRMM )
FUNCTION_QUARK( CORE_dtrsm_quark, TRSM )
FUNCTION_QUARK( CORE_dtrtri_quark, TRTRI )
FUNCTION_QUARK( CORE_dtslqt_quark, TSLQT )
FUNCTION_QUARK( CORE_dtsmlq_quark, TSMLQ )
FUNCTION_QUARK( CORE_dtsmlq_corner_quark, TSMLQ )
FUNCTION_QUARK( CORE_dtsmlq_sytra1_quark, TSMLQ )
FUNCTION_QUARK( CORE_dtsmqr_quark, TSMQR )
FUNCTION_QUARK( CORE_dtsmqr_corner_quark, TSMQR )
FUNCTION_QUARK( CORE_dtsmqr_sytra1_quark, TSMQR )
FUNCTION_QUARK( CORE_dtsqrt_quark, TSQRT )
FUNCTION_QUARK( CORE_dtstrf_quark, TSTRF )
FUNCTION_QUARK( CORE_dttlqt_quark, TTLQT )
FUNCTION_QUARK( CORE_dttmlq_quark, TTMLQ )
FUNCTION_QUARK( CORE_dttmqr_quark, TTMQR )
FUNCTION_QUARK( CORE_dttqrt_quark, TTQRT )
FUNCTION_QUARK( CORE_dormlq_quark, UNMLQ )
FUNCTION_QUARK( CORE_dormqr_quark, UNMQR )

