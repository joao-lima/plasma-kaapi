/**
 *
 * @file core_dblas.h
 *
 *  PLASMA auxiliary routines
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
#ifndef _PLASMA_CORE_DBLAS_H_
#define _PLASMA_CORE_DBLAS_H_
#define REAL

#ifdef __cplusplus
extern "C" {
#endif

/** ****************************************************************************
 *  Declarations of serial kernels - alphabetical order
 **/

int CORE_dlarfx2(int side, int N,
                 double V,
                 double TAU,
                 double *C1, int LDC1,
                 double *C2, int LDC2);
int CORE_dlarfx2c(int uplo,
                  double V,
                  double TAU,
                  double *C1,
                  double *C2,
                  double *C3);
int CORE_dlarfx2ce(int uplo,
                double *V,
                double *TAU,
                double *C1,
                double *C2,
                double *C3);
int CORE_dhbelr(int uplo, int N,
                PLASMA_desc *A,
                double *V,
                double *TAU,
                int st,
                int ed,
                int eltsize);
int CORE_dhbrce(int uplo, int N,
                PLASMA_desc *A,
                double *V,
                double *TAU,
                int st,
                int ed,
                int eltsize);
int CORE_dhblrx(int uplo, int N,
                PLASMA_desc *A,
                double *V,
                double *TAU,
                int st,
                int ed,
                int eltsize);
int CORE_dgbelr(int uplo, int N,
                PLASMA_desc *A,
                double *V,
                double *TAU,
                int st,
                int ed,
                int eltsize);
int CORE_dgbrce(int uplo, int N,
                PLASMA_desc *A,
                double *V,
                double *TAU,
                int st,
                int ed,
                int eltsize);
int CORE_dgblrx(int uplo, int N,
                PLASMA_desc *A,
                double *V,
                double *TAU,
                int st,
                int ed,
                int eltsize);
void CORE_dasum(int storev, int uplo, int M, int N,
                 double *A, int lda, double *work);
void CORE_daxpy(int M, int N,  double alpha,
                double *A, int LDA,
                double *B, int LDB);
void CORE_dbrdalg(PLASMA_enum uplo, int N, int NB, 
                  PLASMA_desc *pA, double *C, double *S,
                  int i, int j, int m, int grsiz);
int  CORE_dgelqt(int M, int N, int IB,
                 double *A, int LDA,
                 double *T, int LDT,
                 double *TAU, double *WORK);
void CORE_dgemm(int transA, int transB,
                int M, int N, int K,
                double alpha, double *A, int LDA,
                double *B, int LDB,
                double beta, double *C, int LDC);
int  CORE_dgeqrt(int M, int N, int IB,
                 double *A, int LDA,
                 double *T, int LDT,
                 double *TAU, double *WORK);
int  CORE_dgessm(int M, int N, int K, int IB,
                 int *IPIV,
                 double *L, int LDL,
                 double *A, int LDA);
int  CORE_dgetrf(int M, int N, 
                 double *A, int LDA,
                 int *IPIV, int *INFO);
int  CORE_dgetrf_incpiv(int M, int N, int IB,
                        double *A, int LDA,
                        int *IPIV, int *INFO);
int  CORE_dgetrf_reclap(const int M, const int N,
                        double *A, const int LDA,
                        int *IPIV, int *info);
int  CORE_dgetrf_rectil(const PLASMA_desc A, int *IPIV, int *info);
void CORE_dgetrip(int m, int n, double *A, 
                  double *work);
#ifdef COMPLEX
void CORE_dsygst(int itype, int uplo, int N,
                 double *A, int LDA,
                 double *B, int LDB, int *INFO);
void CORE_dsymm(int side, int uplo,
                int M, int N,
                double alpha, double *A, int LDA,
                double *B, int LDB,
                double beta, double *C, int LDC);
void CORE_dsyrk(int uplo, int trans,
                int N, int K,
                double alpha, double *A, int LDA,
                double beta, double *C, int LDC);
void CORE_dsyr2k(int uplo, int trans,
                 int N, int K,
                 double alpha, double *A, int LDA,
                 double *B, int LDB,
                 double beta, double *C, int LDC);
int  CORE_dsyrfb(PLASMA_enum uplo, int N, int K, int IB, int NB,
                 double *A, int LDA,
                 double *T, int LDT,
                 double *C, int LDC,
                 double *WORK, int LDWORK);
#endif
void CORE_dlacpy(PLASMA_enum uplo, int M, int N,
                 double *A, int LDA,
                 double *B, int LDB);
void CORE_dlange(int norm, int M, int N,
                 double *A, int LDA,
                 double *work, double *normA);
#ifdef COMPLEX
void CORE_dlansy(int norm, int uplo, int N,
                 double *A, int LDA,
                 double *work, double *normA);
#endif
void CORE_dlansy(int norm, int uplo, int N,
                 double *A, int LDA,
                 double *work, double *normA);
void CORE_dlaset(PLASMA_enum uplo, int n1, int n2, double alpha,
                 double beta, double *tileA, int ldtilea);
void CORE_dlaset2(PLASMA_enum uplo, int n1, int n2, double alpha,
                  double *tileA, int ldtilea);
void CORE_dlaswp(int N, double *A, int LDA, 
                 int I1,  int I2, int *IPIV, int INC);
int  CORE_dlaswp_ontile(PLASMA_desc descA, int i1, int i2, int *ipiv, int inc);
void CORE_dlauum(int uplo, int N, double *A, int LDA);
void CORE_dplgsy(double bump, int m, int n, double *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
void CORE_dplgsy(double bump, int m, int n, double *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
void CORE_dplrnt(int m, int n, double *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
void CORE_dpotrf(int uplo, int N, double *A, int LDA, int *INFO);
void CORE_dshift(int s, int m, int n, int L,
                 double *A);
void CORE_dshiftw(int s, int cl, int m, int n, int L,
                  double *A, double *W);
int  CORE_dssssm(int M1, int N1, int M2, int N2, int K, int IB,
                 double *A1, int LDA1,
                 double *A2, int LDA2,
                 double *L1, int LDL1,
                 double *L2, int LDL2,
                 int *IPIV);
void CORE_dsymm(int side, int uplo,
                int M, int N,
                double alpha, double *A, int LDA,
                double *B, int LDB,
                double beta, double *C, int LDC);
void CORE_dsyrk(int uplo, int trans,
                int N, int K,
                double alpha, double *A, int LDA,
                double beta, double *C, int LDC);
void CORE_dsyr2k(int uplo, int trans,
                 int N, int K,
                 double alpha, double *A, int LDA,
                 double *B, int LDB,
                 double beta, double *C, int LDC);
void CORE_dswpab(int i, int n1, int n2,
                 double *A, double *work);
int  CORE_dswptr_ontile(PLASMA_desc descA, int i1, int i2, int *ipiv, int inc,
                        double *Akk, int ldak);
void CORE_dtrdalg(PLASMA_enum uplo, int N, int NB, 
                  PLASMA_desc *pA, double *C, double *S,
                  int i, int j, int m, int grsiz);
void CORE_dtrmm(int side, int uplo,
                int transA, int diag,
                int M, int N,
                double alpha, double *A, int LDA,
                double *B, int LDB);
void CORE_dtrsm(int side, int uplo,
                int transA, int diag,
                int M, int N,
                double alpha, double *A, int LDA,
                double *B, int LDB);
void CORE_dtrtri(int uplo, int diag, int N, double *A, int LDA, int *info);
int  CORE_dtslqt(int M, int N, int IB,
                 double *A1, int LDA1,
                 double *A2, int LDA2,
                 double *T, int LDT,
                 double *TAU, double *WORK);
int  CORE_dtsmlq(int side, int trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 double *A1, int LDA1,
                 double *A2, int LDA2,
                 double *V, int LDV,
                 double *T, int LDT,
                 double *WORK, int LDWORK);
int CORE_dtsmlq_corner( int m1, int n1, int m2, int n2, int m3, int n3,
                        int k, int ib, int nb,
                        double *A1, int lda1,
                        double *A2, int lda2,
                        double *A3, int lda3,
                        double *V, int ldv,
                        double *T, int ldt,
                        double *WORK, int ldwork);
int CORE_dtsmlq_sytra1( int side, int trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        double *A1, int lda1,
                        double *A2, int lda2,
                        double *V, int ldv,
                        double *T, int ldt,
                        double *WORK, int ldwork);
int  CORE_dtsmqr(int side, int trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 double *A1, int LDA1,
                 double *A2, int LDA2,
                 double *V, int LDV,
                 double *T, int LDT,
                 double *WORK, int LDWORK);
int CORE_dtsmqr_corner( int m1, int n1, int m2, int n2, int m3, int n3,
                        int k, int ib, int nb,
                        double *A1, int lda1,
                        double *A2, int lda2,
                        double *A3, int lda3,
                        double *V, int ldv,
                        double *T, int ldt,
                        double *WORK, int ldwork);
int CORE_dtsmqr_sytra1( int side, int trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        double *A1, int lda1,
                        double *A2, int lda2,
                        double *V, int ldv,
                        double *T, int ldt,
                        double *WORK, int ldwork);
int  CORE_dtsqrt(int M, int N, int IB,
                 double *A1, int LDA1,
                 double *A2, int LDA2,
                 double *T, int LDT,
                 double *TAU, double *WORK);
int  CORE_dtsrfb(int side, int trans, int direct, int storev,
                 int M1, int N1, int M2, int N2, int IB,
                 double *A1, int LDA1,
                 double *A2, int LDA2,
                 double *V, int LDV,
                 double *T, int LDT,
                 double *WORK, int LDWORK);
int  CORE_dtstrf(int M, int N, int IB, int NB,
                 double *U, int LDU,
                 double *A, int LDA,
                 double *L, int LDL,
                 int *IPIV, double *WORK,
                 int LDWORK, int *INFO);
int  CORE_dttmqr(int side, int trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 double *A1, int LDA1,
                 double *A2, int LDA2,
                 double *V, int LDV,
                 double *T, int LDT,
                 double *WORK, int LDWORK);
int  CORE_dttqrt(int M, int N, int IB,
                 double *A1, int LDA1,
                 double *A2, int LDA2,
                 double *T, int LDT,
                 double *TAU,
                 double *WORK);
int  CORE_dttmlq(int side, int trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 double *A1, int LDA1,
                 double *A2, int LDA2,
                 double *V, int LDV,
                 double *T, int LDT,
                 double *WORK, int LDWORK);
int  CORE_dttlqt(int M, int N, int IB,
                 double *A1, int LDA1,
                 double *A2, int LDA2,
                 double *T, int LDT,
                 double *TAU,
                 double *WORK);
int  CORE_dttrfb(int side, int trans, int direct, int storev,
                 int M1, int N1, int M2, int N2, int IB,
                 double *A1, int LDA1,
                 double *A2, int LDA2,
                 double *V, int LDV,
                 double *T, int LDT,
                 double *WORK, int LDWORK);
int  CORE_dormlq(int side, int trans,
                 int M, int N, int IB, int K,
                 double *V, int LDV,
                 double *T, int LDT,
                 double *C, int LDC,
                 double *WORK, int LDWORK);
int  CORE_dormqr(int side, int trans,
                 int M, int N, int K, int IB,
                 double *V, int LDV,
                 double *T, int LDT,
                 double *C, int LDC,
                 double *WORK, int LDWORK);

/** ****************************************************************************
 *  Declarations of QUARK wrappers (called by PLASMA) - alphabetical order
 **/
void QUARK_CORE_dasum(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum storev, PLASMA_enum uplo, int m, int n,
                       double *A, int lda, int szeA,
                       double *work, int szeW);
void QUARK_CORE_dasum_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum storev, PLASMA_enum uplo, int m, int n,
                          double *A, int lda, int szeA,
                          double *work, int szeW,
                          double *fake, int szeF);
void QUARK_CORE_daxpy(Quark *quark, Quark_Task_Flags *task_flags,
                      int m, int n, int nb, double alpha,
                      double *A, int lda,
                      double *B, int ldb);
void QUARK_CORE_dbrdalg(Quark *quark, Quark_Task_Flags *task_flags,
                        int uplo,
                        int N, int NB,
                        PLASMA_desc *A,
                        double *C,
                        double *S,
                        int i, int j, int m, int grsiz, int BAND,
                        int *PCOL, int *ACOL, int *MCOL);
void QUARK_CORE_dgelqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       double *A, int lda,
                       double *T, int ldt);
void QUARK_CORE_dgemm(Quark *quark, Quark_Task_Flags *task_flags,
                      int transA, int transB,
                      int m, int n, int k, int nb,
                      double alpha, double *A, int lda,
                      double *B, int ldb,
                      double beta, double *C, int ldc);
void QUARK_CORE_dgemm2( Quark *quark, Quark_Task_Flags *task_flags,
                        int transA, int transB,
                        int m, int n, int k, int nb,
                        double alpha, double *A, int lda,
                        double *B, int ldb,
                        double beta, double *C, int ldc);
void QUARK_CORE_dgemm_f2(Quark *quark, Quark_Task_Flags *task_flags,
                         int transA, int transB,
                         int m, int n, int k, int nb,
                         double alpha, double *A, int lda,
                         double *B, int ldb,
                         double beta, double *C, int ldc,
                         double *fake1, int szefake1, int flag1,
                         double *fake2, int szefake2, int flag2);
void QUARK_CORE_dgemm_p2(Quark *quark, Quark_Task_Flags *task_flags,
                         int transA, int transB,
                         int m, int n, int k, int nb,
                         double alpha, double *A, int lda,
                         double **B, int ldb,
                         double beta, double *C, int ldc);
void QUARK_CORE_dgemm_p2f1(Quark *quark, Quark_Task_Flags *task_flags,
                           int transA, int transB,
                           int m, int n, int k, int nb,
                           double alpha, double *A, int lda,
                           double **B, int ldb,
                           double beta, double *C, int ldc,
                           double *fake1, int szefake1, int flag1);
void QUARK_CORE_dgemm_p3(Quark *quark, Quark_Task_Flags *task_flags,
                         int transA, int transB,
                         int m, int n, int k, int nb,
                         double alpha, double *A, int lda,
                         double *B, int ldb,
                         double beta, double **C, int ldc);
void QUARK_CORE_dgeqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       double *A, int lda,
                       double *T, int ldt);
void QUARK_CORE_dgessm(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int k, int ib, int nb,
                       int *IPIV,
                       double *L, int ldl,
                       double *A, int lda);
void QUARK_CORE_dgetrf(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int nb,
                       double *A, int lda,
                       int *IPIV,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       PLASMA_bool check_info, int iinfo);
void QUARK_CORE_dgetrf_incpiv(Quark *quark, Quark_Task_Flags *task_flags,
                              int m, int n, int ib, int nb,
                              double *A, int lda,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo);
void QUARK_CORE_dgetrf_reclap(Quark *quark, Quark_Task_Flags *task_flags,
                              int m, int n, int nb,
                              double *A, int lda,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo,
                              int nbthread);
void QUARK_CORE_dgetrf_rectil(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc A, double *Amn, int size,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo,
                              int nbthread);
void QUARK_CORE_dgetrip(Quark *quark, Quark_Task_Flags *task_flags,
                        int m, int n, double *A, int szeA);
void QUARK_CORE_dgetrip_f1(Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n, double *A, int szeA,
                           double *fake, int szeF, int paramF);
void QUARK_CORE_dgetrip_f2(Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n, double *A, int szeA,
                           double *fake1, int szeF1, int paramF1,
                           double *fake2, int szeF2, int paramF2);
void QUARK_CORE_dsymm(Quark *quark, Quark_Task_Flags *task_flags,
                      int side, int uplo,
                      int m, int n, int nb,
                      double alpha, double *A, int lda,
                      double *B, int ldb,
                      double beta, double *C, int ldc);
void QUARK_CORE_dsygst(Quark *quark, Quark_Task_Flags *task_flags,
                       int itype, int uplo, int N,
                       double *A, int LDA,
                       double *B, int LDB,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo);
void QUARK_CORE_dsyrk(Quark *quark, Quark_Task_Flags *task_flags,
                      int uplo, int trans,
                      int n, int k, int nb,
                      double alpha, double *A, int lda,
                      double beta, double *C, int ldc);
void QUARK_CORE_dsyr2k(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int trans,
                       int n, int k, int nb,
                       double alpha, double *A, int lda,
                       double *B, int LDB,
                       double beta, double *C, int ldc);
void QUARK_CORE_dsyrfb(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo,
                       int n, int k, int ib, int nb,
                       double *A, int lda,
                       double *T, int ldt,
                       double *C, int ldc);
void QUARK_CORE_dlacpy(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int m, int n, int mb,
                       double *A, int lda,
                       double *B, int ldb);
void QUARK_CORE_dlange(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, int M, int N,
                       double *A, int LDA, int szeA,
                       int szeW, double *result);
void QUARK_CORE_dlange_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          int norm, int M, int N,
                          double *A, int LDA, int szeA,
                          int szeW, double *result,
                          double *fake, int szeF);
#ifdef COMPLEX
void QUARK_CORE_dlansy(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, int uplo, int N,
                       double *A, int LDA, int szeA,
                       int szeW, double *result);
void QUARK_CORE_dlansy_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          int norm, int uplo, int N,
                          double *A, int LDA, int szeA,
                          int szeW, double *result,
                          double *fake, int szeF);
#endif
void QUARK_CORE_dlansy(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, int uplo, int N,
                       double *A, int LDA, int szeA,
                       int szeW, double *result);
void QUARK_CORE_dlansy_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          int norm, int uplo, int N,
                          double *A, int LDA, int szeA,
                          int szeW, double *result,
                          double *fake, int szeF);
void QUARK_CORE_dlaset(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int n1, int n2, double alpha,
                       double beta, double *tileA, int ldtilea);
void QUARK_CORE_dlaset2(Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum uplo, int n1, int n2, double alpha,
                        double *tileA, int ldtilea);
void QUARK_CORE_dlaswp(Quark *quark, Quark_Task_Flags *task_flags,
                       int n, double *A, int lda, 
                       int i1,  int i2, int *ipiv, int inc);
void QUARK_CORE_dlaswp_f2(Quark *quark, Quark_Task_Flags *task_flags,
                          int n, double *A, int lda, 
                          int i1,  int i2, int *ipiv, int inc,
                          double *fake1, int szefake1, int flag1,
                          double *fake2, int szefake2, int flag2);
void QUARK_CORE_dlaswp_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, double *A, 
                              int i1,  int i2, int *ipiv, int inc, double *fakepanel);
void QUARK_CORE_dlaswp_ontile_f2(Quark *quark, Quark_Task_Flags *task_flags,
                                 PLASMA_desc descA, double *A, 
                                 int i1,  int i2, int *ipiv, int inc, 
                                 double *fake1, int szefake1, int flag1,
                                 double *fake2, int szefake2, int flag2);
void QUARK_CORE_dlauum(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int n, int nb,
                       double *A, int lda);
void QUARK_CORE_dplgsy(Quark *quark, Quark_Task_Flags *task_flags,
                       double bump, int m, int n, double *A, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void QUARK_CORE_dplgsy(Quark *quark, Quark_Task_Flags *task_flags,
                       double bump, int m, int n, double *A, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void QUARK_CORE_dplrnt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, double *A, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void QUARK_CORE_dpotrf(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int n, int nb,
                       double *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo);
void QUARK_CORE_dshift( Quark *quark, Quark_Task_Flags *task_flags,
                        int s, int m, int n, int L,
                        double *A);
void QUARK_CORE_dshiftw(Quark *quark, Quark_Task_Flags *task_flags,
                        int s, int cl, int m, int n, int L,
                        double *A, double *W);
void QUARK_CORE_dssssm(Quark *quark, Quark_Task_Flags *task_flags,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       double *L1, int ldl1,
                       double *L2, int ldl2,
                       int *IPIV);
void QUARK_CORE_dsymm(Quark *quark, Quark_Task_Flags *task_flags,
                      int side, int uplo,
                      int m, int n, int nb,
                      double alpha, double *A, int lda,
                      double *B, int ldb,
                      double beta, double *C, int ldc);
void QUARK_CORE_dsyrk(Quark *quark, Quark_Task_Flags *task_flags,
                      int uplo, int trans,
                      int n, int k, int nb,
                      double alpha, double *A, int lda,
                      double beta, double *C, int ldc);
void QUARK_CORE_dsyr2k(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int trans,
                       int n, int k, int nb,
                       double alpha, double *A, int lda,
                       double *B, int LDB,
                       double beta, double *C, int ldc);
void QUARK_CORE_dswpab(Quark *quark, Quark_Task_Flags *task_flags,
                       int i, int n1, int n2,
                       double *A, int szeA);
void QUARK_CORE_dswptr_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, double *Aij, 
                              int i1,  int i2, int *ipiv, int inc, 
                              double *Akk, int ldak);
void QUARK_CORE_dtrdalg(Quark *quark, Quark_Task_Flags *task_flags,
                        int uplo,
                        int N, int NB,
                        PLASMA_desc *A,
                        double *C,
                        double *S,
                        int i, int j, int m, int grsiz, int BAND,
                        int *PCOL, int *ACOL, int *MCOL);
void QUARK_CORE_dtrmm(Quark *quark, Quark_Task_Flags *task_flags,
                      int side, int uplo, int transA, int diag,
                      int m, int n, int nb,
                      double alpha, double *A, int lda,
                      double *B, int ldb);
void QUARK_CORE_dtrmm_p2(Quark *quark, Quark_Task_Flags *task_flags,
                         int side, int uplo, int transA, int diag,
                         int m, int n, int nb,
                         double alpha, double *A, int lda,
                         double **B, int ldb);
void QUARK_CORE_dtrsm(Quark *quark, Quark_Task_Flags *task_flags,
                      int side, int uplo, int transA, int diag,
                      int m, int n, int nb,
                      double alpha, double *A, int lda,
                      double *B, int ldb);
void QUARK_CORE_dtrtri(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int diag, int n, int nb,
                       double *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo);
void QUARK_CORE_dtslqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       double *T, int ldt);
void QUARK_CORE_dtsmlq(Quark *quark, Quark_Task_Flags *task_flags,
                       int side, int trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       double *V, int ldv,
                       double *T, int ldt);
void QUARK_CORE_dtsmlq_sytra1(Quark *quark, Quark_Task_Flags *task_flags,
                              int side, int trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              double *A1, int lda1,
                              double *A2, int lda2,
                              double *V, int ldv,
                              double *T, int ldt);
void QUARK_CORE_dtsmlq_corner(Quark *quark, Quark_Task_Flags *task_flags,
                              int m1, int n1, int m2, int n2, int m3, int n3, int k, int ib, int nb,
                              double *A1, int lda1,
                              double *A2, int lda2,
                              double *A3, int lda3,
                              double *V, int ldv,
                              double *T, int ldt);
void QUARK_CORE_dtsmqr(Quark *quark, Quark_Task_Flags *task_flags,
                       int side, int trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       double *V, int ldv,
                       double *T, int ldt);
void QUARK_CORE_dtsmqr_sytra1(Quark *quark, Quark_Task_Flags *task_flags,
                              int side, int trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              double *A1, int lda1,
                              double *A2, int lda2,
                              double *V, int ldv,
                              double *T, int ldt);
void QUARK_CORE_dtsmqr_corner(Quark *quark, Quark_Task_Flags *task_flags,
                              int m1, int n1, int m2, int n2, int m3, int n3, int k, int ib, int nb,
                              double *A1, int lda1,
                              double *A2, int lda2,
                              double *A3, int lda3,
                              double *V, int ldv,
                              double *T, int ldt);
void QUARK_CORE_dtsqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       double *T, int ldt);
void QUARK_CORE_dtstrf(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       double *U, int ldu,
                       double *A, int lda,
                       double *L, int ldl,
                       int *IPIV,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       PLASMA_bool check_info, int iinfo);
void QUARK_CORE_dttmqr(Quark *quark, Quark_Task_Flags *task_flags,
                       int side, int trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       double *V, int ldv,
                       double *T, int ldt);
void QUARK_CORE_dttqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       double *T, int ldt);
void QUARK_CORE_dttmlq(Quark *quark, Quark_Task_Flags *task_flags,
                       int side, int trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       double *V, int ldv,
                       double *T, int ldt);
void QUARK_CORE_dttlqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       double *T, int ldt);
void QUARK_CORE_dormlq(Quark *quark, Quark_Task_Flags *task_flags,
                       int side, int trans,
                       int m, int n, int ib,  int nb, int k,
                       double *A, int lda,
                       double *T, int ldt,
                       double *C, int ldc);
void QUARK_CORE_dormqr(Quark *quark, Quark_Task_Flags *task_flags,
                       int side, int trans,
                       int m, int n, int k, int ib, int nb,
                       double *A, int lda,
                       double *T, int ldt,
                       double *C, int ldc);

/** ****************************************************************************
 *  Declarations of QUARK wrappers (called by QUARK) - alphabetical order
 **/
void CORE_dasum_quark(Quark *quark);
void CORE_dasum_f1_quark(Quark *quark);
void CORE_daxpy_quark(Quark *quark);
void CORE_dbrdalg_quark(Quark *quark);
void CORE_dgelqt_quark(Quark *quark);
void CORE_dgemm_quark(Quark *quark);
void CORE_dgeqrt_quark(Quark *quark);
void CORE_dgessm_quark(Quark *quark);
void CORE_dgetrf_quark(Quark *quark);
void CORE_dgetrf_incpiv_quark(Quark *quark);
void CORE_dgetrf_reclap_quark(Quark *quark);
void CORE_dgetrf_rectil_quark(Quark* quark);
void CORE_dgetrip_quark(Quark *quark);
void CORE_dgetrip_f1_quark(Quark *quark);
void CORE_dgetrip_f2_quark(Quark *quark);
#ifdef COMPLEX
void CORE_dsymm_quark(Quark *quark);
void CORE_dsyrk_quark(Quark *quark);
void CORE_dsyr2k_quark(Quark *quark);
#endif
void CORE_dsygst_quark(Quark *quark);
void CORE_dsyrfb_quark(Quark *quark);
void CORE_dlacpy_quark(Quark *quark);
void CORE_dlange_quark(Quark *quark);
void CORE_dlange_f1_quark(Quark *quark);
#ifdef COMPLEX
void CORE_dlansy_quark(Quark *quark);
void CORE_dlansy_f1_quark(Quark *quark);
#endif
void CORE_dlansy_quark(Quark *quark);
void CORE_dlansy_f1_quark(Quark *quark);
void CORE_dlaset_quark(Quark *quark);
void CORE_dlaset2_quark(Quark *quark);
void CORE_dlauum_quark(Quark *quark);
void CORE_dplgsy_quark(Quark *quark);
void CORE_dplgsy_quark(Quark *quark);
void CORE_dplrnt_quark(Quark *quark);
void CORE_dpotrf_quark(Quark *quark);
void CORE_dshift_quark(Quark *quark);
void CORE_dshiftw_quark(Quark *quark);
void CORE_dssssm_quark(Quark *quark);
void CORE_dsymm_quark(Quark *quark);
void CORE_dsyrk_quark(Quark *quark);
void CORE_dsyr2k_quark(Quark *quark);
void CORE_dswpab_quark(Quark *quark);
void CORE_dswptr_ontile_quark(Quark *quark);
void CORE_dtrdalg_quark(Quark *quark);
void CORE_dtrmm_quark(Quark *quark);
void CORE_dtrsm_quark(Quark *quark);
void CORE_dtrtri_quark(Quark *quark);
void CORE_dtslqt_quark(Quark *quark);
void CORE_dtsmlq_quark(Quark *quark);
void CORE_dtsmlq_sytra1_quark(Quark *quark);
void CORE_dtsmlq_corner_quark(Quark *quark);
void CORE_dtsmqr_quark(Quark *quark);
void CORE_dtsmqr_sytra1_quark(Quark *quark);
void CORE_dtsmqr_corner_quark(Quark *quark);
void CORE_dtsqrt_quark(Quark *quark);
void CORE_dtstrf_quark(Quark *quark);
void CORE_dttmqr_quark(Quark *quark);
void CORE_dttqrt_quark(Quark *quark);
void CORE_dttmlq_quark(Quark *quark);
void CORE_dttlqt_quark(Quark *quark);
void CORE_dormlq_quark(Quark *quark);
void CORE_dormqr_quark(Quark *quark);

void CORE_dlaswp_quark(Quark* quark);
void CORE_dlaswp_f2_quark(Quark* quark);
void CORE_dlaswp_ontile_quark(Quark *quark);
void CORE_dlaswp_ontile_f2_quark(Quark *quark);
void CORE_dtrmm_p2_quark(Quark* quark);
void CORE_dgemm_f2_quark(Quark* quark);
void CORE_dgemm_p2_quark(Quark* quark);
void CORE_dgemm_p2f1_quark(Quark* quark);
void CORE_dgemm_p3_quark(Quark* quark);





void CORE_dtrdalg_v2_quark(Quark *quark);


void QUARK_CORE_dtrdalg_v2(Quark *quark, Quark_Task_Flags *task_flags,
                        int uplo,
                        PLASMA_desc *A,
                        double *C,
                        double *S,
                        int grsiz, int s, int id, int blktile);

void CORE_dtrdalg_v2(PLASMA_enum uplo,  
                  PLASMA_desc *pA, double *C, double *S,
                  int grsiz, int s, int id, int blktile);


#ifdef __cplusplus
}
#endif

#endif
