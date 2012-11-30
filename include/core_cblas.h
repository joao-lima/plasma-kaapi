/**
 *
 * @file core_cblas.h
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
 * @generated c Thu Sep 15 12:08:48 2011
 *
 **/
#ifndef _PLASMA_CORE_CBLAS_H_
#define _PLASMA_CORE_CBLAS_H_
#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

/** ****************************************************************************
 *  Declarations of serial kernels - alphabetical order
 **/

int CORE_clarfx2(int side, int N,
                 PLASMA_Complex32_t V,
                 PLASMA_Complex32_t TAU,
                 PLASMA_Complex32_t *C1, int LDC1,
                 PLASMA_Complex32_t *C2, int LDC2);
int CORE_clarfx2c(int uplo,
                  PLASMA_Complex32_t V,
                  PLASMA_Complex32_t TAU,
                  PLASMA_Complex32_t *C1,
                  PLASMA_Complex32_t *C2,
                  PLASMA_Complex32_t *C3);
int CORE_clarfx2ce(int uplo,
                PLASMA_Complex32_t *V,
                PLASMA_Complex32_t *TAU,
                PLASMA_Complex32_t *C1,
                PLASMA_Complex32_t *C2,
                PLASMA_Complex32_t *C3);
int CORE_chbelr(int uplo, int N,
                PLASMA_desc *A,
                PLASMA_Complex32_t *V,
                PLASMA_Complex32_t *TAU,
                int st,
                int ed,
                int eltsize);
int CORE_chbrce(int uplo, int N,
                PLASMA_desc *A,
                PLASMA_Complex32_t *V,
                PLASMA_Complex32_t *TAU,
                int st,
                int ed,
                int eltsize);
int CORE_chblrx(int uplo, int N,
                PLASMA_desc *A,
                PLASMA_Complex32_t *V,
                PLASMA_Complex32_t *TAU,
                int st,
                int ed,
                int eltsize);
int CORE_cgbelr(int uplo, int N,
                PLASMA_desc *A,
                PLASMA_Complex32_t *V,
                PLASMA_Complex32_t *TAU,
                int st,
                int ed,
                int eltsize);
int CORE_cgbrce(int uplo, int N,
                PLASMA_desc *A,
                PLASMA_Complex32_t *V,
                PLASMA_Complex32_t *TAU,
                int st,
                int ed,
                int eltsize);
int CORE_cgblrx(int uplo, int N,
                PLASMA_desc *A,
                PLASMA_Complex32_t *V,
                PLASMA_Complex32_t *TAU,
                int st,
                int ed,
                int eltsize);
void CORE_scasum(int storev, int uplo, int M, int N,
                 PLASMA_Complex32_t *A, int lda, float *work);
void CORE_caxpy(int M, int N,  PLASMA_Complex32_t alpha,
                PLASMA_Complex32_t *A, int LDA,
                PLASMA_Complex32_t *B, int LDB);
void CORE_cbrdalg(PLASMA_enum uplo, int N, int NB, 
                  PLASMA_desc *pA, PLASMA_Complex32_t *C, PLASMA_Complex32_t *S,
                  int i, int j, int m, int grsiz);
int  CORE_cgelqt(int M, int N, int IB,
                 PLASMA_Complex32_t *A, int LDA,
                 PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *TAU, PLASMA_Complex32_t *WORK);
void CORE_cgemm(int transA, int transB,
                int M, int N, int K,
                PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int LDA,
                PLASMA_Complex32_t *B, int LDB,
                PLASMA_Complex32_t beta, PLASMA_Complex32_t *C, int LDC);
int  CORE_cgeqrt(int M, int N, int IB,
                 PLASMA_Complex32_t *A, int LDA,
                 PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *TAU, PLASMA_Complex32_t *WORK);
int  CORE_cgessm(int M, int N, int K, int IB,
                 int *IPIV,
                 PLASMA_Complex32_t *L, int LDL,
                 PLASMA_Complex32_t *A, int LDA);
int  CORE_cgetrf(int M, int N, 
                 PLASMA_Complex32_t *A, int LDA,
                 int *IPIV, int *INFO);
int  CORE_cgetrf_incpiv(int M, int N, int IB,
                        PLASMA_Complex32_t *A, int LDA,
                        int *IPIV, int *INFO);
int  CORE_cgetrf_reclap(const int M, const int N,
                        PLASMA_Complex32_t *A, const int LDA,
                        int *IPIV, int *info);
int  CORE_cgetrf_rectil(const PLASMA_desc A, int *IPIV, int *info);
void CORE_cgetrip(int m, int n, PLASMA_Complex32_t *A, 
                  PLASMA_Complex32_t *work);
#ifdef COMPLEX
void CORE_chegst(int itype, int uplo, int N,
                 PLASMA_Complex32_t *A, int LDA,
                 PLASMA_Complex32_t *B, int LDB, int *INFO);
void CORE_chemm(int side, int uplo,
                int M, int N,
                PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int LDA,
                PLASMA_Complex32_t *B, int LDB,
                PLASMA_Complex32_t beta, PLASMA_Complex32_t *C, int LDC);
void CORE_cherk(int uplo, int trans,
                int N, int K,
                float alpha, PLASMA_Complex32_t *A, int LDA,
                float beta, PLASMA_Complex32_t *C, int LDC);
void CORE_cher2k(int uplo, int trans,
                 int N, int K,
                 PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int LDA,
                 PLASMA_Complex32_t *B, int LDB,
                 float beta, PLASMA_Complex32_t *C, int LDC);
int  CORE_cherfb(PLASMA_enum uplo, int N, int K, int IB, int NB,
                 PLASMA_Complex32_t *A, int LDA,
                 PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *C, int LDC,
                 PLASMA_Complex32_t *WORK, int LDWORK);
#endif
void CORE_clacpy(PLASMA_enum uplo, int M, int N,
                 PLASMA_Complex32_t *A, int LDA,
                 PLASMA_Complex32_t *B, int LDB);
void CORE_clange(int norm, int M, int N,
                 PLASMA_Complex32_t *A, int LDA,
                 float *work, float *normA);
#ifdef COMPLEX
void CORE_clanhe(int norm, int uplo, int N,
                 PLASMA_Complex32_t *A, int LDA,
                 float *work, float *normA);
#endif
void CORE_clansy(int norm, int uplo, int N,
                 PLASMA_Complex32_t *A, int LDA,
                 float *work, float *normA);
void CORE_claset(PLASMA_enum uplo, int n1, int n2, PLASMA_Complex32_t alpha,
                 PLASMA_Complex32_t beta, PLASMA_Complex32_t *tileA, int ldtilea);
void CORE_claset2(PLASMA_enum uplo, int n1, int n2, PLASMA_Complex32_t alpha,
                  PLASMA_Complex32_t *tileA, int ldtilea);
void CORE_claswp(int N, PLASMA_Complex32_t *A, int LDA, 
                 int I1,  int I2, int *IPIV, int INC);
int  CORE_claswp_ontile(PLASMA_desc descA, int i1, int i2, int *ipiv, int inc);
void CORE_clauum(int uplo, int N, PLASMA_Complex32_t *A, int LDA);
void CORE_cplghe(float bump, int m, int n, PLASMA_Complex32_t *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
void CORE_cplgsy(PLASMA_Complex32_t bump, int m, int n, PLASMA_Complex32_t *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
void CORE_cplrnt(int m, int n, PLASMA_Complex32_t *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
void CORE_cpotrf(int uplo, int N, PLASMA_Complex32_t *A, int LDA, int *INFO);
void CORE_cshift(int s, int m, int n, int L,
                 PLASMA_Complex32_t *A);
void CORE_cshiftw(int s, int cl, int m, int n, int L,
                  PLASMA_Complex32_t *A, PLASMA_Complex32_t *W);
int  CORE_cssssm(int M1, int N1, int M2, int N2, int K, int IB,
                 PLASMA_Complex32_t *A1, int LDA1,
                 PLASMA_Complex32_t *A2, int LDA2,
                 PLASMA_Complex32_t *L1, int LDL1,
                 PLASMA_Complex32_t *L2, int LDL2,
                 int *IPIV);
void CORE_csymm(int side, int uplo,
                int M, int N,
                PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int LDA,
                PLASMA_Complex32_t *B, int LDB,
                PLASMA_Complex32_t beta, PLASMA_Complex32_t *C, int LDC);
void CORE_csyrk(int uplo, int trans,
                int N, int K,
                PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int LDA,
                PLASMA_Complex32_t beta, PLASMA_Complex32_t *C, int LDC);
void CORE_csyr2k(int uplo, int trans,
                 int N, int K,
                 PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int LDA,
                 PLASMA_Complex32_t *B, int LDB,
                 PLASMA_Complex32_t beta, PLASMA_Complex32_t *C, int LDC);
void CORE_cswpab(int i, int n1, int n2,
                 PLASMA_Complex32_t *A, PLASMA_Complex32_t *work);
int  CORE_cswptr_ontile(PLASMA_desc descA, int i1, int i2, int *ipiv, int inc,
                        PLASMA_Complex32_t *Akk, int ldak);
void CORE_ctrdalg(PLASMA_enum uplo, int N, int NB, 
                  PLASMA_desc *pA, PLASMA_Complex32_t *C, PLASMA_Complex32_t *S,
                  int i, int j, int m, int grsiz);
void CORE_ctrmm(int side, int uplo,
                int transA, int diag,
                int M, int N,
                PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int LDA,
                PLASMA_Complex32_t *B, int LDB);
void CORE_ctrsm(int side, int uplo,
                int transA, int diag,
                int M, int N,
                PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int LDA,
                PLASMA_Complex32_t *B, int LDB);
void CORE_ctrtri(int uplo, int diag, int N, PLASMA_Complex32_t *A, int LDA, int *info);
int  CORE_ctslqt(int M, int N, int IB,
                 PLASMA_Complex32_t *A1, int LDA1,
                 PLASMA_Complex32_t *A2, int LDA2,
                 PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *TAU, PLASMA_Complex32_t *WORK);
int  CORE_ctsmlq(int side, int trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 PLASMA_Complex32_t *A1, int LDA1,
                 PLASMA_Complex32_t *A2, int LDA2,
                 PLASMA_Complex32_t *V, int LDV,
                 PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *WORK, int LDWORK);
int CORE_ctsmlq_corner( int m1, int n1, int m2, int n2, int m3, int n3,
                        int k, int ib, int nb,
                        PLASMA_Complex32_t *A1, int lda1,
                        PLASMA_Complex32_t *A2, int lda2,
                        PLASMA_Complex32_t *A3, int lda3,
                        PLASMA_Complex32_t *V, int ldv,
                        PLASMA_Complex32_t *T, int ldt,
                        PLASMA_Complex32_t *WORK, int ldwork);
int CORE_ctsmlq_hetra1( int side, int trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        PLASMA_Complex32_t *A1, int lda1,
                        PLASMA_Complex32_t *A2, int lda2,
                        PLASMA_Complex32_t *V, int ldv,
                        PLASMA_Complex32_t *T, int ldt,
                        PLASMA_Complex32_t *WORK, int ldwork);
int  CORE_ctsmqr(int side, int trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 PLASMA_Complex32_t *A1, int LDA1,
                 PLASMA_Complex32_t *A2, int LDA2,
                 PLASMA_Complex32_t *V, int LDV,
                 PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *WORK, int LDWORK);
int CORE_ctsmqr_corner( int m1, int n1, int m2, int n2, int m3, int n3,
                        int k, int ib, int nb,
                        PLASMA_Complex32_t *A1, int lda1,
                        PLASMA_Complex32_t *A2, int lda2,
                        PLASMA_Complex32_t *A3, int lda3,
                        PLASMA_Complex32_t *V, int ldv,
                        PLASMA_Complex32_t *T, int ldt,
                        PLASMA_Complex32_t *WORK, int ldwork);
int CORE_ctsmqr_hetra1( int side, int trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        PLASMA_Complex32_t *A1, int lda1,
                        PLASMA_Complex32_t *A2, int lda2,
                        PLASMA_Complex32_t *V, int ldv,
                        PLASMA_Complex32_t *T, int ldt,
                        PLASMA_Complex32_t *WORK, int ldwork);
int  CORE_ctsqrt(int M, int N, int IB,
                 PLASMA_Complex32_t *A1, int LDA1,
                 PLASMA_Complex32_t *A2, int LDA2,
                 PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *TAU, PLASMA_Complex32_t *WORK);
int  CORE_ctsrfb(int side, int trans, int direct, int storev,
                 int M1, int N1, int M2, int N2, int IB,
                 PLASMA_Complex32_t *A1, int LDA1,
                 PLASMA_Complex32_t *A2, int LDA2,
                 PLASMA_Complex32_t *V, int LDV,
                 PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *WORK, int LDWORK);
int  CORE_ctstrf(int M, int N, int IB, int NB,
                 PLASMA_Complex32_t *U, int LDU,
                 PLASMA_Complex32_t *A, int LDA,
                 PLASMA_Complex32_t *L, int LDL,
                 int *IPIV, PLASMA_Complex32_t *WORK,
                 int LDWORK, int *INFO);
int  CORE_cttmqr(int side, int trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 PLASMA_Complex32_t *A1, int LDA1,
                 PLASMA_Complex32_t *A2, int LDA2,
                 PLASMA_Complex32_t *V, int LDV,
                 PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *WORK, int LDWORK);
int  CORE_cttqrt(int M, int N, int IB,
                 PLASMA_Complex32_t *A1, int LDA1,
                 PLASMA_Complex32_t *A2, int LDA2,
                 PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *TAU,
                 PLASMA_Complex32_t *WORK);
int  CORE_cttmlq(int side, int trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 PLASMA_Complex32_t *A1, int LDA1,
                 PLASMA_Complex32_t *A2, int LDA2,
                 PLASMA_Complex32_t *V, int LDV,
                 PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *WORK, int LDWORK);
int  CORE_cttlqt(int M, int N, int IB,
                 PLASMA_Complex32_t *A1, int LDA1,
                 PLASMA_Complex32_t *A2, int LDA2,
                 PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *TAU,
                 PLASMA_Complex32_t *WORK);
int  CORE_cttrfb(int side, int trans, int direct, int storev,
                 int M1, int N1, int M2, int N2, int IB,
                 PLASMA_Complex32_t *A1, int LDA1,
                 PLASMA_Complex32_t *A2, int LDA2,
                 PLASMA_Complex32_t *V, int LDV,
                 PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *WORK, int LDWORK);
int  CORE_cunmlq(int side, int trans,
                 int M, int N, int IB, int K,
                 PLASMA_Complex32_t *V, int LDV,
                 PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *C, int LDC,
                 PLASMA_Complex32_t *WORK, int LDWORK);
int  CORE_cunmqr(int side, int trans,
                 int M, int N, int K, int IB,
                 PLASMA_Complex32_t *V, int LDV,
                 PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *C, int LDC,
                 PLASMA_Complex32_t *WORK, int LDWORK);

/** ****************************************************************************
 *  Declarations of QUARK wrappers (called by PLASMA) - alphabetical order
 **/
void QUARK_CORE_scasum(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum storev, PLASMA_enum uplo, int m, int n,
                       PLASMA_Complex32_t *A, int lda, int szeA,
                       float *work, int szeW);
void QUARK_CORE_scasum_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum storev, PLASMA_enum uplo, int m, int n,
                          PLASMA_Complex32_t *A, int lda, int szeA,
                          float *work, int szeW,
                          float *fake, int szeF);
void QUARK_CORE_caxpy(Quark *quark, Quark_Task_Flags *task_flags,
                      int m, int n, int nb, PLASMA_Complex32_t alpha,
                      PLASMA_Complex32_t *A, int lda,
                      PLASMA_Complex32_t *B, int ldb);
void QUARK_CORE_cbrdalg(Quark *quark, Quark_Task_Flags *task_flags,
                        int uplo,
                        int N, int NB,
                        PLASMA_desc *A,
                        PLASMA_Complex32_t *C,
                        PLASMA_Complex32_t *S,
                        int i, int j, int m, int grsiz, int BAND,
                        int *PCOL, int *ACOL, int *MCOL);
void QUARK_CORE_cgelqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       PLASMA_Complex32_t *A, int lda,
                       PLASMA_Complex32_t *T, int ldt);
void QUARK_CORE_cgemm(Quark *quark, Quark_Task_Flags *task_flags,
                      int transA, int transB,
                      int m, int n, int k, int nb,
                      PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int lda,
                      PLASMA_Complex32_t *B, int ldb,
                      PLASMA_Complex32_t beta, PLASMA_Complex32_t *C, int ldc);
void QUARK_CORE_cgemm2( Quark *quark, Quark_Task_Flags *task_flags,
                        int transA, int transB,
                        int m, int n, int k, int nb,
                        PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int lda,
                        PLASMA_Complex32_t *B, int ldb,
                        PLASMA_Complex32_t beta, PLASMA_Complex32_t *C, int ldc);
void QUARK_CORE_cgemm_f2(Quark *quark, Quark_Task_Flags *task_flags,
                         int transA, int transB,
                         int m, int n, int k, int nb,
                         PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int lda,
                         PLASMA_Complex32_t *B, int ldb,
                         PLASMA_Complex32_t beta, PLASMA_Complex32_t *C, int ldc,
                         PLASMA_Complex32_t *fake1, int szefake1, int flag1,
                         PLASMA_Complex32_t *fake2, int szefake2, int flag2);
void QUARK_CORE_cgemm_p2(Quark *quark, Quark_Task_Flags *task_flags,
                         int transA, int transB,
                         int m, int n, int k, int nb,
                         PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int lda,
                         PLASMA_Complex32_t **B, int ldb,
                         PLASMA_Complex32_t beta, PLASMA_Complex32_t *C, int ldc);
void QUARK_CORE_cgemm_p2f1(Quark *quark, Quark_Task_Flags *task_flags,
                           int transA, int transB,
                           int m, int n, int k, int nb,
                           PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int lda,
                           PLASMA_Complex32_t **B, int ldb,
                           PLASMA_Complex32_t beta, PLASMA_Complex32_t *C, int ldc,
                           PLASMA_Complex32_t *fake1, int szefake1, int flag1);
void QUARK_CORE_cgemm_p3(Quark *quark, Quark_Task_Flags *task_flags,
                         int transA, int transB,
                         int m, int n, int k, int nb,
                         PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int lda,
                         PLASMA_Complex32_t *B, int ldb,
                         PLASMA_Complex32_t beta, PLASMA_Complex32_t **C, int ldc);
void QUARK_CORE_cgeqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       PLASMA_Complex32_t *A, int lda,
                       PLASMA_Complex32_t *T, int ldt);
void QUARK_CORE_cgessm(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int k, int ib, int nb,
                       int *IPIV,
                       PLASMA_Complex32_t *L, int ldl,
                       PLASMA_Complex32_t *A, int lda);
void QUARK_CORE_cgetrf(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int nb,
                       PLASMA_Complex32_t *A, int lda,
                       int *IPIV,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       PLASMA_bool check_info, int iinfo);
void QUARK_CORE_cgetrf_incpiv(Quark *quark, Quark_Task_Flags *task_flags,
                              int m, int n, int ib, int nb,
                              PLASMA_Complex32_t *A, int lda,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo);
void QUARK_CORE_cgetrf_reclap(Quark *quark, Quark_Task_Flags *task_flags,
                              int m, int n, int nb,
                              PLASMA_Complex32_t *A, int lda,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo,
                              int nbthread);
void QUARK_CORE_cgetrf_rectil(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc A, PLASMA_Complex32_t *Amn, int size,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo,
                              int nbthread);
void QUARK_CORE_cgetrip(Quark *quark, Quark_Task_Flags *task_flags,
                        int m, int n, PLASMA_Complex32_t *A, int szeA);
void QUARK_CORE_cgetrip_f1(Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n, PLASMA_Complex32_t *A, int szeA,
                           PLASMA_Complex32_t *fake, int szeF, int paramF);
void QUARK_CORE_cgetrip_f2(Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n, PLASMA_Complex32_t *A, int szeA,
                           PLASMA_Complex32_t *fake1, int szeF1, int paramF1,
                           PLASMA_Complex32_t *fake2, int szeF2, int paramF2);
void QUARK_CORE_chemm(Quark *quark, Quark_Task_Flags *task_flags,
                      int side, int uplo,
                      int m, int n, int nb,
                      PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int lda,
                      PLASMA_Complex32_t *B, int ldb,
                      PLASMA_Complex32_t beta, PLASMA_Complex32_t *C, int ldc);
void QUARK_CORE_chegst(Quark *quark, Quark_Task_Flags *task_flags,
                       int itype, int uplo, int N,
                       PLASMA_Complex32_t *A, int LDA,
                       PLASMA_Complex32_t *B, int LDB,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo);
void QUARK_CORE_cherk(Quark *quark, Quark_Task_Flags *task_flags,
                      int uplo, int trans,
                      int n, int k, int nb,
                      float alpha, PLASMA_Complex32_t *A, int lda,
                      float beta, PLASMA_Complex32_t *C, int ldc);
void QUARK_CORE_cher2k(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int trans,
                       int n, int k, int nb,
                       PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int lda,
                       PLASMA_Complex32_t *B, int LDB,
                       float beta, PLASMA_Complex32_t *C, int ldc);
void QUARK_CORE_cherfb(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo,
                       int n, int k, int ib, int nb,
                       PLASMA_Complex32_t *A, int lda,
                       PLASMA_Complex32_t *T, int ldt,
                       PLASMA_Complex32_t *C, int ldc);
void QUARK_CORE_clacpy(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int m, int n, int mb,
                       PLASMA_Complex32_t *A, int lda,
                       PLASMA_Complex32_t *B, int ldb);
void QUARK_CORE_clange(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, int M, int N,
                       PLASMA_Complex32_t *A, int LDA, int szeA,
                       int szeW, float *result);
void QUARK_CORE_clange_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          int norm, int M, int N,
                          PLASMA_Complex32_t *A, int LDA, int szeA,
                          int szeW, float *result,
                          float *fake, int szeF);
#ifdef COMPLEX
void QUARK_CORE_clanhe(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, int uplo, int N,
                       PLASMA_Complex32_t *A, int LDA, int szeA,
                       int szeW, float *result);
void QUARK_CORE_clanhe_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          int norm, int uplo, int N,
                          PLASMA_Complex32_t *A, int LDA, int szeA,
                          int szeW, float *result,
                          float *fake, int szeF);
#endif
void QUARK_CORE_clansy(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, int uplo, int N,
                       PLASMA_Complex32_t *A, int LDA, int szeA,
                       int szeW, float *result);
void QUARK_CORE_clansy_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          int norm, int uplo, int N,
                          PLASMA_Complex32_t *A, int LDA, int szeA,
                          int szeW, float *result,
                          float *fake, int szeF);
void QUARK_CORE_claset(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int n1, int n2, PLASMA_Complex32_t alpha,
                       PLASMA_Complex32_t beta, PLASMA_Complex32_t *tileA, int ldtilea);
void QUARK_CORE_claset2(Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum uplo, int n1, int n2, PLASMA_Complex32_t alpha,
                        PLASMA_Complex32_t *tileA, int ldtilea);
void QUARK_CORE_claswp(Quark *quark, Quark_Task_Flags *task_flags,
                       int n, PLASMA_Complex32_t *A, int lda, 
                       int i1,  int i2, int *ipiv, int inc);
void QUARK_CORE_claswp_f2(Quark *quark, Quark_Task_Flags *task_flags,
                          int n, PLASMA_Complex32_t *A, int lda, 
                          int i1,  int i2, int *ipiv, int inc,
                          PLASMA_Complex32_t *fake1, int szefake1, int flag1,
                          PLASMA_Complex32_t *fake2, int szefake2, int flag2);
void QUARK_CORE_claswp_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, PLASMA_Complex32_t *A, 
                              int i1,  int i2, int *ipiv, int inc, PLASMA_Complex32_t *fakepanel);
void QUARK_CORE_claswp_ontile_f2(Quark *quark, Quark_Task_Flags *task_flags,
                                 PLASMA_desc descA, PLASMA_Complex32_t *A, 
                                 int i1,  int i2, int *ipiv, int inc, 
                                 PLASMA_Complex32_t *fake1, int szefake1, int flag1,
                                 PLASMA_Complex32_t *fake2, int szefake2, int flag2);
void QUARK_CORE_clauum(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int n, int nb,
                       PLASMA_Complex32_t *A, int lda);
void QUARK_CORE_cplghe(Quark *quark, Quark_Task_Flags *task_flags,
                       float bump, int m, int n, PLASMA_Complex32_t *A, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void QUARK_CORE_cplgsy(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_Complex32_t bump, int m, int n, PLASMA_Complex32_t *A, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void QUARK_CORE_cplrnt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, PLASMA_Complex32_t *A, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void QUARK_CORE_cpotrf(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int n, int nb,
                       PLASMA_Complex32_t *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo);
void QUARK_CORE_cshift( Quark *quark, Quark_Task_Flags *task_flags,
                        int s, int m, int n, int L,
                        PLASMA_Complex32_t *A);
void QUARK_CORE_cshiftw(Quark *quark, Quark_Task_Flags *task_flags,
                        int s, int cl, int m, int n, int L,
                        PLASMA_Complex32_t *A, PLASMA_Complex32_t *W);
void QUARK_CORE_cssssm(Quark *quark, Quark_Task_Flags *task_flags,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       PLASMA_Complex32_t *A1, int lda1,
                       PLASMA_Complex32_t *A2, int lda2,
                       PLASMA_Complex32_t *L1, int ldl1,
                       PLASMA_Complex32_t *L2, int ldl2,
                       int *IPIV);
void QUARK_CORE_csymm(Quark *quark, Quark_Task_Flags *task_flags,
                      int side, int uplo,
                      int m, int n, int nb,
                      PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int lda,
                      PLASMA_Complex32_t *B, int ldb,
                      PLASMA_Complex32_t beta, PLASMA_Complex32_t *C, int ldc);
void QUARK_CORE_csyrk(Quark *quark, Quark_Task_Flags *task_flags,
                      int uplo, int trans,
                      int n, int k, int nb,
                      PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int lda,
                      PLASMA_Complex32_t beta, PLASMA_Complex32_t *C, int ldc);
void QUARK_CORE_csyr2k(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int trans,
                       int n, int k, int nb,
                       PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int lda,
                       PLASMA_Complex32_t *B, int LDB,
                       PLASMA_Complex32_t beta, PLASMA_Complex32_t *C, int ldc);
void QUARK_CORE_cswpab(Quark *quark, Quark_Task_Flags *task_flags,
                       int i, int n1, int n2,
                       PLASMA_Complex32_t *A, int szeA);
void QUARK_CORE_cswptr_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, PLASMA_Complex32_t *Aij, 
                              int i1,  int i2, int *ipiv, int inc, 
                              PLASMA_Complex32_t *Akk, int ldak);
void QUARK_CORE_ctrdalg(Quark *quark, Quark_Task_Flags *task_flags,
                        int uplo,
                        int N, int NB,
                        PLASMA_desc *A,
                        PLASMA_Complex32_t *C,
                        PLASMA_Complex32_t *S,
                        int i, int j, int m, int grsiz, int BAND,
                        int *PCOL, int *ACOL, int *MCOL);
void QUARK_CORE_ctrmm(Quark *quark, Quark_Task_Flags *task_flags,
                      int side, int uplo, int transA, int diag,
                      int m, int n, int nb,
                      PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int lda,
                      PLASMA_Complex32_t *B, int ldb);
void QUARK_CORE_ctrmm_p2(Quark *quark, Quark_Task_Flags *task_flags,
                         int side, int uplo, int transA, int diag,
                         int m, int n, int nb,
                         PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int lda,
                         PLASMA_Complex32_t **B, int ldb);
void QUARK_CORE_ctrsm(Quark *quark, Quark_Task_Flags *task_flags,
                      int side, int uplo, int transA, int diag,
                      int m, int n, int nb,
                      PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int lda,
                      PLASMA_Complex32_t *B, int ldb);
void QUARK_CORE_ctrtri(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int diag, int n, int nb,
                       PLASMA_Complex32_t *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo);
void QUARK_CORE_ctslqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       PLASMA_Complex32_t *A1, int lda1,
                       PLASMA_Complex32_t *A2, int lda2,
                       PLASMA_Complex32_t *T, int ldt);
void QUARK_CORE_ctsmlq(Quark *quark, Quark_Task_Flags *task_flags,
                       int side, int trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       PLASMA_Complex32_t *A1, int lda1,
                       PLASMA_Complex32_t *A2, int lda2,
                       PLASMA_Complex32_t *V, int ldv,
                       PLASMA_Complex32_t *T, int ldt);
void QUARK_CORE_ctsmlq_hetra1(Quark *quark, Quark_Task_Flags *task_flags,
                              int side, int trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              PLASMA_Complex32_t *A1, int lda1,
                              PLASMA_Complex32_t *A2, int lda2,
                              PLASMA_Complex32_t *V, int ldv,
                              PLASMA_Complex32_t *T, int ldt);
void QUARK_CORE_ctsmlq_corner(Quark *quark, Quark_Task_Flags *task_flags,
                              int m1, int n1, int m2, int n2, int m3, int n3, int k, int ib, int nb,
                              PLASMA_Complex32_t *A1, int lda1,
                              PLASMA_Complex32_t *A2, int lda2,
                              PLASMA_Complex32_t *A3, int lda3,
                              PLASMA_Complex32_t *V, int ldv,
                              PLASMA_Complex32_t *T, int ldt);
void QUARK_CORE_ctsmqr(Quark *quark, Quark_Task_Flags *task_flags,
                       int side, int trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       PLASMA_Complex32_t *A1, int lda1,
                       PLASMA_Complex32_t *A2, int lda2,
                       PLASMA_Complex32_t *V, int ldv,
                       PLASMA_Complex32_t *T, int ldt);
void QUARK_CORE_ctsmqr_hetra1(Quark *quark, Quark_Task_Flags *task_flags,
                              int side, int trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              PLASMA_Complex32_t *A1, int lda1,
                              PLASMA_Complex32_t *A2, int lda2,
                              PLASMA_Complex32_t *V, int ldv,
                              PLASMA_Complex32_t *T, int ldt);
void QUARK_CORE_ctsmqr_corner(Quark *quark, Quark_Task_Flags *task_flags,
                              int m1, int n1, int m2, int n2, int m3, int n3, int k, int ib, int nb,
                              PLASMA_Complex32_t *A1, int lda1,
                              PLASMA_Complex32_t *A2, int lda2,
                              PLASMA_Complex32_t *A3, int lda3,
                              PLASMA_Complex32_t *V, int ldv,
                              PLASMA_Complex32_t *T, int ldt);
void QUARK_CORE_ctsqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       PLASMA_Complex32_t *A1, int lda1,
                       PLASMA_Complex32_t *A2, int lda2,
                       PLASMA_Complex32_t *T, int ldt);
void QUARK_CORE_ctstrf(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       PLASMA_Complex32_t *U, int ldu,
                       PLASMA_Complex32_t *A, int lda,
                       PLASMA_Complex32_t *L, int ldl,
                       int *IPIV,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       PLASMA_bool check_info, int iinfo);
void QUARK_CORE_cttmqr(Quark *quark, Quark_Task_Flags *task_flags,
                       int side, int trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       PLASMA_Complex32_t *A1, int lda1,
                       PLASMA_Complex32_t *A2, int lda2,
                       PLASMA_Complex32_t *V, int ldv,
                       PLASMA_Complex32_t *T, int ldt);
void QUARK_CORE_cttqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       PLASMA_Complex32_t *A1, int lda1,
                       PLASMA_Complex32_t *A2, int lda2,
                       PLASMA_Complex32_t *T, int ldt);
void QUARK_CORE_cttmlq(Quark *quark, Quark_Task_Flags *task_flags,
                       int side, int trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       PLASMA_Complex32_t *A1, int lda1,
                       PLASMA_Complex32_t *A2, int lda2,
                       PLASMA_Complex32_t *V, int ldv,
                       PLASMA_Complex32_t *T, int ldt);
void QUARK_CORE_cttlqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       PLASMA_Complex32_t *A1, int lda1,
                       PLASMA_Complex32_t *A2, int lda2,
                       PLASMA_Complex32_t *T, int ldt);
void QUARK_CORE_cunmlq(Quark *quark, Quark_Task_Flags *task_flags,
                       int side, int trans,
                       int m, int n, int ib,  int nb, int k,
                       PLASMA_Complex32_t *A, int lda,
                       PLASMA_Complex32_t *T, int ldt,
                       PLASMA_Complex32_t *C, int ldc);
void QUARK_CORE_cunmqr(Quark *quark, Quark_Task_Flags *task_flags,
                       int side, int trans,
                       int m, int n, int k, int ib, int nb,
                       PLASMA_Complex32_t *A, int lda,
                       PLASMA_Complex32_t *T, int ldt,
                       PLASMA_Complex32_t *C, int ldc);

/** ****************************************************************************
 *  Declarations of QUARK wrappers (called by QUARK) - alphabetical order
 **/
void CORE_scasum_quark(Quark *quark);
void CORE_scasum_f1_quark(Quark *quark);
void CORE_caxpy_quark(Quark *quark);
void CORE_cbrdalg_quark(Quark *quark);
void CORE_cgelqt_quark(Quark *quark);
void CORE_cgemm_quark(Quark *quark);
void CORE_cgeqrt_quark(Quark *quark);
void CORE_cgessm_quark(Quark *quark);
void CORE_cgetrf_quark(Quark *quark);
void CORE_cgetrf_incpiv_quark(Quark *quark);
void CORE_cgetrf_reclap_quark(Quark *quark);
void CORE_cgetrf_rectil_quark(Quark* quark);
void CORE_cgetrip_quark(Quark *quark);
void CORE_cgetrip_f1_quark(Quark *quark);
void CORE_cgetrip_f2_quark(Quark *quark);
#ifdef COMPLEX
void CORE_chemm_quark(Quark *quark);
void CORE_cherk_quark(Quark *quark);
void CORE_cher2k_quark(Quark *quark);
#endif
void CORE_chegst_quark(Quark *quark);
void CORE_cherfb_quark(Quark *quark);
void CORE_clacpy_quark(Quark *quark);
void CORE_clange_quark(Quark *quark);
void CORE_clange_f1_quark(Quark *quark);
#ifdef COMPLEX
void CORE_clanhe_quark(Quark *quark);
void CORE_clanhe_f1_quark(Quark *quark);
#endif
void CORE_clansy_quark(Quark *quark);
void CORE_clansy_f1_quark(Quark *quark);
void CORE_claset_quark(Quark *quark);
void CORE_claset2_quark(Quark *quark);
void CORE_clauum_quark(Quark *quark);
void CORE_cplghe_quark(Quark *quark);
void CORE_cplgsy_quark(Quark *quark);
void CORE_cplrnt_quark(Quark *quark);
void CORE_cpotrf_quark(Quark *quark);
void CORE_cshift_quark(Quark *quark);
void CORE_cshiftw_quark(Quark *quark);
void CORE_cssssm_quark(Quark *quark);
void CORE_csymm_quark(Quark *quark);
void CORE_csyrk_quark(Quark *quark);
void CORE_csyr2k_quark(Quark *quark);
void CORE_cswpab_quark(Quark *quark);
void CORE_cswptr_ontile_quark(Quark *quark);
void CORE_ctrdalg_quark(Quark *quark);
void CORE_ctrmm_quark(Quark *quark);
void CORE_ctrsm_quark(Quark *quark);
void CORE_ctrtri_quark(Quark *quark);
void CORE_ctslqt_quark(Quark *quark);
void CORE_ctsmlq_quark(Quark *quark);
void CORE_ctsmlq_hetra1_quark(Quark *quark);
void CORE_ctsmlq_corner_quark(Quark *quark);
void CORE_ctsmqr_quark(Quark *quark);
void CORE_ctsmqr_hetra1_quark(Quark *quark);
void CORE_ctsmqr_corner_quark(Quark *quark);
void CORE_ctsqrt_quark(Quark *quark);
void CORE_ctstrf_quark(Quark *quark);
void CORE_cttmqr_quark(Quark *quark);
void CORE_cttqrt_quark(Quark *quark);
void CORE_cttmlq_quark(Quark *quark);
void CORE_cttlqt_quark(Quark *quark);
void CORE_cunmlq_quark(Quark *quark);
void CORE_cunmqr_quark(Quark *quark);

void CORE_claswp_quark(Quark* quark);
void CORE_claswp_f2_quark(Quark* quark);
void CORE_claswp_ontile_quark(Quark *quark);
void CORE_claswp_ontile_f2_quark(Quark *quark);
void CORE_ctrmm_p2_quark(Quark* quark);
void CORE_cgemm_f2_quark(Quark* quark);
void CORE_cgemm_p2_quark(Quark* quark);
void CORE_cgemm_p2f1_quark(Quark* quark);
void CORE_cgemm_p3_quark(Quark* quark);





void CORE_ctrdalg_v2_quark(Quark *quark);


void QUARK_CORE_ctrdalg_v2(Quark *quark, Quark_Task_Flags *task_flags,
                        int uplo,
                        PLASMA_desc *A,
                        PLASMA_Complex32_t *C,
                        PLASMA_Complex32_t *S,
                        int grsiz, int s, int id, int blktile);

void CORE_ctrdalg_v2(PLASMA_enum uplo,  
                  PLASMA_desc *pA, PLASMA_Complex32_t *C, PLASMA_Complex32_t *S,
                  int grsiz, int s, int id, int blktile);


#ifdef __cplusplus
}
#endif

#endif
