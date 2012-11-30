/**
 *
 * @file core_sblas.h
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
 * @generated s Thu Sep 15 12:08:49 2011
 *
 **/
#ifndef _PLASMA_CORE_SBLAS_H_
#define _PLASMA_CORE_SBLAS_H_
#define REAL

#ifdef __cplusplus
extern "C" {
#endif

/** ****************************************************************************
 *  Declarations of serial kernels - alphabetical order
 **/

int CORE_slarfx2(int side, int N,
                 float V,
                 float TAU,
                 float *C1, int LDC1,
                 float *C2, int LDC2);
int CORE_slarfx2c(int uplo,
                  float V,
                  float TAU,
                  float *C1,
                  float *C2,
                  float *C3);
int CORE_slarfx2ce(int uplo,
                float *V,
                float *TAU,
                float *C1,
                float *C2,
                float *C3);
int CORE_shbelr(int uplo, int N,
                PLASMA_desc *A,
                float *V,
                float *TAU,
                int st,
                int ed,
                int eltsize);
int CORE_shbrce(int uplo, int N,
                PLASMA_desc *A,
                float *V,
                float *TAU,
                int st,
                int ed,
                int eltsize);
int CORE_shblrx(int uplo, int N,
                PLASMA_desc *A,
                float *V,
                float *TAU,
                int st,
                int ed,
                int eltsize);
int CORE_sgbelr(int uplo, int N,
                PLASMA_desc *A,
                float *V,
                float *TAU,
                int st,
                int ed,
                int eltsize);
int CORE_sgbrce(int uplo, int N,
                PLASMA_desc *A,
                float *V,
                float *TAU,
                int st,
                int ed,
                int eltsize);
int CORE_sgblrx(int uplo, int N,
                PLASMA_desc *A,
                float *V,
                float *TAU,
                int st,
                int ed,
                int eltsize);
void CORE_sasum(int storev, int uplo, int M, int N,
                 float *A, int lda, float *work);
void CORE_saxpy(int M, int N,  float alpha,
                float *A, int LDA,
                float *B, int LDB);
void CORE_sbrdalg(PLASMA_enum uplo, int N, int NB, 
                  PLASMA_desc *pA, float *C, float *S,
                  int i, int j, int m, int grsiz);
int  CORE_sgelqt(int M, int N, int IB,
                 float *A, int LDA,
                 float *T, int LDT,
                 float *TAU, float *WORK);
void CORE_sgemm(int transA, int transB,
                int M, int N, int K,
                float alpha, float *A, int LDA,
                float *B, int LDB,
                float beta, float *C, int LDC);
int  CORE_sgeqrt(int M, int N, int IB,
                 float *A, int LDA,
                 float *T, int LDT,
                 float *TAU, float *WORK);
int  CORE_sgessm(int M, int N, int K, int IB,
                 int *IPIV,
                 float *L, int LDL,
                 float *A, int LDA);
int  CORE_sgetrf(int M, int N, 
                 float *A, int LDA,
                 int *IPIV, int *INFO);
int  CORE_sgetrf_incpiv(int M, int N, int IB,
                        float *A, int LDA,
                        int *IPIV, int *INFO);
int  CORE_sgetrf_reclap(const int M, const int N,
                        float *A, const int LDA,
                        int *IPIV, int *info);
int  CORE_sgetrf_rectil(const PLASMA_desc A, int *IPIV, int *info);
void CORE_sgetrip(int m, int n, float *A, 
                  float *work);
#ifdef COMPLEX
void CORE_ssygst(int itype, int uplo, int N,
                 float *A, int LDA,
                 float *B, int LDB, int *INFO);
void CORE_ssymm(int side, int uplo,
                int M, int N,
                float alpha, float *A, int LDA,
                float *B, int LDB,
                float beta, float *C, int LDC);
void CORE_ssyrk(int uplo, int trans,
                int N, int K,
                float alpha, float *A, int LDA,
                float beta, float *C, int LDC);
void CORE_ssyr2k(int uplo, int trans,
                 int N, int K,
                 float alpha, float *A, int LDA,
                 float *B, int LDB,
                 float beta, float *C, int LDC);
int  CORE_ssyrfb(PLASMA_enum uplo, int N, int K, int IB, int NB,
                 float *A, int LDA,
                 float *T, int LDT,
                 float *C, int LDC,
                 float *WORK, int LDWORK);
#endif
void CORE_slacpy(PLASMA_enum uplo, int M, int N,
                 float *A, int LDA,
                 float *B, int LDB);
void CORE_slange(int norm, int M, int N,
                 float *A, int LDA,
                 float *work, float *normA);
#ifdef COMPLEX
void CORE_slansy(int norm, int uplo, int N,
                 float *A, int LDA,
                 float *work, float *normA);
#endif
void CORE_slansy(int norm, int uplo, int N,
                 float *A, int LDA,
                 float *work, float *normA);
void CORE_slaset(PLASMA_enum uplo, int n1, int n2, float alpha,
                 float beta, float *tileA, int ldtilea);
void CORE_slaset2(PLASMA_enum uplo, int n1, int n2, float alpha,
                  float *tileA, int ldtilea);
void CORE_slaswp(int N, float *A, int LDA, 
                 int I1,  int I2, int *IPIV, int INC);
int  CORE_slaswp_ontile(PLASMA_desc descA, int i1, int i2, int *ipiv, int inc);
void CORE_slauum(int uplo, int N, float *A, int LDA);
void CORE_splgsy(float bump, int m, int n, float *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
void CORE_splgsy(float bump, int m, int n, float *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
void CORE_splrnt(int m, int n, float *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
void CORE_spotrf(int uplo, int N, float *A, int LDA, int *INFO);
void CORE_sshift(int s, int m, int n, int L,
                 float *A);
void CORE_sshiftw(int s, int cl, int m, int n, int L,
                  float *A, float *W);
int  CORE_sssssm(int M1, int N1, int M2, int N2, int K, int IB,
                 float *A1, int LDA1,
                 float *A2, int LDA2,
                 float *L1, int LDL1,
                 float *L2, int LDL2,
                 int *IPIV);
void CORE_ssymm(int side, int uplo,
                int M, int N,
                float alpha, float *A, int LDA,
                float *B, int LDB,
                float beta, float *C, int LDC);
void CORE_ssyrk(int uplo, int trans,
                int N, int K,
                float alpha, float *A, int LDA,
                float beta, float *C, int LDC);
void CORE_ssyr2k(int uplo, int trans,
                 int N, int K,
                 float alpha, float *A, int LDA,
                 float *B, int LDB,
                 float beta, float *C, int LDC);
void CORE_sswpab(int i, int n1, int n2,
                 float *A, float *work);
int  CORE_sswptr_ontile(PLASMA_desc descA, int i1, int i2, int *ipiv, int inc,
                        float *Akk, int ldak);
void CORE_strdalg(PLASMA_enum uplo, int N, int NB, 
                  PLASMA_desc *pA, float *C, float *S,
                  int i, int j, int m, int grsiz);
void CORE_strmm(int side, int uplo,
                int transA, int diag,
                int M, int N,
                float alpha, float *A, int LDA,
                float *B, int LDB);
void CORE_strsm(int side, int uplo,
                int transA, int diag,
                int M, int N,
                float alpha, float *A, int LDA,
                float *B, int LDB);
void CORE_strtri(int uplo, int diag, int N, float *A, int LDA, int *info);
int  CORE_stslqt(int M, int N, int IB,
                 float *A1, int LDA1,
                 float *A2, int LDA2,
                 float *T, int LDT,
                 float *TAU, float *WORK);
int  CORE_stsmlq(int side, int trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 float *A1, int LDA1,
                 float *A2, int LDA2,
                 float *V, int LDV,
                 float *T, int LDT,
                 float *WORK, int LDWORK);
int CORE_stsmlq_corner( int m1, int n1, int m2, int n2, int m3, int n3,
                        int k, int ib, int nb,
                        float *A1, int lda1,
                        float *A2, int lda2,
                        float *A3, int lda3,
                        float *V, int ldv,
                        float *T, int ldt,
                        float *WORK, int ldwork);
int CORE_stsmlq_sytra1( int side, int trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        float *A1, int lda1,
                        float *A2, int lda2,
                        float *V, int ldv,
                        float *T, int ldt,
                        float *WORK, int ldwork);
int  CORE_stsmqr(int side, int trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 float *A1, int LDA1,
                 float *A2, int LDA2,
                 float *V, int LDV,
                 float *T, int LDT,
                 float *WORK, int LDWORK);
int CORE_stsmqr_corner( int m1, int n1, int m2, int n2, int m3, int n3,
                        int k, int ib, int nb,
                        float *A1, int lda1,
                        float *A2, int lda2,
                        float *A3, int lda3,
                        float *V, int ldv,
                        float *T, int ldt,
                        float *WORK, int ldwork);
int CORE_stsmqr_sytra1( int side, int trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        float *A1, int lda1,
                        float *A2, int lda2,
                        float *V, int ldv,
                        float *T, int ldt,
                        float *WORK, int ldwork);
int  CORE_stsqrt(int M, int N, int IB,
                 float *A1, int LDA1,
                 float *A2, int LDA2,
                 float *T, int LDT,
                 float *TAU, float *WORK);
int  CORE_stsrfb(int side, int trans, int direct, int storev,
                 int M1, int N1, int M2, int N2, int IB,
                 float *A1, int LDA1,
                 float *A2, int LDA2,
                 float *V, int LDV,
                 float *T, int LDT,
                 float *WORK, int LDWORK);
int  CORE_ststrf(int M, int N, int IB, int NB,
                 float *U, int LDU,
                 float *A, int LDA,
                 float *L, int LDL,
                 int *IPIV, float *WORK,
                 int LDWORK, int *INFO);
int  CORE_sttmqr(int side, int trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 float *A1, int LDA1,
                 float *A2, int LDA2,
                 float *V, int LDV,
                 float *T, int LDT,
                 float *WORK, int LDWORK);
int  CORE_sttqrt(int M, int N, int IB,
                 float *A1, int LDA1,
                 float *A2, int LDA2,
                 float *T, int LDT,
                 float *TAU,
                 float *WORK);
int  CORE_sttmlq(int side, int trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 float *A1, int LDA1,
                 float *A2, int LDA2,
                 float *V, int LDV,
                 float *T, int LDT,
                 float *WORK, int LDWORK);
int  CORE_sttlqt(int M, int N, int IB,
                 float *A1, int LDA1,
                 float *A2, int LDA2,
                 float *T, int LDT,
                 float *TAU,
                 float *WORK);
int  CORE_sttrfb(int side, int trans, int direct, int storev,
                 int M1, int N1, int M2, int N2, int IB,
                 float *A1, int LDA1,
                 float *A2, int LDA2,
                 float *V, int LDV,
                 float *T, int LDT,
                 float *WORK, int LDWORK);
int  CORE_sormlq(int side, int trans,
                 int M, int N, int IB, int K,
                 float *V, int LDV,
                 float *T, int LDT,
                 float *C, int LDC,
                 float *WORK, int LDWORK);
int  CORE_sormqr(int side, int trans,
                 int M, int N, int K, int IB,
                 float *V, int LDV,
                 float *T, int LDT,
                 float *C, int LDC,
                 float *WORK, int LDWORK);

/** ****************************************************************************
 *  Declarations of QUARK wrappers (called by PLASMA) - alphabetical order
 **/
void QUARK_CORE_sasum(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum storev, PLASMA_enum uplo, int m, int n,
                       float *A, int lda, int szeA,
                       float *work, int szeW);
void QUARK_CORE_sasum_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum storev, PLASMA_enum uplo, int m, int n,
                          float *A, int lda, int szeA,
                          float *work, int szeW,
                          float *fake, int szeF);
void QUARK_CORE_saxpy(Quark *quark, Quark_Task_Flags *task_flags,
                      int m, int n, int nb, float alpha,
                      float *A, int lda,
                      float *B, int ldb);
void QUARK_CORE_sbrdalg(Quark *quark, Quark_Task_Flags *task_flags,
                        int uplo,
                        int N, int NB,
                        PLASMA_desc *A,
                        float *C,
                        float *S,
                        int i, int j, int m, int grsiz, int BAND,
                        int *PCOL, int *ACOL, int *MCOL);
void QUARK_CORE_sgelqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       float *A, int lda,
                       float *T, int ldt);
void QUARK_CORE_sgemm(Quark *quark, Quark_Task_Flags *task_flags,
                      int transA, int transB,
                      int m, int n, int k, int nb,
                      float alpha, float *A, int lda,
                      float *B, int ldb,
                      float beta, float *C, int ldc);
void QUARK_CORE_sgemm2( Quark *quark, Quark_Task_Flags *task_flags,
                        int transA, int transB,
                        int m, int n, int k, int nb,
                        float alpha, float *A, int lda,
                        float *B, int ldb,
                        float beta, float *C, int ldc);
void QUARK_CORE_sgemm_f2(Quark *quark, Quark_Task_Flags *task_flags,
                         int transA, int transB,
                         int m, int n, int k, int nb,
                         float alpha, float *A, int lda,
                         float *B, int ldb,
                         float beta, float *C, int ldc,
                         float *fake1, int szefake1, int flag1,
                         float *fake2, int szefake2, int flag2);
void QUARK_CORE_sgemm_p2(Quark *quark, Quark_Task_Flags *task_flags,
                         int transA, int transB,
                         int m, int n, int k, int nb,
                         float alpha, float *A, int lda,
                         float **B, int ldb,
                         float beta, float *C, int ldc);
void QUARK_CORE_sgemm_p2f1(Quark *quark, Quark_Task_Flags *task_flags,
                           int transA, int transB,
                           int m, int n, int k, int nb,
                           float alpha, float *A, int lda,
                           float **B, int ldb,
                           float beta, float *C, int ldc,
                           float *fake1, int szefake1, int flag1);
void QUARK_CORE_sgemm_p3(Quark *quark, Quark_Task_Flags *task_flags,
                         int transA, int transB,
                         int m, int n, int k, int nb,
                         float alpha, float *A, int lda,
                         float *B, int ldb,
                         float beta, float **C, int ldc);
void QUARK_CORE_sgeqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       float *A, int lda,
                       float *T, int ldt);
void QUARK_CORE_sgessm(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int k, int ib, int nb,
                       int *IPIV,
                       float *L, int ldl,
                       float *A, int lda);
void QUARK_CORE_sgetrf(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int nb,
                       float *A, int lda,
                       int *IPIV,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       PLASMA_bool check_info, int iinfo);
void QUARK_CORE_sgetrf_incpiv(Quark *quark, Quark_Task_Flags *task_flags,
                              int m, int n, int ib, int nb,
                              float *A, int lda,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo);
void QUARK_CORE_sgetrf_reclap(Quark *quark, Quark_Task_Flags *task_flags,
                              int m, int n, int nb,
                              float *A, int lda,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo,
                              int nbthread);
void QUARK_CORE_sgetrf_rectil(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc A, float *Amn, int size,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo,
                              int nbthread);
void QUARK_CORE_sgetrip(Quark *quark, Quark_Task_Flags *task_flags,
                        int m, int n, float *A, int szeA);
void QUARK_CORE_sgetrip_f1(Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n, float *A, int szeA,
                           float *fake, int szeF, int paramF);
void QUARK_CORE_sgetrip_f2(Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n, float *A, int szeA,
                           float *fake1, int szeF1, int paramF1,
                           float *fake2, int szeF2, int paramF2);
void QUARK_CORE_ssymm(Quark *quark, Quark_Task_Flags *task_flags,
                      int side, int uplo,
                      int m, int n, int nb,
                      float alpha, float *A, int lda,
                      float *B, int ldb,
                      float beta, float *C, int ldc);
void QUARK_CORE_ssygst(Quark *quark, Quark_Task_Flags *task_flags,
                       int itype, int uplo, int N,
                       float *A, int LDA,
                       float *B, int LDB,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo);
void QUARK_CORE_ssyrk(Quark *quark, Quark_Task_Flags *task_flags,
                      int uplo, int trans,
                      int n, int k, int nb,
                      float alpha, float *A, int lda,
                      float beta, float *C, int ldc);
void QUARK_CORE_ssyr2k(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int trans,
                       int n, int k, int nb,
                       float alpha, float *A, int lda,
                       float *B, int LDB,
                       float beta, float *C, int ldc);
void QUARK_CORE_ssyrfb(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo,
                       int n, int k, int ib, int nb,
                       float *A, int lda,
                       float *T, int ldt,
                       float *C, int ldc);
void QUARK_CORE_slacpy(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int m, int n, int mb,
                       float *A, int lda,
                       float *B, int ldb);
void QUARK_CORE_slange(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, int M, int N,
                       float *A, int LDA, int szeA,
                       int szeW, float *result);
void QUARK_CORE_slange_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          int norm, int M, int N,
                          float *A, int LDA, int szeA,
                          int szeW, float *result,
                          float *fake, int szeF);
#ifdef COMPLEX
void QUARK_CORE_slansy(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, int uplo, int N,
                       float *A, int LDA, int szeA,
                       int szeW, float *result);
void QUARK_CORE_slansy_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          int norm, int uplo, int N,
                          float *A, int LDA, int szeA,
                          int szeW, float *result,
                          float *fake, int szeF);
#endif
void QUARK_CORE_slansy(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, int uplo, int N,
                       float *A, int LDA, int szeA,
                       int szeW, float *result);
void QUARK_CORE_slansy_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          int norm, int uplo, int N,
                          float *A, int LDA, int szeA,
                          int szeW, float *result,
                          float *fake, int szeF);
void QUARK_CORE_slaset(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int n1, int n2, float alpha,
                       float beta, float *tileA, int ldtilea);
void QUARK_CORE_slaset2(Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum uplo, int n1, int n2, float alpha,
                        float *tileA, int ldtilea);
void QUARK_CORE_slaswp(Quark *quark, Quark_Task_Flags *task_flags,
                       int n, float *A, int lda, 
                       int i1,  int i2, int *ipiv, int inc);
void QUARK_CORE_slaswp_f2(Quark *quark, Quark_Task_Flags *task_flags,
                          int n, float *A, int lda, 
                          int i1,  int i2, int *ipiv, int inc,
                          float *fake1, int szefake1, int flag1,
                          float *fake2, int szefake2, int flag2);
void QUARK_CORE_slaswp_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, float *A, 
                              int i1,  int i2, int *ipiv, int inc, float *fakepanel);
void QUARK_CORE_slaswp_ontile_f2(Quark *quark, Quark_Task_Flags *task_flags,
                                 PLASMA_desc descA, float *A, 
                                 int i1,  int i2, int *ipiv, int inc, 
                                 float *fake1, int szefake1, int flag1,
                                 float *fake2, int szefake2, int flag2);
void QUARK_CORE_slauum(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int n, int nb,
                       float *A, int lda);
void QUARK_CORE_splgsy(Quark *quark, Quark_Task_Flags *task_flags,
                       float bump, int m, int n, float *A, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void QUARK_CORE_splgsy(Quark *quark, Quark_Task_Flags *task_flags,
                       float bump, int m, int n, float *A, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void QUARK_CORE_splrnt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, float *A, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void QUARK_CORE_spotrf(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int n, int nb,
                       float *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo);
void QUARK_CORE_sshift( Quark *quark, Quark_Task_Flags *task_flags,
                        int s, int m, int n, int L,
                        float *A);
void QUARK_CORE_sshiftw(Quark *quark, Quark_Task_Flags *task_flags,
                        int s, int cl, int m, int n, int L,
                        float *A, float *W);
void QUARK_CORE_sssssm(Quark *quark, Quark_Task_Flags *task_flags,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       float *A1, int lda1,
                       float *A2, int lda2,
                       float *L1, int ldl1,
                       float *L2, int ldl2,
                       int *IPIV);
void QUARK_CORE_ssymm(Quark *quark, Quark_Task_Flags *task_flags,
                      int side, int uplo,
                      int m, int n, int nb,
                      float alpha, float *A, int lda,
                      float *B, int ldb,
                      float beta, float *C, int ldc);
void QUARK_CORE_ssyrk(Quark *quark, Quark_Task_Flags *task_flags,
                      int uplo, int trans,
                      int n, int k, int nb,
                      float alpha, float *A, int lda,
                      float beta, float *C, int ldc);
void QUARK_CORE_ssyr2k(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int trans,
                       int n, int k, int nb,
                       float alpha, float *A, int lda,
                       float *B, int LDB,
                       float beta, float *C, int ldc);
void QUARK_CORE_sswpab(Quark *quark, Quark_Task_Flags *task_flags,
                       int i, int n1, int n2,
                       float *A, int szeA);
void QUARK_CORE_sswptr_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, float *Aij, 
                              int i1,  int i2, int *ipiv, int inc, 
                              float *Akk, int ldak);
void QUARK_CORE_strdalg(Quark *quark, Quark_Task_Flags *task_flags,
                        int uplo,
                        int N, int NB,
                        PLASMA_desc *A,
                        float *C,
                        float *S,
                        int i, int j, int m, int grsiz, int BAND,
                        int *PCOL, int *ACOL, int *MCOL);
void QUARK_CORE_strmm(Quark *quark, Quark_Task_Flags *task_flags,
                      int side, int uplo, int transA, int diag,
                      int m, int n, int nb,
                      float alpha, float *A, int lda,
                      float *B, int ldb);
void QUARK_CORE_strmm_p2(Quark *quark, Quark_Task_Flags *task_flags,
                         int side, int uplo, int transA, int diag,
                         int m, int n, int nb,
                         float alpha, float *A, int lda,
                         float **B, int ldb);
void QUARK_CORE_strsm(Quark *quark, Quark_Task_Flags *task_flags,
                      int side, int uplo, int transA, int diag,
                      int m, int n, int nb,
                      float alpha, float *A, int lda,
                      float *B, int ldb);
void QUARK_CORE_strtri(Quark *quark, Quark_Task_Flags *task_flags,
                       int uplo, int diag, int n, int nb,
                       float *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo);
void QUARK_CORE_stslqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       float *A1, int lda1,
                       float *A2, int lda2,
                       float *T, int ldt);
void QUARK_CORE_stsmlq(Quark *quark, Quark_Task_Flags *task_flags,
                       int side, int trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       float *A1, int lda1,
                       float *A2, int lda2,
                       float *V, int ldv,
                       float *T, int ldt);
void QUARK_CORE_stsmlq_sytra1(Quark *quark, Quark_Task_Flags *task_flags,
                              int side, int trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              float *A1, int lda1,
                              float *A2, int lda2,
                              float *V, int ldv,
                              float *T, int ldt);
void QUARK_CORE_stsmlq_corner(Quark *quark, Quark_Task_Flags *task_flags,
                              int m1, int n1, int m2, int n2, int m3, int n3, int k, int ib, int nb,
                              float *A1, int lda1,
                              float *A2, int lda2,
                              float *A3, int lda3,
                              float *V, int ldv,
                              float *T, int ldt);
void QUARK_CORE_stsmqr(Quark *quark, Quark_Task_Flags *task_flags,
                       int side, int trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       float *A1, int lda1,
                       float *A2, int lda2,
                       float *V, int ldv,
                       float *T, int ldt);
void QUARK_CORE_stsmqr_sytra1(Quark *quark, Quark_Task_Flags *task_flags,
                              int side, int trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              float *A1, int lda1,
                              float *A2, int lda2,
                              float *V, int ldv,
                              float *T, int ldt);
void QUARK_CORE_stsmqr_corner(Quark *quark, Quark_Task_Flags *task_flags,
                              int m1, int n1, int m2, int n2, int m3, int n3, int k, int ib, int nb,
                              float *A1, int lda1,
                              float *A2, int lda2,
                              float *A3, int lda3,
                              float *V, int ldv,
                              float *T, int ldt);
void QUARK_CORE_stsqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       float *A1, int lda1,
                       float *A2, int lda2,
                       float *T, int ldt);
void QUARK_CORE_ststrf(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       float *U, int ldu,
                       float *A, int lda,
                       float *L, int ldl,
                       int *IPIV,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       PLASMA_bool check_info, int iinfo);
void QUARK_CORE_sttmqr(Quark *quark, Quark_Task_Flags *task_flags,
                       int side, int trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       float *A1, int lda1,
                       float *A2, int lda2,
                       float *V, int ldv,
                       float *T, int ldt);
void QUARK_CORE_sttqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       float *A1, int lda1,
                       float *A2, int lda2,
                       float *T, int ldt);
void QUARK_CORE_sttmlq(Quark *quark, Quark_Task_Flags *task_flags,
                       int side, int trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       float *A1, int lda1,
                       float *A2, int lda2,
                       float *V, int ldv,
                       float *T, int ldt);
void QUARK_CORE_sttlqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       float *A1, int lda1,
                       float *A2, int lda2,
                       float *T, int ldt);
void QUARK_CORE_sormlq(Quark *quark, Quark_Task_Flags *task_flags,
                       int side, int trans,
                       int m, int n, int ib,  int nb, int k,
                       float *A, int lda,
                       float *T, int ldt,
                       float *C, int ldc);
void QUARK_CORE_sormqr(Quark *quark, Quark_Task_Flags *task_flags,
                       int side, int trans,
                       int m, int n, int k, int ib, int nb,
                       float *A, int lda,
                       float *T, int ldt,
                       float *C, int ldc);

/** ****************************************************************************
 *  Declarations of QUARK wrappers (called by QUARK) - alphabetical order
 **/
void CORE_sasum_quark(Quark *quark);
void CORE_sasum_f1_quark(Quark *quark);
void CORE_saxpy_quark(Quark *quark);
void CORE_sbrdalg_quark(Quark *quark);
void CORE_sgelqt_quark(Quark *quark);
void CORE_sgemm_quark(Quark *quark);
void CORE_sgeqrt_quark(Quark *quark);
void CORE_sgessm_quark(Quark *quark);
void CORE_sgetrf_quark(Quark *quark);
void CORE_sgetrf_incpiv_quark(Quark *quark);
void CORE_sgetrf_reclap_quark(Quark *quark);
void CORE_sgetrf_rectil_quark(Quark* quark);
void CORE_sgetrip_quark(Quark *quark);
void CORE_sgetrip_f1_quark(Quark *quark);
void CORE_sgetrip_f2_quark(Quark *quark);
#ifdef COMPLEX
void CORE_ssymm_quark(Quark *quark);
void CORE_ssyrk_quark(Quark *quark);
void CORE_ssyr2k_quark(Quark *quark);
#endif
void CORE_ssygst_quark(Quark *quark);
void CORE_ssyrfb_quark(Quark *quark);
void CORE_slacpy_quark(Quark *quark);
void CORE_slange_quark(Quark *quark);
void CORE_slange_f1_quark(Quark *quark);
#ifdef COMPLEX
void CORE_slansy_quark(Quark *quark);
void CORE_slansy_f1_quark(Quark *quark);
#endif
void CORE_slansy_quark(Quark *quark);
void CORE_slansy_f1_quark(Quark *quark);
void CORE_slaset_quark(Quark *quark);
void CORE_slaset2_quark(Quark *quark);
void CORE_slauum_quark(Quark *quark);
void CORE_splgsy_quark(Quark *quark);
void CORE_splgsy_quark(Quark *quark);
void CORE_splrnt_quark(Quark *quark);
void CORE_spotrf_quark(Quark *quark);
void CORE_sshift_quark(Quark *quark);
void CORE_sshiftw_quark(Quark *quark);
void CORE_sssssm_quark(Quark *quark);
void CORE_ssymm_quark(Quark *quark);
void CORE_ssyrk_quark(Quark *quark);
void CORE_ssyr2k_quark(Quark *quark);
void CORE_sswpab_quark(Quark *quark);
void CORE_sswptr_ontile_quark(Quark *quark);
void CORE_strdalg_quark(Quark *quark);
void CORE_strmm_quark(Quark *quark);
void CORE_strsm_quark(Quark *quark);
void CORE_strtri_quark(Quark *quark);
void CORE_stslqt_quark(Quark *quark);
void CORE_stsmlq_quark(Quark *quark);
void CORE_stsmlq_sytra1_quark(Quark *quark);
void CORE_stsmlq_corner_quark(Quark *quark);
void CORE_stsmqr_quark(Quark *quark);
void CORE_stsmqr_sytra1_quark(Quark *quark);
void CORE_stsmqr_corner_quark(Quark *quark);
void CORE_stsqrt_quark(Quark *quark);
void CORE_ststrf_quark(Quark *quark);
void CORE_sttmqr_quark(Quark *quark);
void CORE_sttqrt_quark(Quark *quark);
void CORE_sttmlq_quark(Quark *quark);
void CORE_sttlqt_quark(Quark *quark);
void CORE_sormlq_quark(Quark *quark);
void CORE_sormqr_quark(Quark *quark);

void CORE_slaswp_quark(Quark* quark);
void CORE_slaswp_f2_quark(Quark* quark);
void CORE_slaswp_ontile_quark(Quark *quark);
void CORE_slaswp_ontile_f2_quark(Quark *quark);
void CORE_strmm_p2_quark(Quark* quark);
void CORE_sgemm_f2_quark(Quark* quark);
void CORE_sgemm_p2_quark(Quark* quark);
void CORE_sgemm_p2f1_quark(Quark* quark);
void CORE_sgemm_p3_quark(Quark* quark);





void CORE_strdalg_v2_quark(Quark *quark);


void QUARK_CORE_strdalg_v2(Quark *quark, Quark_Task_Flags *task_flags,
                        int uplo,
                        PLASMA_desc *A,
                        float *C,
                        float *S,
                        int grsiz, int s, int id, int blktile);

void CORE_strdalg_v2(PLASMA_enum uplo,  
                  PLASMA_desc *pA, float *C, float *S,
                  int grsiz, int s, int id, int blktile);


#ifdef __cplusplus
}
#endif

#endif
