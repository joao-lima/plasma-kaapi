/**
 *
 * @file core_dsblas.h
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated ds Thu Sep 15 12:08:49 2011
 *
 **/
#ifndef _PLASMA_CORE_DSBLAS_H_
#define _PLASMA_CORE_DSBLAS_H_
#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

/** ****************************************************************************
 *  Declarations of serial kernels - alphabetical order
 **/
void CORE_slag2d(int m, int n, 
                 float *A, int lda, 
                 double *B, int ldb);
void CORE_dlag2s(int m, int n, 
                 double *A, int lda, 
                 float *B, int ldb, int *info);

/** ****************************************************************************
 *  Declarations of QUARK wrappers (called by PLASMA) - alphabetical order
 **/
void QUARK_CORE_slag2d(Quark *quark, Quark_Task_Flags *task_flags,
                      int m, int n, int nb,
                      float *A, int lda,
                       double *B, int ldb);
void QUARK_CORE_dlag2s(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int nb,
                       double *A, int lda,
                       float *B, int ldb,
                       PLASMA_sequence *sequence, PLASMA_request *request);

/** ****************************************************************************
 *  Declarations of QUARK wrappers (called by QUARK) - alphabetical order
 **/
void CORE_slag2d_quark(Quark *quark);
void CORE_dlag2s_quark(Quark *quark);

#ifdef __cplusplus
}
#endif

#endif
