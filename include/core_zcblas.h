/**
 *
 * @file core_zcblas.h
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
 * @precisions mixed zc -> ds
 *
 **/
#ifndef _PLASMA_CORE_ZCBLAS_H_
#define _PLASMA_CORE_ZCBLAS_H_
#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

/** ****************************************************************************
 *  Declarations of serial kernels - alphabetical order
 **/
void CORE_clag2z(int m, int n, 
                 PLASMA_Complex32_t *A, int lda, 
                 PLASMA_Complex64_t *B, int ldb);
void CORE_zlag2c(int m, int n, 
                 PLASMA_Complex64_t *A, int lda, 
                 PLASMA_Complex32_t *B, int ldb, int *info);

/** ****************************************************************************
 *  Declarations of QUARK wrappers (called by PLASMA) - alphabetical order
 **/
void QUARK_CORE_clag2z(Quark *quark, Quark_Task_Flags *task_flags,
                      int m, int n, int nb,
                      PLASMA_Complex32_t *A, int lda,
                       PLASMA_Complex64_t *B, int ldb);
void QUARK_CORE_zlag2c(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int nb,
                       PLASMA_Complex64_t *A, int lda,
                       PLASMA_Complex32_t *B, int ldb,
                       PLASMA_sequence *sequence, PLASMA_request *request);

/** ****************************************************************************
 *  Declarations of QUARK wrappers (called by QUARK) - alphabetical order
 **/
void CORE_clag2z_quark(Quark *quark);
void CORE_zlag2c_quark(Quark *quark);

#ifdef __cplusplus
}
#endif

#endif
