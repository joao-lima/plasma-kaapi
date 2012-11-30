/**
 *
 * @file core_blas.h
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @date 2010-11-15
 *
 **/
#ifndef _PLASMA_CORE_BLAS_H_
#define _PLASMA_CORE_BLAS_H_

#include <cblas.h>
#include "quark.h"

#include "core_zblas.h"
#include "core_dblas.h"
#include "core_cblas.h"
#include "core_sblas.h"
#include "core_zcblas.h"
#include "core_dsblas.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* 
   * Coreblas Error 
   */
#define coreblas_error(k, str) fprintf(stderr, "%s: Parameter %d / %s\n", __func__, k, str)

  /* 
   * CBlas enum 
   */
#define CBLAS_TRANSPOSE enum CBLAS_TRANSPOSE
#define CBLAS_UPLO      enum CBLAS_UPLO
#define CBLAS_DIAG      enum CBLAS_DIAG
#define CBLAS_SIDE      enum CBLAS_SIDE

/** ****************************************************************************
 *  LAPACK Constants
 **/
extern char *plasma_lapack_constants[];
#define lapack_const(plasma_const) plasma_lapack_constants[plasma_const][0]
  
  /* 
   * Functions which don't depend on precision
   */
void CORE_free_quark(Quark *quark);
void CORE_foo_quark(Quark *quark);
void CORE_foo2_quark(Quark *quark);

void QUARK_CORE_free(Quark *quark, Quark_Task_Flags *task_flags,
                     void *A, int szeA);

#ifdef __cplusplus
}
#endif

#endif
