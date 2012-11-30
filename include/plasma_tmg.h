/**
 *
 * @file lapack.h
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @date 2010-11-15
 *
 **/
#ifndef _PLASMA_LAPACK_H_
#define _PLASMA_LAPACK_H_

#if defined( _WIN32 ) || defined( _WIN64 )
#include <float.h>
#define isnan _isnan
#endif

#ifdef ADD_
    #define slagsy slagsy_
    #define dlagsy dlagsy_
    #define clagsy clagsy_
    #define zlagsy zlagsy_
    #define claghe claghe_
    #define zlaghe zlaghe_
#endif

#ifdef UPCASE
    #define slagsy SLAGSY
    #define dlagsy DLAGSY
    #define clagsy CLAGSY
    #define zlagsy ZLAGSY
    #define claghe CLAGHE
    #define zlaghe ZLAGHE
#endif

#ifdef __cplusplus
extern "C" {
#endif

#if 0
void slagsy(int*, int*, float*,  float*,               int*, int*, float*,               int*);
void dlagsy(int*, int*, double*, double*,              int*, int*, double*,              int*);
void clagsy(int*, int*, float*,  PLASMA_Complex32_t *, int*, int*, PLASMA_Complex32_t *, int*);
void zlagsy(int*, int*, double*, PLASMA_Complex64_t *, int*, int*, PLASMA_Complex64_t *, int*);
void claghe(int*, int*, float*,  PLASMA_Complex32_t *, int*, int*, PLASMA_Complex32_t *, int*);
void zlaghe(int*, int*, double*, PLASMA_Complex64_t *, int*, int*, PLASMA_Complex64_t *, int*);
#endif

#ifdef __cplusplus
}
#endif

#endif
