/**
 *
 * @file testing_dgemm.c
 *
 *  PLASMA testing routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Emmanuel Agullo
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Thu Sep 15 12:09:32 2011
 *
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <plasma.h>
#include <cblas.h>
#include <lapacke.h>
#include <core_blas.h>
#include "testing_dmain.h"

#undef COMPLEX
#define REAL

int   norm[4]    = { PlasmaMaxNorm, PlasmaOneNorm, PlasmaInfNorm, PlasmaFrobeniusNorm };
char *normstr[4] = { "Max", "One", "Inf", "Fro" };


int testing_dlange(int argc, char **argv)
{
    /* Check for number of arguments*/
    if ( argc != 3) {
        USAGE("LANGE", "M N LDA",
              "   - M      : number of rows of matrices A and C\n"
              "   - N      : number of columns of matrices B and C\n"
              "   - LDA    : leading dimension of matrix A\n");
        return -1;
    }

    int M     = atoi(argv[0]);
    int N     = atoi(argv[1]);
    int LDA   = atoi(argv[2]);
    int LDAxN = LDA*N;
    int n, u;
    double eps;

    double *A    = (double *)malloc(LDAxN*sizeof(double));
    double             *work = (double*) malloc(max(M,N)*sizeof(double));
    double normplasma, normlapack;

    eps = LAPACKE_dlamch_work('e');

    printf("\n");
    printf("------ TESTS FOR PLASMA DLANGE ROUTINE -------  \n");
    printf("            Size of the Matrix %d by %d\n", M, N);
    printf("\n");
    printf(" The matrix A is randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n",eps);
    printf(" Computational tests pass if scaled residuals are less than 10.\n");

    /*----------------------------------------------------------
     *  TESTING DLANGE
     */

    /* Initialize A, B, C */
    LAPACKE_dlarnv_work(IONE, ISEED, LDAxN, A);

    /* PLASMA DLANGE */
    for(n=0; n<3; n++) {
        normplasma = PLASMA_dlange(norm[n], M, N, A, LDA, work);
        normlapack = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(norm[n]), M, N, A, LDA, work);
        
        printf("Lapack %e, Plasma %e\n", normlapack, normplasma);
        
        printf("***************************************************\n");
        if ( abs(normlapack-normplasma) < eps ) {
            printf(" ---- TESTING DLANGE (%s)............... PASSED !\n", normstr[n]);
        }
        else {
            printf(" - TESTING DLANGE (%s)... FAILED !\n", normstr[n]);
        }
        printf("***************************************************\n");
    }      

    /* PLASMA DLANSY */
    for(n=0; n<3; n++) {
        for(u=0; u<2; u++) {
            normplasma = PLASMA_dlansy(norm[n], uplo[u], min(M,N), A, LDA, work);
            normlapack = LAPACKE_dlansy_work(LAPACK_COL_MAJOR, lapack_const(norm[n]), lapack_const(uplo[u]), min(M,N), A, LDA, work);
        
            printf("Lapack %e, Plasma %e\n", normlapack, normplasma);
            
            printf("***************************************************\n");
            if ( abs(normlapack-normplasma) < eps ) {
                printf(" ---- TESTING DLANSY (%s, %s)......... PASSED !\n", normstr[n], uplostr[u]);
            }
            else {
                printf(" - TESTING DLANSY (%s, %s)... FAILED !\n", normstr[n], uplostr[u]);
            }
            printf("***************************************************\n");
        }      
    }

#ifdef COMPLEX
    /* PLASMA DLANSY */
    {
      int j;
      for (j=0; j<min(M,N); j++) {
        A[j*LDA+j] -= I*cimag(A[j*LDA+j]);
      }
    }

    for(n=0; n<3; n++) {
        for(u=0; u<2; u++) {
            normplasma = PLASMA_dlansy(norm[n], uplo[u], min(M,N), A, LDA, work);
            normlapack = LAPACKE_dlansy_work(LAPACK_COL_MAJOR, lapack_const(norm[n]), lapack_const(uplo[u]), min(M,N), A, LDA, work);
        
            printf("Lapack %e, Plasma %e\n", normlapack, normplasma);
            
            printf("***************************************************\n");
            if ( abs(normlapack-normplasma) < eps ) {
                printf(" ---- TESTING DLANSY (%s, %s)......... PASSED !\n", normstr[n], uplostr[u]);
            }
            else {
                printf(" - TESTING DLANSY (%s, %s)... FAILED !\n", normstr[n], uplostr[u]);
            }
            printf("***************************************************\n");
        }      
    }
#endif

    free(A); free(work);
    return 0;
}
