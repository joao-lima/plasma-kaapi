/**
 *
 * @file testing_zsyrk.c
 *
 *  PLASMA testing routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <plasma.h>
#include <cblas.h>
#include <lapacke.h>
#include <plasma_tmg.h>
#include <core_blas.h>
#include "testing_zmain.h"

static int check_solution(PLASMA_enum uplo, PLASMA_enum trans, int N, int K,
                          PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA,
                          PLASMA_Complex64_t beta,  PLASMA_Complex64_t *Cref, PLASMA_Complex64_t *Cplasma, int LDC);


int testing_zsyrk(int argc, char **argv)
{
    /* Check for number of arguments*/
    if ( argc != 6){
        USAGE("SYRK", "alpha beta M N LDA LDC",
              "   - alpha : alpha coefficient\n"
              "   - beta : beta coefficient\n"
              "   - N : number of columns and rows of matrix C and number of row of matrix A\n"
              "   - K : number of columns of matrix A\n"
              "   - LDA : leading dimension of matrix A\n"
              "   - LDC : leading dimension of matrix C\n");
        return -1;
    }

    PLASMA_Complex64_t alpha = (PLASMA_Complex64_t) atol(argv[0]);
    PLASMA_Complex64_t beta  = (PLASMA_Complex64_t) atol(argv[1]);
    int N     = atoi(argv[2]);
    int K     = atoi(argv[3]);
    int LDA   = atoi(argv[4]);
    int LDC   = atoi(argv[5]);
    int NKmax = max(N, K);
    int NminusOne = N - 1;

    double eps;
    int info_solution;
    int info, u, t;
    size_t LDAxK = LDA*NKmax;
    size_t LDCxN = LDC*N;

    PLASMA_Complex64_t *A      = (PLASMA_Complex64_t *)malloc(LDAxK*sizeof(PLASMA_Complex64_t));
    PLASMA_Complex64_t *C      = (PLASMA_Complex64_t *)malloc(LDCxN*sizeof(PLASMA_Complex64_t));
    PLASMA_Complex64_t *Cinit  = (PLASMA_Complex64_t *)malloc(LDCxN*sizeof(PLASMA_Complex64_t));
    PLASMA_Complex64_t *Cfinal = (PLASMA_Complex64_t *)malloc(LDCxN*sizeof(PLASMA_Complex64_t));
    PLASMA_Complex64_t *WORK   = (PLASMA_Complex64_t *)malloc(2*LDC*sizeof(PLASMA_Complex64_t));
    double             *D      = (double *)            malloc(LDC  *sizeof(double));

    /* Check if unable to allocate memory */
    if ( (!A) || (!Cinit) || (!Cfinal) || (!D) ){
        printf("Out of Memory \n ");
        return -2;
    }

    eps = LAPACKE_dlamch_work('e');

    printf("\n");
    printf("------ TESTS FOR PLASMA ZSYRK ROUTINE -------  \n");
    printf("            Size of the Matrix A %d by %d\n", N, K);
    printf("\n");
    printf(" The matrix A is randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n",eps);
    printf(" Computational tests pass if scaled residuals are less than 10.\n");

    /*----------------------------------------------------------
    *  TESTING ZSYRK
    */

    /* Initialize A */
    LAPACKE_zlarnv_work(IONE, ISEED, LDAxK, A);

    /* Initialize C */
    LAPACKE_dlarnv_work(IONE, ISEED, LDC, D);
    zlagsy(&N, &NminusOne, D, C, &LDC, ISEED, WORK, &info);
    free(D); free(WORK);

    for (u=0; u<2; u++) {
        for (t=0; t<2; t++) {
            memcpy(Cinit,  C, LDCxN*sizeof(PLASMA_Complex64_t));
            memcpy(Cfinal, C, LDCxN*sizeof(PLASMA_Complex64_t));

            /* PLASMA ZSYRK */
            PLASMA_zsyrk(uplo[u], trans[t], N, K, alpha, A, LDA, beta, Cfinal, LDC);

            /* Check the solution */
            info_solution = check_solution(uplo[u], trans[t], N, K,
                                           alpha, A, LDA, beta, Cinit, Cfinal, LDC);

            if (info_solution == 0) {
                printf("***************************************************\n");
                printf(" ---- TESTING ZSYRK (%5s, %s) ........... PASSED !\n", uplostr[u], transstr[t]);
                printf("***************************************************\n");
            }
            else {
                printf("************************************************\n");
                printf(" - TESTING ZSYRK (%5s, %s) ... FAILED !\n", uplostr[u], transstr[t]);
                printf("************************************************\n");
            }
        }
    }

    free(A); free(C);
    free(Cinit); free(Cfinal);

    return 0;
}

/*--------------------------------------------------------------
 * Check the solution
 */

static int check_solution(PLASMA_enum uplo, PLASMA_enum trans, int N, int K,
                          PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA,
                          PLASMA_Complex64_t beta,  PLASMA_Complex64_t *Cref, PLASMA_Complex64_t *Cplasma, int LDC)
{
    int info_solution;
    double Anorm, Cinitnorm, Cplasmanorm, Clapacknorm, Rnorm;
    double eps;
    PLASMA_Complex64_t beta_const;
    double result;
    double *work = (double *)malloc(max(N, K)* sizeof(double));

    beta_const  = -1.0;
    Anorm       = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm),
                                (trans == PlasmaNoTrans) ? N : K,
                                (trans == PlasmaNoTrans) ? K : N, A, LDA, work);
    Cinitnorm   = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, N, Cref,    LDC, work);
    Cplasmanorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, N, Cplasma, LDC, work);

    cblas_zsyrk(CblasColMajor, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                N, K, CBLAS_SADDR(alpha), A, LDA, CBLAS_SADDR(beta), Cref, LDC);

    Clapacknorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, N, Cref, LDC, work);

    cblas_zaxpy(LDC*N, CBLAS_SADDR(beta_const), Cplasma, 1, Cref, 1);

    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, N, Cref, LDC, work);

    eps = LAPACKE_dlamch_work('e');

    printf("Rnorm %e, Anorm %e, Cinitnorm %e, Cplasmanorm %e, Clapacknorm %e\n",
           Rnorm, Anorm, Cinitnorm, Cplasmanorm, Clapacknorm);

    result = Rnorm / ((Anorm + Cinitnorm) * N * eps);

    printf("============\n");
    printf("Checking the norm of the difference against reference ZSYRK \n");
    printf("-- ||Cplasma - Clapack||_oo/((||A||_oo+||C||_oo).N.eps) = %e \n", result);

    if ( isinf(Clapacknorm) || isinf(Cplasmanorm) || isnan(result) || isinf(result) || (result > 10.0) ) {
        printf("-- The solution is suspicious ! \n");
        info_solution = 1;
    }
    else {
        printf("-- The solution is CORRECT ! \n");
        info_solution= 0 ;
    }

    free(work);

    return info_solution;
}
