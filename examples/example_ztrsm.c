/**
 *
 * @file example_ztrsm.c
 *
 *  PLASMA testing routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @brief Example for solving a system of linear equations using 
 *        Cholesky factorization and PLASMA_ztrsm routine
 *
 * @version 2.4.2
 * @author Bilel Hadri
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

int check_solution(int, int, PLASMA_Complex64_t*, int, PLASMA_Complex64_t*, PLASMA_Complex64_t*, int);

int IONE=1;
int ISEED[4] = {0,0,0,1};   /* initial seed for zlarnv() */

int main ()
{

    int cores = 2;
    int N     = 10 ;
    int LDA   = 10 ;
    int NRHS  = 5 ;
    int LDB   = 10 ;
    int info;
    int info_solution;
    int i,j;
    int NminusOne = N-1;
    int LDBxNRHS = LDB*NRHS;

    PLASMA_Complex64_t *A1   = (PLASMA_Complex64_t *)malloc(LDA*N*sizeof(PLASMA_Complex64_t));
    PLASMA_Complex64_t *A2   = (PLASMA_Complex64_t *)malloc(LDA*N*sizeof(PLASMA_Complex64_t));
    PLASMA_Complex64_t *B1   = (PLASMA_Complex64_t *)malloc(LDB*NRHS*sizeof(PLASMA_Complex64_t));
    PLASMA_Complex64_t *B2   = (PLASMA_Complex64_t *)malloc(LDB*NRHS*sizeof(PLASMA_Complex64_t));
    PLASMA_Complex64_t *WORK = (PLASMA_Complex64_t *)malloc(2*LDA*sizeof(PLASMA_Complex64_t));
    double *D                = (double *)malloc(LDA*sizeof(double));

    /* Check if unable to allocate memory */
    if ((!A1)||(!A2)||(!B1)||(!B2)){
        printf("Out of Memory \n ");
        exit(0);
    }

    /* Plasma Initialize */
    PLASMA_Init(cores);
    printf("-- PLASMA is initialized to run on %d cores. \n",cores);

    /* Initialize A1 and A2 for Symmetric Positive Matrix */
    LAPACKE_dlarnv_work(IONE, ISEED, LDA, D);
    zlaghe(&N, &NminusOne, D, A1, &LDA, ISEED, WORK, &info);
    for ( i = 0; i < N; i++)
        for (  j = 0; j < N; j++)
            A2[LDA*j+i] = A1[LDA*j+i];

    for ( i = 0; i < N; i++){
      A1[LDA*i+i] = A1[LDA*i+i]+ (PLASMA_Complex64_t)N ;
       A2[LDA*i+i] = A1[LDA*i+i];
    }

    /* Initialize B1 and B2 */
    LAPACKE_zlarnv_work(IONE, ISEED, LDBxNRHS, B1);
    for ( i = 0; i < N; i++)
        for ( j = 0; j < NRHS; j++)
            B2[LDB*j+i] = B1[LDB*j+i];

    /* PLASMA routines */
    info = PLASMA_zpotrf(PlasmaLower, N, A2, LDA);
    info = PLASMA_ztrsm(PlasmaLeft, PlasmaLower, PlasmaNoTrans, PlasmaNonUnit, 
                        N, NRHS, (PLASMA_Complex64_t)1.0, A2, LDA, B2, LDB);
    info = PLASMA_ztrsm(PlasmaLeft, PlasmaLower, PlasmaConjTrans, PlasmaNonUnit,
                        N, NRHS, (PLASMA_Complex64_t)1.0, A2, LDA, B2, LDB);

    /* Check the solution */
    info_solution = check_solution(N, NRHS, A1, LDA, B1, B2, LDB);

    if ((info_solution != 0)|(info != 0))
       printf("-- Error in ZTRSM example ! \n");
    else
       printf("-- Run of ZTRSM example successful ! \n");

    free(A1); free(A2); free(B1); free(B2); free(WORK); free(D);

    PLASMA_Finalize();

    exit(0);
}


/*------------------------------------------------------------------------
 *  Check the accuracy of the solution of the linear system
 */

int check_solution(int N, int NRHS, PLASMA_Complex64_t *A1, int LDA, PLASMA_Complex64_t *B1, PLASMA_Complex64_t *B2, int LDB)
{
    int info_solution;
    double Rnorm, Anorm, Xnorm, Bnorm;
    PLASMA_Complex64_t alpha, beta;
    double *work = (double *)malloc(N*sizeof(double));
    double eps;

    eps = LAPACKE_dlamch_work('e');

    alpha = 1.0;
    beta  = -1.0;

    Xnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, NRHS, B2, LDB, work);
    Anorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, N, A1, LDA, work);
    Bnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, NRHS, B1, LDB, work);

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, NRHS, N, CBLAS_SADDR(alpha), A1, LDA, B2, LDB, CBLAS_SADDR(beta), B1, LDB);
    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, NRHS, B1, LDB, work);

    printf("============\n");
    printf("Checking the Residual of the solution \n");
    printf("-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps) = %e \n",Rnorm/((Anorm*Xnorm+Bnorm)*N*eps));

    if (Rnorm/((Anorm*Xnorm+Bnorm)*N*eps) > 10.0){
        printf("-- The solution is suspicious ! \n");
        info_solution = 1;
     }
    else{
        printf("-- The solution is CORRECT ! \n");
        info_solution = 0;
    }

    free(work);

    return info_solution;
}
