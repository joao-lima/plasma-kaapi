/**
 *
 * @file example_dgelqs.c
 *
 *  PLASMA testing routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @brief Example for solving undetermined linear systems
 *
 * @version 2.4.2
 * @author Bilel Hadri
 * @date 2010-11-15
 * @generated d Thu Sep 15 12:09:37 2011
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

#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif

int check_solution(int, int, int, double*, int, double*, double*, int);

int IONE=1;
int ISEED[4] = {0,0,0,1};   /* initial seed for dlarnv() */

int main ()
{

    int cores = 2;
    int M     = 10;
    int N     = 15;
    int LDA   = 10;
    int NRHS  = 5;
    int LDB   = 15;

    int info;
    int info_solution;
    int i,j;
    int LDAxN = LDA*N;
    int LDBxNRHS = LDB*NRHS;

    double *A1 = (double *)malloc(LDA*N*sizeof(double));
    double *A2 = (double *)malloc(LDA*N*sizeof(double));
    double *B1 = (double *)malloc(LDB*NRHS*sizeof(double));
    double *B2 = (double *)malloc(LDB*NRHS*sizeof(double));
    double *T;

    /* Check if unable to allocate memory */
    if ((!A1)||(!A2)||(!B1)||(!B2)){
        printf("Out of Memory \n ");
        exit(0);
    }

    /* Plasma Initialization */
    PLASMA_Init(cores);
    printf("-- PLASMA is initialized to run on %d cores. \n",cores);

    /* Allocate T */
    PLASMA_Alloc_Workspace_dgeqrf(M, N, &T);

    /* Initialize A1 and A2 */
    LAPACKE_dlarnv_work(IONE, ISEED, LDAxN, A1);
    for (i = 0; i < M; i++)
        for (j = 0; j < N; j++)
            A2[LDA*j+i] = A1[LDA*j+i] ;

    /* Initialize B1 and B2 */
    LAPACKE_dlarnv_work(IONE, ISEED, LDBxNRHS, B1);
    for (i = 0; i < M; i++)
        for (j = 0; j < NRHS; j++)
             B2[LDB*j+i] = B1[LDB*j+i] ;

    /* Factorization QR of the matrix A2 */
    info = PLASMA_dgelqf(M, N, A2, LDA, T);

    /* Solve the problem */
    info = PLASMA_dgelqs(M, N, NRHS, A2, LDA, T, B2, LDB);

    /* Check the solution */
    info_solution = check_solution(M, N, NRHS, A1, LDA, B1, B2, LDB);

    if ((info_solution != 0)|(info != 0))
       printf("-- Error in DGELQS example ! \n");
    else
       printf("-- Run of DGELQS example successful ! \n");

    free(A1); free(A2); free(B1); free(B2); free(T);

    PLASMA_Finalize();

    exit(0);
}

/*--------------------------------------------------------------
 * Check the solution
 */

int check_solution(int M, int N, int NRHS, double *A1, int LDA, double *B1, double *B2, int LDB)
{
    int info_solution;
    double Rnorm, Anorm, Xnorm, Bnorm;
    double alpha, beta;
    double *work = (double *)malloc(max(M, N)* sizeof(double));
    double eps;

    eps = LAPACKE_dlamch_work('e');

    alpha = 1.0;
    beta  = -1.0;

    Anorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, N, A1, LDA, work);
    Xnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, NRHS, B2, LDB, work);
    Bnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, NRHS, B1, LDB, work);

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, NRHS, N, (alpha), A1, LDA, B2, LDB, (beta), B1, LDB);

    if (M >= N) {
       double *Residual = (double *)malloc(M*NRHS*sizeof(double));
       memset((void*)Residual, 0, M*NRHS*sizeof(double));
       cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, N, NRHS, M, (alpha), A1, LDA, B1, LDB, (beta), Residual, M);
       Rnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, NRHS, Residual, M, work);
       free(Residual);
    }
    else {
       double *Residual = (double *)malloc(N*NRHS*sizeof(double));
       memset((void*)Residual, 0, N*NRHS*sizeof(double));
       cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, N, NRHS, M, (alpha), A1, LDA, B1, LDB, (beta), Residual, N);
       Rnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, NRHS, Residual, N, work);
       free(Residual);
    }

    printf("============\n");
    printf("Checking the Residual of the solution \n");
    printf("-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||)_oo.N.eps) = %e \n",Rnorm/((Anorm*Xnorm+Bnorm)*N*eps));

    if (isnan(Rnorm / ((Anorm * Xnorm + Bnorm) * N * eps)) || (Rnorm / ((Anorm * Xnorm + Bnorm) * N * eps) > 10.0) ) {
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
