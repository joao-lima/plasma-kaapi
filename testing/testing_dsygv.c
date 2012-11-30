/**
 *
 * @file testing_dsygv.c
 *
 *  PLASMA testing routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @generated d Thu Sep 15 12:09:31 2011
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
#include "testing_dmain.h"

#undef COMPLEX
#define REAL

static int check_orthogonality(int, int, double*, int, double);
static int check_reduction(int, int, int, int, double*, double*, int, double*, int, double*, double);
#if 0
static int check_solution(int, int, int, double*, int, double*, double*, int, double);
#endif
static int check_solution(int, double*, double*, double);

int testing_dsygv(int argc, char **argv)
{
    /* Check for number of arguments*/
    if (argc != 3) {
        USAGE("HEGV", "N LDA LDB",
              "   - N    : size of the matrices A and B\n"
              "   - LDA  : leading dimension of the matrix A\n"
              "   - LDB  : leading dimension of the matrix B\n");
        return -1;
    }

    double      eps = LAPACKE_dlamch_work('e');
    PLASMA_enum vec = PlasmaNoVec;
    int    N        = atoi(argv[0]);
    int    LDA      = atoi(argv[1]);
    int    LDB      = atoi(argv[2]);
    int    LDQ      = LDA;
    int    LDAxN    = LDA*N;
    int    LDBxN    = LDB*N;
    int    LDQxN    = LDQ*N;

    int info_ortho     = 0;
    int info_solution  = 0;
    int info_reduction = 0;
    int i, u;

    double *A1    = (double *)malloc(LDAxN*sizeof(double));
    double *A2    = (double *)malloc(LDAxN*sizeof(double));
    double *B1    = (double *)malloc(LDBxN*sizeof(double));
    double *B2    = (double *)malloc(LDBxN*sizeof(double));
    double *Q     = (double *)malloc(LDQxN*sizeof(double));
    double *Ainit = (double *)malloc(LDAxN*sizeof(double));
    double *Binit = (double *)malloc(LDBxN*sizeof(double));
    double *W1 = (double *)malloc(N*sizeof(double));
    double *W2 = (double *)malloc(N*sizeof(double));
    double *work = (double *)malloc(3*N* sizeof(double));
    PLASMA_desc *T;

    /* Check if unable to allocate memory */
    if ((!A1)||(!A2)||(!B1)||(!B2)||(!Q)||(!Ainit)||(!Binit)){
        printf("Out of Memory \n ");
        return -2;
    }

    /*
    PLASMA_Disable(PLASMA_AUTOTUNING);
    PLASMA_Set(PLASMA_TILE_SIZE, 120);
    PLASMA_Set(PLASMA_INNER_BLOCK_SIZE, 20);
    */
    
    PLASMA_Enable(PLASMA_WARNINGS);
    PLASMA_Enable(PLASMA_ERRORS);
    PLASMA_Alloc_Workspace_dsygv(N, N, &T);

    /*----------------------------------------------------------
    *  TESTING DSYGV
    */

    /* Initialize A1 and Ainit */
    PLASMA_dplgsy(0., N, A1, LDA, 5198);
    LAPACKE_dlacpy_work(LAPACK_COL_MAJOR, 'A', N, N, A1, LDA, Ainit, LDA);

    /* Initialize B1 and Binit */
    PLASMA_dplgsy((double)N, N, B1, LDB, 4321 );
    LAPACKE_dlacpy_work(LAPACK_COL_MAJOR, 'A', N, N, B1, LDB, Binit, LDB);

    printf("\n");
    printf("------ TESTS FOR PLASMA DSYGV ROUTINE -------  \n");
    printf("        Size of the Matrix %d by %d\n", N, N);
    printf("\n");
    printf(" The matrix A is randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n",eps);
    printf(" Computational tests pass if scaled residuals are less than 60.\n");

    /*----------------------------------------------------------
     *  TESTING DSYGV
     */

    for (i=0; i<3; i++) {
        for (u=0; u<2; u++) {
            LAPACKE_dlaset_work(LAPACK_COL_MAJOR, 'A', LDA, N, 0., 1., Q, LDA);

            memcpy(A2, Ainit, LDAxN*sizeof(double));
            memcpy(B2, Binit, LDBxN*sizeof(double));

            PLASMA_dsygv(itype[i], vec, uplo[u], N, A2, LDA, B2, LDB, W2, T, Q, LDQ);

            /* Check the orthogonality, reduction and the eigen solutions */
            if (vec == PlasmaVec)
                info_ortho = check_orthogonality(N, N, Q, LDA, eps);
            /* 
             * WARNING: For now, Q is associated to Band tridiagonal reduction and 
             * not to the final tridiagonal reduction, so we can not call the check
             */
            if (0)
                info_reduction = check_reduction(itype[i], uplo[u], N, 1, A1, A2, LDA, B2, LDB, Q, eps);

            memcpy(A1, Ainit, LDAxN*sizeof(double));
            memcpy(B1, Binit, LDBxN*sizeof(double));

            LAPACKE_dsygv( LAPACK_COL_MAJOR, 
                     itype[i], lapack_const(vec), lapack_const(uplo[u]),
                     N, A1, LDA, B1, LDB, W1);

            /*info_solution  = check_solution(N, N, N, A1, LDA, B1, B2, LDB, eps);*/
            info_solution = check_solution(N, W1, W2, eps);
         
            if ( (info_ortho == 0) & (info_reduction == 0) & (info_solution == 0)) {
                printf("***************************************************\n");
                printf(" ---- TESTING DSYGV (%s, %s) ...................... PASSED !\n", itypestr[i], uplostr[u]);
                printf("***************************************************\n");
            }
            else {
                printf("************************************************\n");
                printf(" - TESTING DSYGV (%s, %s) ... FAILED !\n", itypestr[i], uplostr[u]);
                printf("************************************************\n");
            } 
        }
    }

    PLASMA_Dealloc_Handle_Tile(&T);
    free(A1); 
    free(A2); 
    free(B1); 
    free(B2); 
    free(Q); 
    free(Ainit); 
    free(Binit); 
    free(W1);
    free(W2);
    free(work);

    return 0;
}

/*-------------------------------------------------------------------
 * Check the orthogonality of Q
 */
static int check_orthogonality(int M, int N, double *Q, int LDQ, double eps)
{
    double  alpha =  1.0;
    double  beta  = -1.0;
    double  normQ, result;
    int     info_ortho;
    int     minMN = min(M, N);
    double *work = (double *)malloc(minMN*sizeof(double));

    /* Build the idendity matrix */
    double *Id = (double *) malloc(minMN*minMN*sizeof(double));
    LAPACKE_dlaset_work(LAPACK_COL_MAJOR, 'A', minMN, minMN, 0., 1., Id, minMN);

    /* Perform Id - Q'Q */
    if (M >= N)
        cblas_dsyrk(CblasColMajor, CblasUpper, CblasTrans, N, M, alpha, Q, LDQ, beta, Id, N);
    else
        cblas_dsyrk(CblasColMajor, CblasUpper, CblasNoTrans,   M, N, alpha, Q, LDQ, beta, Id, M);

    normQ = LAPACKE_dlansy_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), 'U', minMN, Id, minMN, work);

    result = normQ / (minMN * eps);
    printf("============\n");
    printf("Checking the orthogonality of Q \n");
    printf("||Id-Q'*Q||_oo / (minMN*eps) = %e \n", result);

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        printf("-- Orthogonality is suspicious ! \n");
        info_ortho=1;
    }
    else {
        printf("-- Orthogonality is CORRECT ! \n");
        info_ortho=0;
    }

    free(work); free(Id);

    return info_ortho;
}

/*------------------------------------------------------------
 *  Check the tridiagonal reduction
 */

static int check_reduction(int itype, int uplo, int N, int bw, 
                           double *A1, double *A2, int LDA, 
                           double *B2, int LDB, double *Q, double eps )
{
    double alpha = 1.0;
    double beta  = 0.0;
    double Anorm, Rnorm, result;
    int info_reduction;
    int i, j;
    char *str;
    
    double *Aorig    = (double *)malloc(N*N*sizeof(double));
    double *Residual = (double *)malloc(N*N*sizeof(double));
    double *T        = (double *)malloc(N*N*sizeof(double));
    double *work                 = (double *)malloc(N*sizeof(double));

    memset((void*)T, 0, N*N*sizeof(double));
    LAPACKE_dlacpy_work(LAPACK_COL_MAJOR, ' ', N, N, A1, LDA, Residual, N);

    /* Rebuild the T */
    LAPACKE_dlacpy_work(LAPACK_COL_MAJOR, lapack_const(uplo), N, N, A2, LDA, T, N);

    if (uplo == PlasmaLower) {
        /* Set the reflectors to 0 */
        for (i = bw+1; i < N; i++)
            for (j = 0 ; (j < N) && (j < i-bw); j++)
                T[j*N+i] = 0.;

        /* Copy the lower part to the upper part to rebuild the symmetry */
        for (i = 0; i < N; i++)
            for (j = 0 ; j < i; j++)
                T[i*N+j] = (T[j*N+i]);
    } else {
        /* Set the reflectors to 0 */
        for (j = bw+1; j < N; j++)
            for (i = 0 ; (i < N) && (i < j-bw); i++)
                T[j*N+i] = 0.;
        /* Copy the upper part to the lower part to rebuild the symmetry */
        for (i = 0; i < N; i++)
            for (j = i+1 ; j < N; j++)
                T[i*N+j] = (T[j*N+i]);
    }

    memset((void*)Aorig, 0, N*N*sizeof(double));

    if (itype == 1) {
        if (uplo == PlasmaLower) {
            str = "L*Q*T*Q'*L'";
   
            /* Compute Aorig=Q*T*Q' */
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,   N, N, N, (alpha), Q,     LDA, T, N,   (beta), Aorig, N);
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, N, N, N, (alpha), Aorig, N,   Q, LDA, (beta), T,     N);

            LAPACKE_dlacpy_work(LAPACK_COL_MAJOR, PlasmaUpperLower, N, N, T, N, Aorig, N);

            cblas_dtrmm(CblasColMajor, CblasLeft,  CblasLower, CblasNoTrans,   CblasNonUnit, N, N, (alpha), B2, LDB, Aorig, N);   
            cblas_dtrmm(CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit, N, N, (alpha), B2, LDB, Aorig, N); 

        }
        else {
            str = "U'*Q*T*Q'*U";

            /* Compute Aorig=Q'*T*Q */
            cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, N, N, N, (alpha), Q,     LDA, T, N,   (beta), Aorig, N);
            cblas_dgemm(CblasColMajor, CblasNoTrans,   CblasNoTrans, N, N, N, (alpha), Aorig, N,   Q, LDA, (beta), T,     N);

            LAPACKE_dlacpy_work(LAPACK_COL_MAJOR, 'A', N, N, T, N, Aorig, N);

            cblas_dtrmm(CblasColMajor, CblasLeft,  CblasUpper, CblasTrans, CblasNonUnit, N, N, (alpha), B2, LDB, Aorig, N);   
            cblas_dtrmm(CblasColMajor, CblasRight, CblasUpper, CblasNoTrans,   CblasNonUnit, N, N, (alpha), B2, LDB, Aorig, N);  
        }
    }
    else {
        if (uplo == PlasmaLower) {
            str = "inv(L')*Q*A2*Q'*inv(L)";
            
            /* Compute Aorig=Q*T*Q' */
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,   N, N, N, (alpha), Q,     LDA, T, N,   (beta), Aorig, N);
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, N, N, N, (alpha), Aorig, N,   Q, LDA, (beta), T, N   );

            LAPACKE_dlacpy_work(LAPACK_COL_MAJOR, 'A', N, N, T, N, Aorig, N);

            cblas_dtrsm(CblasColMajor, CblasLeft,  CblasLower, CblasTrans, CblasNonUnit, N, N, (alpha), B2, LDB, Aorig, N);   
            cblas_dtrsm(CblasColMajor, CblasRight, CblasLower, CblasNoTrans,   CblasNonUnit, N, N, (alpha), B2, LDB, Aorig, N);   
        }
        else{
            str = "inv(U)*Q*A2*Q'*inv(U')";

            /* Compute Aorig=Q'*T*Q */
            cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, N, N, N, (alpha), Q,     LDA, T, N,   (beta), Aorig, N);
            cblas_dgemm(CblasColMajor, CblasNoTrans,   CblasNoTrans, N, N, N, (alpha), Aorig, N,   Q, LDA, (beta), T,     N);
      
            LAPACKE_dlacpy_work(LAPACK_COL_MAJOR, 'A', N, N, T, N, Aorig, N);
          
            cblas_dtrsm(CblasColMajor, CblasLeft,  CblasUpper, CblasNoTrans,   CblasNonUnit, N, N, (alpha), B2, LDB, Aorig, N);   
            cblas_dtrsm(CblasColMajor, CblasRight, CblasUpper, CblasTrans, CblasNonUnit, N, N, (alpha), B2, LDB, Aorig, N);   
        }
    }
   
    /* Compute the Residual */
    for (i = 0; i < N; i++)
        for (j = 0 ; j < N; j++)
            Residual[j*N+i] = A1[j*LDA+i]-Aorig[j*N+i];
    
    Rnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, N, Residual, N,   work);
    Anorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, N, A2,       LDA, work);
    
    result = Rnorm / (Anorm * N *eps);
    printf("============\n");
    printf("Checking the tridiagonal reduction \n");
    printf("-- ||A-%s||_oo/(||A||_oo.N.eps) = %e \n", str, result );
    
    if (isnan(result) || isinf(result) || (result > 60.0) ) {
        printf("-- Reduction is suspicious ! \n");
        info_reduction = 1;
    }
    else {
        printf("-- Reduction is CORRECT ! \n");
        info_reduction = 0;
    }
    
    free(Aorig); free(Residual); free(T); free(work);
    
    return info_reduction;
}

/*--------------------------------------------------------------
 * Check the solution
 */
static int check_solution(int N, double *E1, double *E2, double eps)
{
    int info_solution, i;
    double resid;
    double maxtmp;
    double maxel = fabs( fabs(E1[0]) - fabs(E2[0]) );
    double maxeig = max( fabs(E1[0]), fabs(E2[0]) );
    for (i = 1; i < N; i++){
        resid   = fabs(fabs(E1[i])-fabs(E2[i]));
        maxtmp  = max(fabs(E1[i]), fabs(E2[i]));

        /* Update */
        maxeig = max(maxtmp, maxeig);
        maxel  = max(resid,  maxel );
    }

    maxel = maxel / (maxeig * N * eps);
    printf(" ======================================================\n");
    printf(" | D - eigcomputed | / (|D| * N * eps) : %15.3E \n",  maxel );
    printf(" ======================================================\n");

    printf("============\n");
    printf("Checking the eigenvalues of A\n");
    if ( isnan(maxel) || isinf(maxel) || (maxel > 100) ) {
        printf("-- The eigenvalues are suspicious ! \n");
        info_solution = 1;
    }
    else{
        printf("-- The eigenvalues are CORRECT ! \n");
        info_solution = 0;
    }

    return info_solution;
}

#if 0
static int check_solution(int M, int N, int NRHS, double *A1, int LDA, 
                          double *B1, double *B2, int LDB, double eps)
{
    double alpha = 1.0;
    double beta  = 0.0;
    double *Residual;
    int info_solution;
    int maxMN = max(M, N);
    double Rnorm, Anorm, Xnorm, Bnorm, result;
    double *work = (double *)malloc( max(maxMN, NRHS)* sizeof(double));

    Anorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, N,    A1, LDA, work);
    Xnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, NRHS, B2, LDB, work);
    Bnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, NRHS, B1, LDB, work);

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, NRHS, N, 
                (alpha), A1, LDA, 
                                    B2, LDB, 
                (beta),  B1, LDB);

    Residual = (double *)malloc(maxMN*NRHS*sizeof(double));
    memset((void*)Residual, 0, maxMN*NRHS*sizeof(double));
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, N, NRHS, M, 
                (alpha), A1, LDA, 
                                    B1, LDB, 
                (beta),  Residual, maxMN);
    
    Rnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), maxMN, NRHS, Residual, maxMN, work);
    free(Residual);
    free(work);

    if (getenv("PLASMA_TESTING_VERBOSE"))
       printf( "||A||_oo=%f\n||X||_oo=%f\n||B||_oo=%f\n||A X - B||_oo=%e\n", Anorm, Xnorm, Bnorm, Rnorm );

    result = Rnorm / ( (Anorm*Xnorm+Bnorm)*N*eps ) ;
    printf("============\n");
    printf("Checking the Residual of the solution \n");
    printf("-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps) = %e \n", result);
    
    if (  isnan(Xnorm) || isinf(Xnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        printf("-- The solution is suspicious ! \n");
        info_solution = 1;
    }
    else{
        printf("-- The solution is CORRECT ! \n");
        info_solution = 0;
    }
    return info_solution;
}
#endif
