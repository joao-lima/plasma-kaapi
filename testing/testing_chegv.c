/**
 *
 * @file testing_chegv.c
 *
 *  PLASMA testing routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @generated c Thu Sep 15 12:09:31 2011
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
#include "testing_cmain.h"

#undef REAL
#define COMPLEX

static int check_orthogonality(int, int, PLASMA_Complex32_t*, int, float);
static int check_reduction(int, int, int, int, PLASMA_Complex32_t*, PLASMA_Complex32_t*, int, PLASMA_Complex32_t*, int, PLASMA_Complex32_t*, float);
#if 0
static int check_solution(int, int, int, PLASMA_Complex32_t*, int, PLASMA_Complex32_t*, PLASMA_Complex32_t*, int, float);
#endif
static int check_solution(int, float*, float*, float);

int testing_chegv(int argc, char **argv)
{
    /* Check for number of arguments*/
    if (argc != 3) {
        USAGE("HEGV", "N LDA LDB",
              "   - N    : size of the matrices A and B\n"
              "   - LDA  : leading dimension of the matrix A\n"
              "   - LDB  : leading dimension of the matrix B\n");
        return -1;
    }

    float      eps = LAPACKE_slamch_work('e');
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

    PLASMA_Complex32_t *A1    = (PLASMA_Complex32_t *)malloc(LDAxN*sizeof(PLASMA_Complex32_t));
    PLASMA_Complex32_t *A2    = (PLASMA_Complex32_t *)malloc(LDAxN*sizeof(PLASMA_Complex32_t));
    PLASMA_Complex32_t *B1    = (PLASMA_Complex32_t *)malloc(LDBxN*sizeof(PLASMA_Complex32_t));
    PLASMA_Complex32_t *B2    = (PLASMA_Complex32_t *)malloc(LDBxN*sizeof(PLASMA_Complex32_t));
    PLASMA_Complex32_t *Q     = (PLASMA_Complex32_t *)malloc(LDQxN*sizeof(PLASMA_Complex32_t));
    PLASMA_Complex32_t *Ainit = (PLASMA_Complex32_t *)malloc(LDAxN*sizeof(PLASMA_Complex32_t));
    PLASMA_Complex32_t *Binit = (PLASMA_Complex32_t *)malloc(LDBxN*sizeof(PLASMA_Complex32_t));
    float *W1 = (float *)malloc(N*sizeof(float));
    float *W2 = (float *)malloc(N*sizeof(float));
    PLASMA_Complex32_t *work = (PLASMA_Complex32_t *)malloc(3*N* sizeof(PLASMA_Complex32_t));
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
    PLASMA_Alloc_Workspace_chegv(N, N, &T);

    /*----------------------------------------------------------
    *  TESTING CHEGV
    */

    /* Initialize A1 and Ainit */
    PLASMA_cplghe(0., N, A1, LDA, 5198);
    LAPACKE_clacpy_work(LAPACK_COL_MAJOR, 'A', N, N, A1, LDA, Ainit, LDA);

    /* Initialize B1 and Binit */
    PLASMA_cplghe((float)N, N, B1, LDB, 4321 );
    LAPACKE_clacpy_work(LAPACK_COL_MAJOR, 'A', N, N, B1, LDB, Binit, LDB);

    printf("\n");
    printf("------ TESTS FOR PLASMA CHEGV ROUTINE -------  \n");
    printf("        Size of the Matrix %d by %d\n", N, N);
    printf("\n");
    printf(" The matrix A is randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n",eps);
    printf(" Computational tests pass if scaled residuals are less than 60.\n");

    /*----------------------------------------------------------
     *  TESTING CHEGV
     */

    for (i=0; i<3; i++) {
        for (u=0; u<2; u++) {
            LAPACKE_claset_work(LAPACK_COL_MAJOR, 'A', LDA, N, 0., 1., Q, LDA);

            memcpy(A2, Ainit, LDAxN*sizeof(PLASMA_Complex32_t));
            memcpy(B2, Binit, LDBxN*sizeof(PLASMA_Complex32_t));

            PLASMA_chegv(itype[i], vec, uplo[u], N, A2, LDA, B2, LDB, W2, T, Q, LDQ);

            /* Check the orthogonality, reduction and the eigen solutions */
            if (vec == PlasmaVec)
                info_ortho = check_orthogonality(N, N, Q, LDA, eps);
            /* 
             * WARNING: For now, Q is associated to Band tridiagonal reduction and 
             * not to the final tridiagonal reduction, so we can not call the check
             */
            if (0)
                info_reduction = check_reduction(itype[i], uplo[u], N, 1, A1, A2, LDA, B2, LDB, Q, eps);

            memcpy(A1, Ainit, LDAxN*sizeof(PLASMA_Complex32_t));
            memcpy(B1, Binit, LDBxN*sizeof(PLASMA_Complex32_t));

            LAPACKE_chegv( LAPACK_COL_MAJOR, 
                     itype[i], lapack_const(vec), lapack_const(uplo[u]),
                     N, A1, LDA, B1, LDB, W1);

            /*info_solution  = check_solution(N, N, N, A1, LDA, B1, B2, LDB, eps);*/
            info_solution = check_solution(N, W1, W2, eps);
         
            if ( (info_ortho == 0) & (info_reduction == 0) & (info_solution == 0)) {
                printf("***************************************************\n");
                printf(" ---- TESTING CHEGV (%s, %s) ...................... PASSED !\n", itypestr[i], uplostr[u]);
                printf("***************************************************\n");
            }
            else {
                printf("************************************************\n");
                printf(" - TESTING CHEGV (%s, %s) ... FAILED !\n", itypestr[i], uplostr[u]);
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
static int check_orthogonality(int M, int N, PLASMA_Complex32_t *Q, int LDQ, float eps)
{
    float  alpha =  1.0;
    float  beta  = -1.0;
    float  normQ, result;
    int     info_ortho;
    int     minMN = min(M, N);
    float *work = (float *)malloc(minMN*sizeof(float));

    /* Build the idendity matrix */
    PLASMA_Complex32_t *Id = (PLASMA_Complex32_t *) malloc(minMN*minMN*sizeof(PLASMA_Complex32_t));
    LAPACKE_claset_work(LAPACK_COL_MAJOR, 'A', minMN, minMN, 0., 1., Id, minMN);

    /* Perform Id - Q'Q */
    if (M >= N)
        cblas_cherk(CblasColMajor, CblasUpper, CblasConjTrans, N, M, alpha, Q, LDQ, beta, Id, N);
    else
        cblas_cherk(CblasColMajor, CblasUpper, CblasNoTrans,   M, N, alpha, Q, LDQ, beta, Id, M);

    normQ = LAPACKE_clansy_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), 'U', minMN, Id, minMN, work);

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
                           PLASMA_Complex32_t *A1, PLASMA_Complex32_t *A2, int LDA, 
                           PLASMA_Complex32_t *B2, int LDB, PLASMA_Complex32_t *Q, float eps )
{
    PLASMA_Complex32_t alpha = 1.0;
    PLASMA_Complex32_t beta  = 0.0;
    float Anorm, Rnorm, result;
    int info_reduction;
    int i, j;
    char *str;
    
    PLASMA_Complex32_t *Aorig    = (PLASMA_Complex32_t *)malloc(N*N*sizeof(PLASMA_Complex32_t));
    PLASMA_Complex32_t *Residual = (PLASMA_Complex32_t *)malloc(N*N*sizeof(PLASMA_Complex32_t));
    PLASMA_Complex32_t *T        = (PLASMA_Complex32_t *)malloc(N*N*sizeof(PLASMA_Complex32_t));
    float *work                 = (float *)malloc(N*sizeof(float));

    memset((void*)T, 0, N*N*sizeof(PLASMA_Complex32_t));
    LAPACKE_clacpy_work(LAPACK_COL_MAJOR, ' ', N, N, A1, LDA, Residual, N);

    /* Rebuild the T */
    LAPACKE_clacpy_work(LAPACK_COL_MAJOR, lapack_const(uplo), N, N, A2, LDA, T, N);

    if (uplo == PlasmaLower) {
        /* Set the reflectors to 0 */
        for (i = bw+1; i < N; i++)
            for (j = 0 ; (j < N) && (j < i-bw); j++)
                T[j*N+i] = 0.;

        /* Copy the lower part to the upper part to rebuild the symmetry */
        for (i = 0; i < N; i++)
            for (j = 0 ; j < i; j++)
                T[i*N+j] = conjf(T[j*N+i]);
    } else {
        /* Set the reflectors to 0 */
        for (j = bw+1; j < N; j++)
            for (i = 0 ; (i < N) && (i < j-bw); i++)
                T[j*N+i] = 0.;
        /* Copy the upper part to the lower part to rebuild the symmetry */
        for (i = 0; i < N; i++)
            for (j = i+1 ; j < N; j++)
                T[i*N+j] = conjf(T[j*N+i]);
    }

    memset((void*)Aorig, 0, N*N*sizeof(PLASMA_Complex32_t));

    if (itype == 1) {
        if (uplo == PlasmaLower) {
            str = "L*Q*T*Q'*L'";
   
            /* Compute Aorig=Q*T*Q' */
            cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,   N, N, N, CBLAS_SADDR(alpha), Q,     LDA, T, N,   CBLAS_SADDR(beta), Aorig, N);
            cblas_cgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, N, N, N, CBLAS_SADDR(alpha), Aorig, N,   Q, LDA, CBLAS_SADDR(beta), T,     N);

            LAPACKE_clacpy_work(LAPACK_COL_MAJOR, PlasmaUpperLower, N, N, T, N, Aorig, N);

            cblas_ctrmm(CblasColMajor, CblasLeft,  CblasLower, CblasNoTrans,   CblasNonUnit, N, N, CBLAS_SADDR(alpha), B2, LDB, Aorig, N);   
            cblas_ctrmm(CblasColMajor, CblasRight, CblasLower, CblasConjTrans, CblasNonUnit, N, N, CBLAS_SADDR(alpha), B2, LDB, Aorig, N); 

        }
        else {
            str = "U'*Q*T*Q'*U";

            /* Compute Aorig=Q'*T*Q */
            cblas_cgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, N, N, N, CBLAS_SADDR(alpha), Q,     LDA, T, N,   CBLAS_SADDR(beta), Aorig, N);
            cblas_cgemm(CblasColMajor, CblasNoTrans,   CblasNoTrans, N, N, N, CBLAS_SADDR(alpha), Aorig, N,   Q, LDA, CBLAS_SADDR(beta), T,     N);

            LAPACKE_clacpy_work(LAPACK_COL_MAJOR, 'A', N, N, T, N, Aorig, N);

            cblas_ctrmm(CblasColMajor, CblasLeft,  CblasUpper, CblasConjTrans, CblasNonUnit, N, N, CBLAS_SADDR(alpha), B2, LDB, Aorig, N);   
            cblas_ctrmm(CblasColMajor, CblasRight, CblasUpper, CblasNoTrans,   CblasNonUnit, N, N, CBLAS_SADDR(alpha), B2, LDB, Aorig, N);  
        }
    }
    else {
        if (uplo == PlasmaLower) {
            str = "inv(L')*Q*A2*Q'*inv(L)";
            
            /* Compute Aorig=Q*T*Q' */
            cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,   N, N, N, CBLAS_SADDR(alpha), Q,     LDA, T, N,   CBLAS_SADDR(beta), Aorig, N);
            cblas_cgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, N, N, N, CBLAS_SADDR(alpha), Aorig, N,   Q, LDA, CBLAS_SADDR(beta), T, N   );

            LAPACKE_clacpy_work(LAPACK_COL_MAJOR, 'A', N, N, T, N, Aorig, N);

            cblas_ctrsm(CblasColMajor, CblasLeft,  CblasLower, CblasConjTrans, CblasNonUnit, N, N, CBLAS_SADDR(alpha), B2, LDB, Aorig, N);   
            cblas_ctrsm(CblasColMajor, CblasRight, CblasLower, CblasNoTrans,   CblasNonUnit, N, N, CBLAS_SADDR(alpha), B2, LDB, Aorig, N);   
        }
        else{
            str = "inv(U)*Q*A2*Q'*inv(U')";

            /* Compute Aorig=Q'*T*Q */
            cblas_cgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, N, N, N, CBLAS_SADDR(alpha), Q,     LDA, T, N,   CBLAS_SADDR(beta), Aorig, N);
            cblas_cgemm(CblasColMajor, CblasNoTrans,   CblasNoTrans, N, N, N, CBLAS_SADDR(alpha), Aorig, N,   Q, LDA, CBLAS_SADDR(beta), T,     N);
      
            LAPACKE_clacpy_work(LAPACK_COL_MAJOR, 'A', N, N, T, N, Aorig, N);
          
            cblas_ctrsm(CblasColMajor, CblasLeft,  CblasUpper, CblasNoTrans,   CblasNonUnit, N, N, CBLAS_SADDR(alpha), B2, LDB, Aorig, N);   
            cblas_ctrsm(CblasColMajor, CblasRight, CblasUpper, CblasConjTrans, CblasNonUnit, N, N, CBLAS_SADDR(alpha), B2, LDB, Aorig, N);   
        }
    }
   
    /* Compute the Residual */
    for (i = 0; i < N; i++)
        for (j = 0 ; j < N; j++)
            Residual[j*N+i] = A1[j*LDA+i]-Aorig[j*N+i];
    
    Rnorm = LAPACKE_clange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, N, Residual, N,   work);
    Anorm = LAPACKE_clange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, N, A2,       LDA, work);
    
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
static int check_solution(int N, float *E1, float *E2, float eps)
{
    int info_solution, i;
    float resid;
    float maxtmp;
    float maxel = fabs( fabs(E1[0]) - fabs(E2[0]) );
    float maxeig = max( fabs(E1[0]), fabs(E2[0]) );
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
static int check_solution(int M, int N, int NRHS, PLASMA_Complex32_t *A1, int LDA, 
                          PLASMA_Complex32_t *B1, PLASMA_Complex32_t *B2, int LDB, float eps)
{
    PLASMA_Complex32_t alpha = 1.0;
    PLASMA_Complex32_t beta  = 0.0;
    PLASMA_Complex32_t *Residual;
    int info_solution;
    int maxMN = max(M, N);
    float Rnorm, Anorm, Xnorm, Bnorm, result;
    float *work = (float *)malloc( max(maxMN, NRHS)* sizeof(float));

    Anorm = LAPACKE_clange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, N,    A1, LDA, work);
    Xnorm = LAPACKE_clange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, NRHS, B2, LDB, work);
    Bnorm = LAPACKE_clange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, NRHS, B1, LDB, work);

    cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, NRHS, N, 
                CBLAS_SADDR(alpha), A1, LDA, 
                                    B2, LDB, 
                CBLAS_SADDR(beta),  B1, LDB);

    Residual = (PLASMA_Complex32_t *)malloc(maxMN*NRHS*sizeof(PLASMA_Complex32_t));
    memset((void*)Residual, 0, maxMN*NRHS*sizeof(PLASMA_Complex32_t));
    cblas_cgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, N, NRHS, M, 
                CBLAS_SADDR(alpha), A1, LDA, 
                                    B1, LDB, 
                CBLAS_SADDR(beta),  Residual, maxMN);
    
    Rnorm = LAPACKE_clange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), maxMN, NRHS, Residual, maxMN, work);
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
