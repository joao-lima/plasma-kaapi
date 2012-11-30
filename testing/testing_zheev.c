/**
 *
 * @file testing_zheev.c
 *
 *  PLASMA testing routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Hatem Ltaief
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

#undef REAL
#define COMPLEX

static int check_orthogonality(int, int, PLASMA_Complex64_t*, int, double);
#if 0
static int check_reduction(int, int, int, PLASMA_Complex64_t*, PLASMA_Complex64_t*, int, PLASMA_Complex64_t*, double);
#endif
static int check_solution(int, double*, double*, double);

int testing_zheev(int argc, char **argv)
{
    /* Check for number of arguments*/
    if (argc != 2) {
        USAGE("HEEV", "N LDA",
              "   - N    : size of the matrix A\n"
              "   - LDA  : leading dimension of the matrix A\n");
        return -1;
    }

    int    N     = atoi(argv[0]);
    int    LDA   = atoi(argv[1]);
    int    LDQ   = LDA;
    int    mode  = 4;
    double eps   = LAPACKE_dlamch_work('e');
    double dmax  = 1.0;
    double rcond = 1.0e6;

    PLASMA_enum uplo = PlasmaLower;
    PLASMA_enum vec  = PlasmaNoVec;
    int info_ortho     = 0;
    int info_solution  = 0;
    int info_reduction = 0;
    int i;
    int LDAxN = LDA*N;
    int LDQxN = LDQ*N;

    PLASMA_Complex64_t *A1   = NULL;
    PLASMA_Complex64_t *A2   = (PLASMA_Complex64_t *)malloc(LDAxN*sizeof(PLASMA_Complex64_t));
    PLASMA_Complex64_t *Q    = NULL;
    double             *W1   = (double *)malloc(N*sizeof(double));
    double             *W2   = (double *)malloc(N*sizeof(double));
    PLASMA_Complex64_t *work = (PLASMA_Complex64_t *)malloc(3*N* sizeof(PLASMA_Complex64_t));
    PLASMA_desc *T;

    /* Check if unable to allocate memory */
    if ( (!A2) || (!W1) || (!W2) ){
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
    PLASMA_Alloc_Workspace_zheev(N, N, &T);

    /*----------------------------------------------------------
    *  TESTING ZHEEV
    */

    /* Initialize A1 */
    LAPACKE_zlatms_work( LAPACK_COL_MAJOR, N, N,
                         lapack_const(PlasmaDistSymmetric), ISEED,
                         lapack_const(PlasmaHermGeev), W1, mode, rcond,
                         dmax, N, N,
                         lapack_const(PlasmaNoPacking), A2, LDA, work );
    
    /*
     * Sort the eigenvalue because when computing the tridiag
     * and then the eigenvalue of the DSTQR are sorted.
     * So to avoid testing fail when having good results W1 should be sorted
    */
    LAPACKE_dlasrt_work( 'I', N, W1 );

    if (getenv("PLASMA_TESTING_VERBOSE"))
    {
        printf("Eigenvalues original\n");
        for (i = 0; i < min(N,25); i++){
            printf("%f \n", W1[i]);
        }
        printf("\n");
    }

    if ( vec == PlasmaVec ) {
        Q  = (PLASMA_Complex64_t *)malloc(LDQxN*sizeof(PLASMA_Complex64_t));
        A1 = (PLASMA_Complex64_t *)malloc(LDAxN*sizeof(PLASMA_Complex64_t));    

        /* Copy A2 into A1 */
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', N, N, A2, LDA, A1, LDA);
    }

    /* PLASMA ZHEEV */
    PLASMA_zheev(vec, uplo, N, A2, LDA, W2, T, Q, LDQ);

    if (getenv("PLASMA_TESTING_VERBOSE"))
    {
        printf("Eigenvalues computed\n");
        for (i = 0; i < min(N,25); i++){
            printf("%f \n", W2[i]);
        }
        printf("\n");
    }

    printf("\n");
    printf("------ TESTS FOR PLASMA ZHEEV ROUTINE -------  \n");
    printf("        Size of the Matrix %d by %d\n", N, N);
    printf("\n");
    printf(" The matrix A is randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n",eps);
    printf(" Computational tests pass if scaled residuals are less than 60.\n");

    /* Check the orthogonality, reduction and the eigen solutions */
    if (vec == PlasmaVec) {
        info_ortho = check_orthogonality(N, N, Q, LDQ, eps);
        /* 
         * WARNING: For now, Q is associated to Band tridiagonal reduction and 
         * not to the final tridiagonal reduction, so we can not call the check
         */
        /*info_reduction = check_reduction(uplo, N, 1, A1, A2, LDA, Q, eps);*/
    }
    info_solution = check_solution(N, W1, W2, eps);

    if ( (info_solution == 0) & (info_ortho == 0) & (info_reduction == 0) ) {
        printf("***************************************************\n");
        printf(" ---- TESTING ZHEEV ...................... PASSED !\n");
        printf("***************************************************\n");
    }
    else {
        printf("************************************************\n");
        printf(" - TESTING ZHEEV ... FAILED !\n");
        printf("************************************************\n");
    }

    PLASMA_Dealloc_Handle_Tile(&T);
    free(A2);
    free(W1);
    free(W2);
    free(work);
    if (Q  != NULL) free(Q);
    if (A1 != NULL) free(A1);

    return 0;
}

/*-------------------------------------------------------------------
 * Check the orthogonality of Q
 */

static int check_orthogonality(int M, int N, PLASMA_Complex64_t *Q, int LDQ, double eps)
{
    double  alpha =  1.0;
    double  beta  = -1.0;
    double  normQ, result;
    int     info_ortho;
    int     minMN = min(M, N);
    double *work = (double *)malloc(minMN*sizeof(double));

    /* Build the idendity matrix */
    PLASMA_Complex64_t *Id = (PLASMA_Complex64_t *) malloc(minMN*minMN*sizeof(PLASMA_Complex64_t));
    LAPACKE_zlaset_work(LAPACK_COL_MAJOR, 'A', minMN, minMN, 0., 1., Id, minMN);

    /* Perform Id - Q'Q */
    if (M >= N)
        cblas_zherk(CblasColMajor, CblasUpper, CblasConjTrans, N, M, alpha, Q, LDQ, beta, Id, N);
    else
        cblas_zherk(CblasColMajor, CblasUpper, CblasNoTrans,   M, N, alpha, Q, LDQ, beta, Id, M);

    normQ = LAPACKE_zlansy_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), 'U', minMN, Id, minMN, work);

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
 *  Check the reduction 
 */
#if 0
static int check_reduction(int uplo, int N, int bw, PLASMA_Complex64_t *A1, PLASMA_Complex64_t *A2, int LDA, PLASMA_Complex64_t *Q, double eps )
{
    PLASMA_Complex64_t alpha = 1.0;
    PLASMA_Complex64_t beta  = 0.0;
    double Anorm, Rnorm, result;
    int info_reduction;
    int i, j;

    PLASMA_Complex64_t *Aorig    = (PLASMA_Complex64_t *)malloc(N*N*sizeof(PLASMA_Complex64_t));
    PLASMA_Complex64_t *Residual = (PLASMA_Complex64_t *)malloc(N*N*sizeof(PLASMA_Complex64_t));
    PLASMA_Complex64_t *T        = (PLASMA_Complex64_t *)malloc(N*N*sizeof(PLASMA_Complex64_t));
    double *work = (double *)malloc(N*sizeof(double));

    memset((void*)T, 0, N*N*sizeof(PLASMA_Complex64_t));

    /* Rebuild the T */
    LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, lapack_const(uplo), N, N, A2, LDA, T, N);

    printf("============\n");
    printf("Checking the tridiagonal reduction \n");

    if (uplo == PlasmaLower) 
    {
        /* Copy the lower part to the upper part to rebuild the hermitian/symmetry */
        for (i = 0; i < N; i++)
            for (j = max(0, i-bw); j < i; j++)
                T[i*N+j] = conj(T[j*N+i]);
        
        /* Compute Aorig = Q * T * Q^H */
        memset((void*)Aorig, 0, N*N*sizeof(PLASMA_Complex64_t));
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,   N, N, N, CBLAS_SADDR(alpha), Q,    LDA, T, N,   CBLAS_SADDR(beta), Aorig, N);
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, N, N, N, CBLAS_SADDR(alpha), Aorig, N,  Q, LDA, CBLAS_SADDR(beta), T,     N);
    } 
    else
    {
        /* Copy the upper part to the lower part to rebuild the symmetry */
        for (i = 0; i < N; i++)
            for (j = i+1 ; j < min(i+bw, N); j++)
                T[i*N+j] = conj(T[j*N+i]);

        /* Compute Aorig = Q^H * T * Q */
        memset((void*)Aorig, 0, N*N*sizeof(PLASMA_Complex64_t));
        cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, N, N, N, CBLAS_SADDR(alpha), Q,     LDA, T, N,   CBLAS_SADDR(beta), Aorig, N);
        cblas_zgemm(CblasColMajor, CblasNoTrans,   CblasNoTrans, N, N, N, CBLAS_SADDR(alpha), Aorig, N,   Q, LDA, CBLAS_SADDR(beta), T,     N);
    }

    /* Compute the Residual */
    for (i = 0; i < N; i++)
        for (j = 0 ; j < N; j++)
            Residual[j*N+i] = A1[j*LDA+i] - T[j*N+i];
    
    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, N, Residual, N,   work);
    Anorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, N, A2,       LDA, work);
    
    result = Rnorm / ( Anorm * N * eps);
    if ( uplo == PlasmaLower )
        printf("-- ||A-Q*T*Q'||_oo/(||A||_oo.N.eps) = %e \n", result);
    else 
        printf("-- ||A-Q'*T*Q||_oo/(||A||_oo.N.eps) = %e \n", result);

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        printf("-- Reduction is suspicious ! \n");
        info_reduction = 1;
    }
    else {
        printf("-- Reduction is CORRECT ! \n");
        info_reduction = 0;
    }

    free(Aorig); free(Residual); free(T);

    return info_reduction;
}
#endif

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
