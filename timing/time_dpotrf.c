/**
 *
 * @generated d Thu Sep 15 12:09:45 2011
 *
 **/
#define _TYPE  double
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_dpotrf"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_POTRF( n )
#define _FADDS FADDS_POTRF( n )

#include "./timing.c"

static void generate_matrix( double* A, size_t N )
{
  size_t i, j;
  for (i = 0; i< N; i++) {
    A[i*N+i] = A[i*N+i] + 1.*N; 
    for (j = 0; j < i; j++)
      A[i*N+j] = A[j*N+i];
  }
}

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    double *A, *Acpy = NULL, *b, *x;
    real_Double_t       t;
    int n     = iparam[TIMING_N];
    int nrhs  = iparam[TIMING_NRHS];
    int check = iparam[TIMING_CHECK];
    int lda = n;
    int ldb = n;
    int uplo = PlasmaLower;

    /* Allocate Data */
    A = (double *)malloc(lda*n*   sizeof(double));

    /* Check if unable to allocate memory */
    if ( !A ) {
        printf("Out of Memory \n ");
        exit(0);
    }
    
    /* TODO kaapi_init */
    /* Initialize Plasma */ 
    PLASMA_Init( iparam[TIMING_THRDNBR] );
    PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_DYNAMIC_SCHEDULING );

#if defined(PLASMA_CUDA)
    cudaHostRegister(A, lda*n*sizeof(double), cudaHostRegisterPortable);
#endif
  
#if defined(PLASMA_CUDA)
    core_cublas_init();
#endif

    /*if ( !iparam[TIMING_AUTOTUNING] ) {*/
    PLASMA_Disable(PLASMA_AUTOTUNING);
    PLASMA_Set(PLASMA_TILE_SIZE, iparam[TIMING_NB] );
    /* } */

    /* Initialiaze Data */
    PLASMA_dplgsy( (double)n, n, A, lda, 51 );
    /* 
    int ISEED[4] = {0,0,0,1};
    LAPACKE_dlarnv(1, ISEED, n*n, A);
    generate_matrix( A, n );
    */

    /* Save A and b  */
    if (check) {
        Acpy = (double *)malloc(lda*n*sizeof(double));
        LAPACKE_dlacpy_work(LAPACK_COL_MAJOR,' ', n, n, A, lda, Acpy, lda);
    }

    /* PLASMA DPOSV */
    t = -cWtime();
    PLASMA_dpotrf(uplo, n, A, lda);
    t += cWtime();
    *t_ = t;

    /* Check the solution */
    if (check)
      {
        b = (double *)malloc(ldb*nrhs*sizeof(double));
        x = (double *)malloc(ldb*nrhs*sizeof(double));
        LAPACKE_dlarnv_work(1, ISEED, n*nrhs, x);
        LAPACKE_dlacpy_work(LAPACK_COL_MAJOR, 'A', n, nrhs, x, ldb, b, ldb);

        PLASMA_dpotrs(uplo, n, nrhs, A, lda, x, ldb);

        dparam[TIMING_RES] = d_check_solution(n, n, nrhs, Acpy, lda, b, x, ldb,
                                             &(dparam[TIMING_ANORM]), 
                                             &(dparam[TIMING_BNORM]),
                                             &(dparam[TIMING_XNORM]));

        free(Acpy); free(b); free(x);
      }

#if defined(PLASMA_CUDA)
    cudaHostUnregister(A);
#endif
    free(A);

    PLASMA_Finalize();

    return 0;
}
