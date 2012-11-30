/**
 *
 * @generated d Thu Sep 15 12:09:45 2011
 *
 **/
#define _TYPE  double
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_dgemm"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GEMM(n, n, n)
#define _FADDS FADDS_GEMM(n, n, n)

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    double *A, *B, *C1;
    double *C2 = NULL;
    double alpha, beta;
    real_Double_t       t;
    int n       = iparam[TIMING_N];
    int check   = iparam[TIMING_CHECK];
    int lda     = n;
    
    /* Allocate Data */
    A  = (double *)malloc(lda*n*   sizeof(double));
    B  = (double *)malloc(lda*n*   sizeof(double));
    C1 = (double *)malloc(lda*n*   sizeof(double));

    LAPACKE_dlarnv_work(1, ISEED, 1,  &alpha);
    LAPACKE_dlarnv_work(1, ISEED, 1,  &beta);

    /* Check if unable to allocate memory */
    if ( (!A) || (!B) || (!C1) ) {
        printf("Out of Memory \n ");
        exit(0);
    }
    
    /* Initialize Plasma */ 
    PLASMA_Init( iparam[TIMING_THRDNBR] );
    if ( iparam[TIMING_SCHEDULER] )
        PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_DYNAMIC_SCHEDULING );
    else
        PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_STATIC_SCHEDULING );

    /*if ( !iparam[TIMING_AUTOTUNING] ) {*/
        PLASMA_Disable(PLASMA_AUTOTUNING);
        PLASMA_Set(PLASMA_TILE_SIZE, iparam[TIMING_NB] );
    /* } */

     /* Initialiaze Data */
    LAPACKE_dlarnv_work(1, ISEED, n*lda,  A);
    LAPACKE_dlarnv_work(1, ISEED, n*lda,  B);
    LAPACKE_dlarnv_work(1, ISEED, n*lda,  C1);

    if (check)
      {
          C2 = (double *)malloc(lda*n*   sizeof(double));
          memcpy(C2, C1, lda*n*sizeof(double));
      }

    t = -cWtime();
    PLASMA_dgemm( PlasmaNoTrans, PlasmaNoTrans, n, n, n, alpha, A, lda, B, lda, beta, C1, lda );
    t += cWtime();
    *t_ = t;
    
    /* Check the solution */
    if (check)
      {
          dparam[TIMING_RES] = d_check_gemm( PlasmaNoTrans, PlasmaNoTrans, n, n, n, 
                                            alpha, A, lda, B, lda, beta, C1, C2, lda,
                                            &(dparam[TIMING_ANORM]), &(dparam[TIMING_BNORM]), 
                                            &(dparam[TIMING_XNORM]));
          free(C2);
      }

    free( A );
    free( B );
    free( C1 );

    PLASMA_Finalize();

    return 0;
}
