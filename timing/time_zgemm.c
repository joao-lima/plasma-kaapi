/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_zgemm"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GEMM(n, n, n)
#define _FADDS FADDS_GEMM(n, n, n)

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    PLASMA_Complex64_t *A, *B, *C1;
    PLASMA_Complex64_t *C2 = NULL;
    PLASMA_Complex64_t alpha, beta;
    real_Double_t       t;
    int n       = iparam[TIMING_N];
    int check   = iparam[TIMING_CHECK];
    int lda     = n;
    
    /* Allocate Data */
    A  = (PLASMA_Complex64_t *)malloc(lda*n*   sizeof(PLASMA_Complex64_t));
    B  = (PLASMA_Complex64_t *)malloc(lda*n*   sizeof(PLASMA_Complex64_t));
    C1 = (PLASMA_Complex64_t *)malloc(lda*n*   sizeof(PLASMA_Complex64_t));

    LAPACKE_zlarnv_work(1, ISEED, 1,  &alpha);
    LAPACKE_zlarnv_work(1, ISEED, 1,  &beta);

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
    LAPACKE_zlarnv_work(1, ISEED, n*lda,  A);
    LAPACKE_zlarnv_work(1, ISEED, n*lda,  B);
    LAPACKE_zlarnv_work(1, ISEED, n*lda,  C1);

    if (check)
      {
          C2 = (PLASMA_Complex64_t *)malloc(lda*n*   sizeof(PLASMA_Complex64_t));
          memcpy(C2, C1, lda*n*sizeof(PLASMA_Complex64_t));
      }

    t = -cWtime();
    PLASMA_zgemm( PlasmaNoTrans, PlasmaNoTrans, n, n, n, alpha, A, lda, B, lda, beta, C1, lda );
    t += cWtime();
    *t_ = t;
    
    /* Check the solution */
    if (check)
      {
          dparam[TIMING_RES] = z_check_gemm( PlasmaNoTrans, PlasmaNoTrans, n, n, n, 
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
