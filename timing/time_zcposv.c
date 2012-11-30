/**
 *
 * @precisions mixed zc -> ds
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_zposv"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_POTRF( n ) + FMULS_POTRS( n, nrhs ))
#define _FADDS (FADDS_POTRF( n ) + FADDS_POTRS( n, nrhs ))

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    PLASMA_Complex64_t *A, *b, *x;
    real_Double_t       t;
    int n     = iparam[TIMING_N];
    int nrhs  = iparam[TIMING_NRHS];
    int check = iparam[TIMING_CHECK];
    int lda = n;
    int ldb = n;
    int iter;

    /* Allocate Data */
    A = (PLASMA_Complex64_t *)malloc(lda*n*   sizeof(PLASMA_Complex64_t));
    b = (PLASMA_Complex64_t *)malloc(ldb*nrhs*sizeof(PLASMA_Complex64_t));
    x = (PLASMA_Complex64_t *)malloc(ldb*nrhs*sizeof(PLASMA_Complex64_t));

    /* Check if unable to allocate memory */
    if ( (!A) || (!b) || (!x) ) {
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
    PLASMA_zplghe((double)n, n, A, lda, 51 );
    LAPACKE_zlarnv_work(1, ISEED, n*nrhs, b);

    /* PLASMA ZCPOSV */
    t = -cWtime();
    PLASMA_zcposv(PlasmaUpper, n, nrhs, A, lda, b, ldb, x, ldb, &iter);
    t += cWtime();
    *t_ = t;

    /* Check the solution */
    if (check)
      {
        dparam[TIMING_RES] = z_check_solution(n, n, nrhs, A, lda, b, x, ldb,
                                             &(dparam[TIMING_ANORM]), &(dparam[TIMING_BNORM]), 
                                             &(dparam[TIMING_XNORM]));
      }

    free(A); free(b); free(x); 

    PLASMA_Finalize();

    return 0;
}
