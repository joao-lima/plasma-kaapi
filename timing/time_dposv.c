/**
 *
 * @generated d Thu Sep 15 12:09:45 2011
 *
 **/
#define _TYPE  double
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_dposv"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_POTRF( n ) + FMULS_POTRS( n, nrhs ))
#define _FADDS (FADDS_POTRF( n ) + FADDS_POTRS( n, nrhs ))

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    double *A, *Acpy = NULL, *b = NULL, *x;
    real_Double_t       t;
    int n     = iparam[TIMING_N];
    int nrhs  = iparam[TIMING_NRHS];
    int check = iparam[TIMING_CHECK];
    int lda = n;
    int ldb = n;

    /* Allocate Data */
    A = (double *)malloc(lda*n*   sizeof(double));
    x = (double *)malloc(ldb*nrhs*sizeof(double));

    /* Check if unable to allocate memory */
    if ( (!A) || (!x) ) {
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
    PLASMA_dplgsy((double)n, n, A, lda, 51 );
    LAPACKE_dlarnv_work(1, ISEED, n*nrhs, x);

    /* Save A and b  */
    if (check) {
        Acpy = (double *)malloc(lda*n*   sizeof(double));
        b    = (double *)malloc(ldb*nrhs*sizeof(double));
        LAPACKE_dlacpy_work(LAPACK_COL_MAJOR,' ', n, n,    A, lda, Acpy, lda);
        LAPACKE_dlacpy_work(LAPACK_COL_MAJOR,' ', n, nrhs, x, ldb, b,    ldb);
      }

    /* PLASMA DPOSV */
    t = -cWtime();
    PLASMA_dposv(PlasmaUpper, n, nrhs, A, lda, x, ldb);
    t += cWtime();
    *t_ = t;

    /* Check the solution */
    if (check)
      {
        dparam[TIMING_RES] = d_check_solution(n, n, nrhs, Acpy, lda, b, x, ldb,
                                             &(dparam[TIMING_ANORM]), &(dparam[TIMING_BNORM]), 
                                             &(dparam[TIMING_XNORM]));
        free(Acpy); free(b);
      }

    free(A); free(x); 

    PLASMA_Finalize();

    return 0;
}
