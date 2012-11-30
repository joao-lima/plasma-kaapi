/**
 *
 * @generated s Thu Sep 15 12:09:46 2011
 *
 **/
#define _TYPE  float
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_sgels"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_GEQRF( n, n ) + FMULS_GEQRS( n, n, nrhs ))
#define _FADDS (FADDS_GEQRF( n, n ) + FADDS_GEQRS( n, n, nrhs ))

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    float *A, *x, *T;
    float *Acpy = NULL;
    float *b = NULL;
    real_Double_t       t;
    int n       = iparam[TIMING_N];
    int nrhs    = iparam[TIMING_NRHS];
    int check   = iparam[TIMING_CHECK];
    int lda     = n;
    int ldb     = n;
    
    /* Allocate Data */
    A = (float *)malloc(lda*n*   sizeof(float));
    x = (float *)malloc(ldb*nrhs*sizeof(float));

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
        PLASMA_Set(PLASMA_TILE_SIZE,        iparam[TIMING_NB] );
        PLASMA_Set(PLASMA_INNER_BLOCK_SIZE, iparam[TIMING_IB] );
    /* } */

     /* Initialiaze Data */
    LAPACKE_slarnv_work(1, ISEED, n*lda,  A);
    LAPACKE_slarnv_work(1, ISEED, n*nrhs, x);

    PLASMA_Alloc_Workspace_sgels(n, n, &T);

    /* Save A and b  */
    if (check) {
        Acpy = (float *)malloc(lda*n*   sizeof(float));
        b    = (float *)malloc(ldb*nrhs*sizeof(float));
        LAPACKE_slacpy_work(LAPACK_COL_MAJOR,' ', n, n,    A, lda, Acpy, lda);
        LAPACKE_slacpy_work(LAPACK_COL_MAJOR,' ', n, nrhs, x, ldb, b,    ldb);
      }

    t = -cWtime();
    PLASMA_sgels( PlasmaNoTrans, n, n, nrhs, A, lda, T, x, ldb );
    t += cWtime();
    *t_ = t;
    
    /* Check the solution */
    if (check)
      {
        dparam[TIMING_RES] = s_check_solution(n, n, nrhs, Acpy, lda, b, x, ldb,
                                             &(dparam[TIMING_ANORM]), &(dparam[TIMING_BNORM]), 
                                             &(dparam[TIMING_XNORM]));
        free(Acpy); free(b);
      }

    free( T );
    free( A );
    free( x );

    PLASMA_Finalize();

    return 0;
}
