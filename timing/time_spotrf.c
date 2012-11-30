/**
 *
 * @generated s Thu Sep 15 12:09:45 2011
 *
 **/
#define _TYPE  float
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_spotrf"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_POTRF( n )
#define _FADDS FADDS_POTRF( n )

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    float *A, *Acpy = NULL, *b, *x;
    real_Double_t       t;
    int n     = iparam[TIMING_N];
    int nrhs  = iparam[TIMING_NRHS];
    int check = iparam[TIMING_CHECK];
    int lda = n;
    int ldb = n;
    int uplo = PlasmaLower;

    /* Allocate Data */
    A = (float *)malloc(lda*n*   sizeof(float));

    /* Check if unable to allocate memory */
    if ( !A ) {
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
    PLASMA_splgsy( (float)n, n, A, lda, 51 );

    /* Save A and b  */
    if (check) {
        Acpy = (float *)malloc(lda*n*sizeof(float));
        LAPACKE_slacpy_work(LAPACK_COL_MAJOR,' ', n, n, A, lda, Acpy, lda);
    }

    /* PLASMA SPOSV */
    t = -cWtime();
    PLASMA_spotrf(uplo, n, A, lda);
    t += cWtime();
    *t_ = t;

    /* Check the solution */
    if (check)
      {
        b = (float *)malloc(ldb*nrhs*sizeof(float));
        x = (float *)malloc(ldb*nrhs*sizeof(float));
        LAPACKE_slarnv_work(1, ISEED, n*nrhs, x);
        LAPACKE_slacpy_work(LAPACK_COL_MAJOR, 'A', n, nrhs, x, ldb, b, ldb);

        PLASMA_spotrs(uplo, n, nrhs, A, lda, x, ldb);

        dparam[TIMING_RES] = s_check_solution(n, n, nrhs, Acpy, lda, b, x, ldb,
                                             &(dparam[TIMING_ANORM]), 
                                             &(dparam[TIMING_BNORM]),
                                             &(dparam[TIMING_XNORM]));

        free(Acpy); free(b); free(x);
      }

    free(A);

    PLASMA_Finalize();

    return 0;
}
