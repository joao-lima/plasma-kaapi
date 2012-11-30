/**
 *
 * @generated s Thu Sep 15 12:09:46 2011
 *
 **/
#define _TYPE  float
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_sgetrf_incpiv_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GETRF(n, n)
#define _FADDS FADDS_GETRF(n, n)

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    float *A, *Acpy = NULL, *L, *b, *x;
    real_Double_t       t;
    int                *piv;
    int n     = iparam[TIMING_N];
    int nrhs  = iparam[TIMING_NRHS];
    int check = iparam[TIMING_CHECK];
    int lda      = n;
    int ldb      = n;

    /* Allocate Data */
    A = (float *)malloc(lda*n*sizeof(float));
    
    /* Check if unable to allocate memory */
    if ( !A ){
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
    /* } else { */
    /*     PLASMA_Get(PLASMA_TILE_SIZE,        &iparam[TIMING_NB] ); */
    /*     PLASMA_Get(PLASMA_INNER_BLOCK_SIZE, &iparam[TIMING_IB] ); */
    /* }  */
    
    /* Initialiaze Data */
    LAPACKE_slarnv_work(1, ISEED, n*lda, A);

    /* Allocate Workspace */
    PLASMA_Alloc_Workspace_sgesv_incpiv(n, &L, &piv);

    /* Save AT in lapack layout for check */
    if ( check ) {
        Acpy = (float *)malloc(lda*n*sizeof(float));
        LAPACKE_slacpy_work(LAPACK_COL_MAJOR,' ', n, n, A, lda, Acpy, lda);
    }

    t = -cWtime();
    PLASMA_sgetrf_incpiv( n, n, A, lda, L, piv );
    t += cWtime();
    *t_ = t;
    
    /* Check the solution */
    if ( check )
      {
        b  = (float *)malloc(ldb*nrhs *sizeof(float));
        x  = (float *)malloc(ldb*nrhs *sizeof(float));

        LAPACKE_slarnv_work(1, ISEED, ldb*nrhs, x);
        LAPACKE_slacpy_work(LAPACK_COL_MAJOR,' ', n, nrhs, x, ldb, b, ldb);

        PLASMA_sgetrs_incpiv( PlasmaNoTrans, n, nrhs, A, lda, L, piv, x, ldb );

        dparam[TIMING_RES] = s_check_solution(n, n, nrhs, Acpy, lda, b, x, ldb,
                                             &(dparam[TIMING_ANORM]), &(dparam[TIMING_BNORM]), 
                                             &(dparam[TIMING_XNORM]));

        free( Acpy ); free( b ); free( x );
      }

    free( A );
    free( L );
    free( piv );
    PLASMA_Finalize();

    return 0;
}
