/**
 *
 * @generated c Thu Sep 15 12:09:46 2011
 *
 **/
#define _TYPE  PLASMA_Complex32_t
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_cgetrf_incpiv_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GETRF(n, n)
#define _FADDS FADDS_GETRF(n, n)

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    PLASMA_Complex32_t *A, *Acpy = NULL, *L, *b, *x;
    real_Double_t       t;
    int                *piv;
    int n     = iparam[TIMING_N];
    int nrhs  = iparam[TIMING_NRHS];
    int check = iparam[TIMING_CHECK];
    int lda      = n;
    int ldb      = n;

    /* Allocate Data */
    A = (PLASMA_Complex32_t *)malloc(lda*n*sizeof(PLASMA_Complex32_t));
    
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
    LAPACKE_clarnv_work(1, ISEED, n*lda, A);

    /* Allocate Workspace */
    PLASMA_Alloc_Workspace_cgesv_incpiv(n, &L, &piv);

    /* Save AT in lapack layout for check */
    if ( check ) {
        Acpy = (PLASMA_Complex32_t *)malloc(lda*n*sizeof(PLASMA_Complex32_t));
        LAPACKE_clacpy_work(LAPACK_COL_MAJOR,' ', n, n, A, lda, Acpy, lda);
    }

    t = -cWtime();
    PLASMA_cgetrf_incpiv( n, n, A, lda, L, piv );
    t += cWtime();
    *t_ = t;
    
    /* Check the solution */
    if ( check )
      {
        b  = (PLASMA_Complex32_t *)malloc(ldb*nrhs *sizeof(PLASMA_Complex32_t));
        x  = (PLASMA_Complex32_t *)malloc(ldb*nrhs *sizeof(PLASMA_Complex32_t));

        LAPACKE_clarnv_work(1, ISEED, ldb*nrhs, x);
        LAPACKE_clacpy_work(LAPACK_COL_MAJOR,' ', n, nrhs, x, ldb, b, ldb);

        PLASMA_cgetrs_incpiv( PlasmaNoTrans, n, nrhs, A, lda, L, piv, x, ldb );

        dparam[TIMING_RES] = c_check_solution(n, n, nrhs, Acpy, lda, b, x, ldb,
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
