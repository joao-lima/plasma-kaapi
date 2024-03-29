/**
 *
 * @generated s Thu Sep 15 12:09:45 2011
 *
 **/
#define _TYPE  float
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_sgetrf_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GETRF(m, n)
#define _FADDS FADDS_GETRF(m, n)

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    float *A, *Acpy = NULL, *b, *x;
    real_Double_t       t;
    int                *piv;
    int m     = iparam[TIMING_M];
    int n     = iparam[TIMING_N];
    int nrhs  = iparam[TIMING_NRHS];
    int check = iparam[TIMING_CHECK];
    int lda   = m;
    int ldb   = m;

    /* Allocate Data */
    A   = (float *)malloc(lda*n*sizeof(float));
    piv = (int *)malloc( min(m, n) * sizeof(int));
    
    /* Check if unable to allocate memory */
    if ( !A || !piv ){
        printf("Out of Memory \n ");
        return -1;
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
    
    /* Initialize Data */
    /*LAPACKE_slarnv_work(1, ISEED, n*lda, A);*/
    PLASMA_splrnt(m, n, A, lda, 3456);

    /* Save AT in lapack layout for check */
    if ( check && (m == n) ) {
        Acpy = (float *)malloc(lda*n*sizeof(float));
        LAPACKE_slacpy_work(LAPACK_COL_MAJOR, 'A', m, n, A, lda, Acpy, lda);
    }

    t = -cWtime();
    PLASMA_sgetrf( m, n, A, lda, piv );
    t += cWtime();
    *t_ = t;
    
    /* Check the solution */
    if ( check && (m == n) )
      {
        b  = (float *)malloc(ldb*nrhs *sizeof(float));
        x  = (float *)malloc(ldb*nrhs *sizeof(float));

        LAPACKE_slarnv_work(1, ISEED, ldb*nrhs, x);
        LAPACKE_slacpy_work(LAPACK_COL_MAJOR, 'A', n, nrhs, x, ldb, b, ldb);

        PLASMA_sgetrs( PlasmaNoTrans, n, nrhs, A, lda, piv, x, ldb );

        dparam[TIMING_RES] = s_check_solution(m, n, nrhs, Acpy, lda, b, x, ldb,
                                              &(dparam[TIMING_ANORM]), &(dparam[TIMING_BNORM]), 
                                              &(dparam[TIMING_XNORM]));

        free( Acpy ); free( b ); free( x );
      }

    free( A );
    free( piv );
    PLASMA_Finalize();

    return 0;
}
