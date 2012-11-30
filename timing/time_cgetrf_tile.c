/**
 *
 * @generated c Thu Sep 15 12:09:46 2011
 *
 **/
#define _TYPE  PLASMA_Complex32_t
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_cgetrf_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GETRF(m, n)
#define _FADDS FADDS_GETRF(m, n)

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    PLASMA_Complex32_t *A = NULL, *AT, *b = NULL, *bT, *x;
    PLASMA_desc        *descA, *descB;
    real_Double_t       t;
    int                *piv;
    int m     = iparam[TIMING_M];
    int n     = iparam[TIMING_N];
    int nb    = iparam[TIMING_NB];
    int nrhs  = iparam[TIMING_NRHS];
    int check = iparam[TIMING_CHECK];
    int lda   = m;
    int ldb   = m;

    /* Allocate Data */
    AT  = (PLASMA_Complex32_t *)malloc(lda*n*sizeof(PLASMA_Complex32_t));
    piv = (int *)malloc( min(m, n) * sizeof(int));

    /* Check if unable to allocate memory */
    if ( !AT || !piv ){
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
    
    /* Initialiaze Data */
    PLASMA_Desc_Create(&descA, AT, PlasmaComplexFloat, nb, nb, nb*nb, lda, n, 0, 0, m, n);
    /*LAPACKE_clarnv_work(1, ISEED, n*lda, AT);*/
    PLASMA_cplrnt_Tile(descA, 3456);

    /* Save AT in lapack layout for check */
    if ( check && (m == n) ) {
        A = (PLASMA_Complex32_t *)malloc(lda*n    *sizeof(PLASMA_Complex32_t));
        PLASMA_cTile_to_Lapack(descA, (void*)A, lda);
    }

    t = -cWtime();
    PLASMA_cgetrf_Tile( descA, piv );
    t += cWtime();
    *t_ = t;
    
    /* Check the solution */
    if ( check && (m == n) )
      {
        b  = (PLASMA_Complex32_t *)malloc(ldb*nrhs*sizeof(PLASMA_Complex32_t));
        bT = (PLASMA_Complex32_t *)malloc(ldb*nrhs*sizeof(PLASMA_Complex32_t));
        x  = (PLASMA_Complex32_t *)malloc(ldb*nrhs*sizeof(PLASMA_Complex32_t));

        LAPACKE_clarnv_work(1, ISEED, ldb*nrhs, b);

        PLASMA_Desc_Create(&descB, bT, PlasmaComplexFloat, nb, nb, nb*nb, ldb, nrhs, 0, 0, n, nrhs);
        PLASMA_cLapack_to_Tile((void*)b, ldb, descB);

        PLASMA_cgetrs_Tile( PlasmaNoTrans, descA, piv, descB );

        PLASMA_cTile_to_Lapack(descB, (void*)x, ldb);
        PLASMA_Desc_Destroy(&descB);

        dparam[TIMING_RES] = c_check_solution(m, n, nrhs, A, lda, b, x, ldb,
                                              &(dparam[TIMING_ANORM]), &(dparam[TIMING_BNORM]), 
                                              &(dparam[TIMING_XNORM]));
        
        free( A ); free( b ); free( bT ); free( x );
      }

    PLASMA_Desc_Destroy(&descA);

    free( AT );
    free( piv );
    PLASMA_Finalize();

    return 0;
}
