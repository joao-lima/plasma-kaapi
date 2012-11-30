/**
 *
 * @generated s Thu Sep 15 12:09:46 2011
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
    float *A = NULL, *AT, *b = NULL, *bT, *x;
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
    AT  = (float *)malloc(lda*n*sizeof(float));
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
    PLASMA_Desc_Create(&descA, AT, PlasmaRealFloat, nb, nb, nb*nb, lda, n, 0, 0, m, n);
    /*LAPACKE_slarnv_work(1, ISEED, n*lda, AT);*/
    PLASMA_splrnt_Tile(descA, 3456);

    /* Save AT in lapack layout for check */
    if ( check && (m == n) ) {
        A = (float *)malloc(lda*n    *sizeof(float));
        PLASMA_sTile_to_Lapack(descA, (void*)A, lda);
    }

    t = -cWtime();
    PLASMA_sgetrf_Tile( descA, piv );
    t += cWtime();
    *t_ = t;
    
    /* Check the solution */
    if ( check && (m == n) )
      {
        b  = (float *)malloc(ldb*nrhs*sizeof(float));
        bT = (float *)malloc(ldb*nrhs*sizeof(float));
        x  = (float *)malloc(ldb*nrhs*sizeof(float));

        LAPACKE_slarnv_work(1, ISEED, ldb*nrhs, b);

        PLASMA_Desc_Create(&descB, bT, PlasmaRealFloat, nb, nb, nb*nb, ldb, nrhs, 0, 0, n, nrhs);
        PLASMA_sLapack_to_Tile((void*)b, ldb, descB);

        PLASMA_sgetrs_Tile( PlasmaNoTrans, descA, piv, descB );

        PLASMA_sTile_to_Lapack(descB, (void*)x, ldb);
        PLASMA_Desc_Destroy(&descB);

        dparam[TIMING_RES] = s_check_solution(m, n, nrhs, A, lda, b, x, ldb,
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
