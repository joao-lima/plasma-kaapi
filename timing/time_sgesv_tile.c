/**
 *
 * @generated s Thu Sep 15 12:09:46 2011
 *
 **/
#define _TYPE  float
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_sgesv_Tile"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_GETRF( n, n ) + FMULS_GETRS( n, nrhs ))
#define _FADDS (FADDS_GETRF( n, n ) + FADDS_GETRS( n, nrhs ))

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    float *AT, *bT, *x;
    float *A = NULL;
    float *b = NULL;
    PLASMA_desc        *descA, *descB;
    real_Double_t       t;
    int                *piv;
    int n     = iparam[TIMING_N];
    int nb    = iparam[TIMING_NB];
    int nrhs  = iparam[TIMING_NRHS];
    int check = iparam[TIMING_CHECK];
    int lda = n;
    int ldb = n;

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
    
    /* Allocate Data */
    AT  = (float *)malloc(lda*n   *sizeof(float));
    bT  = (float *)malloc(ldb*nrhs*sizeof(float));
    piv = (int *)malloc( n*sizeof(int));

    /* Check if unable to allocate memory */
    if ( (!AT) || (!bT) || (!piv) ) {
        printf("Out of Memory \n ");
        return -1;
    }

    /* Initialize AT and bT for Symmetric Positif Matrix */
    PLASMA_Desc_Create(&descA, AT, PlasmaRealFloat, nb, nb, nb*nb, lda, n,    0, 0, n, n);
    PLASMA_Desc_Create(&descB, bT, PlasmaRealFloat, nb, nb, nb*nb, ldb, nrhs, 0, 0, n, nrhs);
    LAPACKE_slarnv_work(1, ISEED, lda*n,    AT);
    LAPACKE_slarnv_work(1, ISEED, ldb*nrhs, bT);

    /* Save AT and bT in lapack layout for check */
    if ( check ) {
        A = (float *)malloc(lda*n   *sizeof(float));
        b = (float *)malloc(ldb*nrhs*sizeof(float));
        PLASMA_sTile_to_Lapack(descA, (void*)A, lda);
        PLASMA_sTile_to_Lapack(descB, (void*)b, ldb);
    }

    t = -cWtime();
    PLASMA_sgesv_Tile( descA, piv, descB );
    t += cWtime();
    *t_ = t;
    
    /* Check the solution */
    if ( check )
      {
        x = (float *)malloc(ldb*nrhs *sizeof(float));
        PLASMA_sTile_to_Lapack(descB, (void*)x, n);

        dparam[TIMING_RES] = s_check_solution(n, n, nrhs, A, lda, b, x, ldb,
                                             &(dparam[TIMING_ANORM]), &(dparam[TIMING_BNORM]), 
                                             &(dparam[TIMING_XNORM]));
        free(A); free(b); free(x);
      }

    PLASMA_Desc_Destroy(&descA);
    PLASMA_Desc_Destroy(&descB);

    free( AT ); free( bT );
    free( piv );
    PLASMA_Finalize();

    return 0;
}
