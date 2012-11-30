/**
 *
 * @generated ds Thu Sep 15 12:09:46 2011
 *
 **/
#define _TYPE  double
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_dgesv_Tile"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_GETRF( n, n ) + FMULS_GETRS( n, nrhs ))
#define _FADDS (FADDS_GETRF( n, n ) + FADDS_GETRS( n, nrhs ))

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    double *A = NULL, *AT, *b = NULL, *bT, *x, *xT;
    PLASMA_desc        *descA, *descB, *descX;
    real_Double_t       t;
    int                *piv;
    int n     = iparam[TIMING_N];
    int nb    = iparam[TIMING_NB];
    int nrhs  = iparam[TIMING_NRHS];
    int check = iparam[TIMING_CHECK];
    int lda = n;
    int ldb = n;
    int iter;

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
    AT  = (double *)malloc(lda*n   *sizeof(double));
    bT  = (double *)malloc(ldb*nrhs*sizeof(double));
    xT  = (double *)malloc(ldb*nrhs*sizeof(double));
    piv = (int *)malloc( n*sizeof(int));

    /* Check if unable to allocate memory */
    if ( (!AT) || (!bT) || (!xT) || (!piv) ) {
        printf("Out of Memory \n ");
        return -1;
    }

    /* Initialize AT and bT for Symmetric Positif Matrix */
    PLASMA_Desc_Create(&descA, AT, PlasmaRealDouble, nb, nb, nb*nb, lda, n,    0, 0, n, n);
    PLASMA_Desc_Create(&descB, bT, PlasmaRealDouble, nb, nb, nb*nb, ldb, nrhs, 0, 0, n, nrhs);
    PLASMA_Desc_Create(&descX, xT, PlasmaRealDouble, nb, nb, nb*nb, ldb, nrhs, 0, 0, n, nrhs);
    LAPACKE_dlarnv_work(1, ISEED, lda*n,    AT);
    LAPACKE_dlarnv_work(1, ISEED, ldb*nrhs, bT);

    /* Save AT and bT in lapack layout for check */
    if ( check ) {
        A = (double *)malloc(lda*n   *sizeof(double));
        b = (double *)malloc(ldb*nrhs*sizeof(double));
        PLASMA_zTile_to_Lapack(descA, (void*)A, lda);
        PLASMA_zTile_to_Lapack(descB, (void*)b, ldb);
    }

    t = -cWtime();
    PLASMA_dsgesv_Tile( descA, piv, descB, descX, &iter );
    t += cWtime();
    *t_ = t;
    
    /* Check the solution */
    if ( check )
      {
        x = (double *)malloc(ldb*nrhs *sizeof(double));
        PLASMA_zTile_to_Lapack(descX, (void*)x, n);

        dparam[TIMING_RES] = d_check_solution(n, n, nrhs, A, lda, b, x, ldb,
                                             &(dparam[TIMING_ANORM]), &(dparam[TIMING_BNORM]), 
                                             &(dparam[TIMING_XNORM]));
        free(A); free(b); free(x);
      }

    PLASMA_Desc_Destroy(&descA);
    PLASMA_Desc_Destroy(&descB);

    free( AT ); free( bT ); free(xT);
    free( piv );
    PLASMA_Finalize();

    return 0;
}
