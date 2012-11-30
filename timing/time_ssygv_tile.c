/**
 *
 * @generated s Thu Sep 15 12:09:47 2011
 *
 **/
#define _TYPE  float
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_ssygv_Tile"
/* See Lawn 41 page 120 */
/* cholesky + 2 trsm's + trd */
#define _FMULS ((2. / 3.) * ((float)n * (float)n * (float)n)) + (n * (1.0 / 6.0 * n + 0.5) * n) + 2 * (n * n * (float)((n + 1) / 2.0) )
#define _FADDS ((2. / 3.) * ((float)n * (float)n * (float)n)) + (n * (1.0 / 6.0 * n )      * n) + 2 * (n * n * (float)((n + 1) / 2.0) )

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    float *AT, *BT, *Q = NULL;
    float *W;
    PLASMA_desc *descA, *descB, *descQ, *descT;
    real_Double_t       t;
    int nb, nb2, nt;
    int n     = iparam[TIMING_N];
    int check = iparam[TIMING_CHECK];
    int lda = n;
    int itype = 1;
    int vec   = PlasmaNoVec;
    int uplo  = PlasmaUpper;

    /* Initialize Plasma */ 
    PLASMA_Init( iparam[TIMING_THRDNBR] );
    //    if ( iparam[TIMING_SCHEDULER] )
        PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_DYNAMIC_SCHEDULING );
    /* else */
    /*     PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_STATIC_SCHEDULING ); */

    /*if ( !iparam[TIMING_AUTOTUNING] ) {*/
        PLASMA_Disable(PLASMA_AUTOTUNING);
        PLASMA_Set(PLASMA_TILE_SIZE,        iparam[TIMING_NB] );
        PLASMA_Set(PLASMA_INNER_BLOCK_SIZE, iparam[TIMING_IB] );
    /* } else { */
    /*     PLASMA_Get(PLASMA_TILE_SIZE,        &iparam[TIMING_NB] ); */
    /*     PLASMA_Get(PLASMA_INNER_BLOCK_SIZE, &iparam[TIMING_IB] ); */
    /* }  */
    nb  = iparam[TIMING_NB];
    nb2 = nb * nb;
    nt  = n / nb + ((n % nb == 0) ? 0 : 1);
    
    /* Allocate Data */
    AT  = (float *)malloc(nt*nt*nb2*sizeof(float));
    BT  = (float *)malloc(nt*nt*nb2*sizeof(float));
    W  = (float *)malloc(n*sizeof(float));
    if (vec == PlasmaVec){
       Q = (float *)malloc(nt*nt*nb2*sizeof(float));
       if ( (!Q) ) {
          printf("Out of Memory -Q-\n ");
          exit(0);
       }
    }
       
    /* Check if unable to allocate memory */
    if ( (!AT) || (!BT) || (!W) ) {
        printf("Out of Memory -\n ");
        exit(0);
    }

    /* Initialiaze Data */
    PLASMA_Desc_Create(&descA, AT, PlasmaRealFloat, nb, nb, nb*nb, lda, n, 0, 0, n, n);
    PLASMA_splgsy_Tile((float)0.0, descA, 51 );

    PLASMA_Desc_Create(&descB, BT, PlasmaRealFloat, nb, nb, nb*nb, lda, n, 0, 0, n, n);
    PLASMA_splgsy_Tile((float)n, descB, 51 );

    PLASMA_Desc_Create(&descQ, Q, PlasmaRealFloat, nb, nb, nb*nb, lda, n, 0, 0, n, n);

    /* Save AT and bT in lapack layout for check */
    if ( check ) {
    }

    /* Allocate Workspace */
    PLASMA_Alloc_Workspace_ssygv(n, n, &descT);

    t = -cWtime();
    PLASMA_ssygv_Tile( itype, vec, uplo, descA, descB, W, descT, descQ );
    t += cWtime();
    *t_ = t;
    
    /* Check the solution */
    if ( check )
      {
      }

    /* DeAllocate Workspace */
    PLASMA_Dealloc_Handle_Tile(&descT);

    PLASMA_Desc_Destroy(&descA);
    PLASMA_Desc_Destroy(&descB);
    PLASMA_Desc_Destroy(&descQ);

    if (vec == PlasmaVec)
       free( Q );
    free( AT );
    free( W );
    PLASMA_Finalize();

    return 0;
}
