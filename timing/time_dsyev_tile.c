/**
 *
 * @generated d Thu Sep 15 12:09:47 2011
 *
 **/
#define _TYPE  double
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_dsyev_Tile"
/* See Lawn 41 page 120 */
#define _FMULS ((2. / 3.) * ((double)n * (double)n * (double)n))
#define _FADDS ((2. / 3.) * ((double)n * (double)n * (double)n))

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    double *AT, *Q = NULL;
    double *W;
    PLASMA_desc *descA = NULL;
    PLASMA_desc *descQ = NULL;
    PLASMA_desc *descT = NULL;
    real_Double_t       t;
    int nb, nb2, nt;
    int n     = iparam[TIMING_N];
    int check = iparam[TIMING_CHECK];
    int lda = n;
    int uplo = PlasmaUpper;
    int vec  = PlasmaNoVec;

    /* Initialize Plasma */ 
    PLASMA_Init( iparam[TIMING_THRDNBR] );
    PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_DYNAMIC_SCHEDULING );

    PLASMA_Disable(PLASMA_AUTOTUNING);
    PLASMA_Set(PLASMA_TILE_SIZE,        iparam[TIMING_NB] );
    PLASMA_Set(PLASMA_INNER_BLOCK_SIZE, iparam[TIMING_IB] );

    nb  = iparam[TIMING_NB];
    nb2 = nb * nb;
    nt  = n / nb + ((n % nb == 0) ? 0 : 1);
    
    /* Allocate Data */
    AT  = (double *)malloc(lda*n*sizeof(double));
    W   = (double *)malloc(n*sizeof(double));
    if (vec == PlasmaVec){
       Q = (double *)malloc(lda*n*sizeof(double));
       if ( (!Q) ) {
          printf("Out of Memory -Q-\n ");
          return -2;
       }
    }
       
    /* Check if unable to allocate memory */
    if ( (!AT) || (!W) ) {
        printf("Out of Memory -\n ");
        return -2;
    }

    /* Initialize Data */
    PLASMA_Desc_Create(&descA, AT, PlasmaRealDouble, nb, nb, nb*nb, lda, n, 0, 0, n, n);
    PLASMA_dplgsy_Tile((double)0.0, descA, 51 );

    if (vec == PlasmaVec)
      PLASMA_Desc_Create(&descQ, Q, PlasmaRealDouble, nb, nb, nb*nb, lda, n, 0, 0, n, n);

    /* Save AT and bT in lapack layout for check */
    if ( check ) {
    }

    /* Allocate Workspace */
    PLASMA_Alloc_Workspace_dsyev(n, n, &descT);

    t = -cWtime();
    PLASMA_dsyev_Tile( vec, uplo, descA, W, descT, descQ );
    t += cWtime();
    *t_ = t;
    
    /* Check the solution */
    if ( check )
      {
      }

    /* DeAllocate Workspace */
    PLASMA_Dealloc_Handle_Tile(&descT);

    PLASMA_Desc_Destroy(&descA);
 
    if (vec == PlasmaVec) {
      PLASMA_Desc_Destroy(&descQ);
      free( Q );
    }
    free( AT );
    free( W );
    PLASMA_Finalize();

    return 0;
}
