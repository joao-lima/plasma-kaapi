/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_zgebrd_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GEBRD( n, n )
#define _FADDS FADDS_GEBRD( n, n )

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    PLASMA_Complex64_t *AT, *Q = NULL;
    double *W;
    PLASMA_desc *descA, *descQ, *descT;
    real_Double_t       t;
    int nb, nb2, nt;
    int n     = iparam[TIMING_N];
    int check = iparam[TIMING_CHECK];
    int lda = n;
    int vec  = PlasmaNoVec;

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
    AT  = (PLASMA_Complex64_t *)malloc(nt*nt*nb2*sizeof(PLASMA_Complex64_t));
    W  = (double *)malloc(n*sizeof(double));
    if (vec == PlasmaVec)
       Q = (PLASMA_Complex64_t *)malloc(nt*nt*nb2*sizeof(PLASMA_Complex64_t));

    /* Check if unable to allocate memory */
    if ( (!AT) || (!Q) || (!W) ) {
        printf("Out of Memory \n ");
        exit(0);
    }

    /* Initialiaze Data */
    PLASMA_Desc_Create(&descA, AT, PlasmaComplexDouble, nb, nb, nb*nb, lda, n, 0, 0, n, n);
    PLASMA_zplghe_Tile((double)0.0, descA, 51 );
    PLASMA_Desc_Create(&descQ, Q, PlasmaComplexDouble, nb, nb, nb*nb, lda, n, 0, 0, n, n);

    /* Save AT and bT in lapack layout for check */
    if ( check ) {
    }

    /* Allocate Workspace */
    PLASMA_Alloc_Workspace_zgesvd(n, n, &descT);

    t = -cWtime();
    PLASMA_zgesvd_Tile( vec, vec, descA, W, descQ, descQ, descT );
    t += cWtime();
    *t_ = t;
    
    /* Check the solution */
    if ( check )
      {
      }

    /* DeAllocate Workspace */
    PLASMA_Dealloc_Handle_Tile(&descT);

    PLASMA_Desc_Destroy(&descA);

    if (vec == PlasmaVec)
       free( Q );
    free( AT );
    free( W );
    PLASMA_Finalize();

    return 0;
}
