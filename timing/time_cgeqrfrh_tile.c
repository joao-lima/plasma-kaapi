/**
 *
 * @generated c Thu Sep 15 12:09:47 2011
 *
 **/
#define _TYPE  PLASMA_Complex32_t
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_cgeqrfrh_Tile"
/* See Lawn 41 page 120 */
/* Matrix is n x m in timing.c */
#define _FMULS (m * (m * (n - 1.0/3*m + 0.5 ) + n + 23.0/6) )
#define _FADDS (m * (m * (n - 1.0/3*m + 0.5 ) + 5.0/6))

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    PLASMA_Complex32_t *AT;
    PLASMA_desc        *descA, *descT;
    real_Double_t       t;
    int nb;
    int M     = iparam[TIMING_N];
    int N     = iparam[TIMING_M];
    //int N = M/nrhs; //RUN WITH NRHS = 10 or 20 (ALSO USED IN TIMING.C)
    int lda = M;

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
    nb  = iparam[TIMING_NB];

    /* Householder mode */
    //PLASMA_Set(PLASMA_HOUSEHOLDER_MODE, PLASMA_FLAT_HOUSEHOLDER);
    PLASMA_Set(PLASMA_HOUSEHOLDER_MODE, PLASMA_TREE_HOUSEHOLDER);
    PLASMA_Set(PLASMA_HOUSEHOLDER_SIZE, 4);
    
    /* Allocate Data */
    AT  = (PLASMA_Complex32_t *)malloc(lda*N*sizeof(PLASMA_Complex32_t));

    /* Check if unable to allocate memory */
    if ( !AT ){
        printf("Out of Memory \n ");
        exit(0);
    }

    /* Initialiaze Data */
    PLASMA_Desc_Create(&descA, AT, PlasmaComplexFloat, nb, nb, nb*nb, M, N, 0, 0, M, N);
    LAPACKE_clarnv_work(1, ISEED, lda*N, AT);

    /* Allocate Workspace */
    PLASMA_Alloc_Workspace_cgels_Tile(M, N, &descT);

    t = -cWtime();
    PLASMA_cgeqrf_Tile( descA, descT );
    t += cWtime();
    *t_ = t;
    
    /* Allocate Workspace */
    PLASMA_Dealloc_Handle_Tile(&descT);

    PLASMA_Desc_Destroy(&descA);

    free( AT );
    PLASMA_Finalize();

    return 0;
}
