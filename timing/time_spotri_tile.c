/**
 *
 * @generated s Thu Sep 15 12:09:46 2011
 *
 **/
#define _TYPE  float
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_spotri_Tile"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_POTRF( n ) + FMULS_POTRI( n ))
#define _FADDS (FADDS_POTRF( n ) + FADDS_POTRI( n ))

//#define POTRI_SYNC

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    float *AT;
    real_Double_t       t;
    PLASMA_desc        *descA;
    int nb, nb2, nt;
    int n     = iparam[TIMING_N];
    int check = iparam[TIMING_CHECK];
    PLASMA_enum uplo = PlasmaLower;

    /* Initialize Plasma */ 
    PLASMA_Init( iparam[TIMING_THRDNBR] );

    /*if ( !iparam[TIMING_AUTOTUNING] ) {*/
        PLASMA_Disable(PLASMA_AUTOTUNING);
        PLASMA_Set(PLASMA_TILE_SIZE,        iparam[TIMING_NB] );
    /* } else { */
    /*     PLASMA_Get(PLASMA_TILE_SIZE,        &iparam[TIMING_NB] ); */
    /* }  */
    nb  = iparam[TIMING_NB];
    nb2 = nb * nb;
    nt  = n / nb + ((n % nb == 0) ? 0 : 1);
    
    /* Allocate Data */
    AT = (float *)malloc(nt*nt*nb2*sizeof(float));

    /* Check if unable to allocate memory */
    if ( !AT ){
        printf("Out of Memory \n ");
        exit(0);
    }

    /* 
     * Initialize Data 
     * It's done in static to avoid having the same sequence than one 
     * the function we want to trace
     */
    PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_STATIC_SCHEDULING );
    PLASMA_Desc_Create(&descA, AT, PlasmaRealFloat, nb, nb, nb*nb, n, n, 0, 0, n, n);
    PLASMA_splgsy_Tile( (float)n, descA, 51 );

    if ( iparam[TIMING_SCHEDULER] )
        PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_DYNAMIC_SCHEDULING );
    else
        PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_STATIC_SCHEDULING );

    /* Save AT in lapack layout for check */
    if ( check ) {
    }

    /* PLASMA SPOTRF / STRTRI / SLAUUM  */
    /*
     * Example of the different way to combine several asynchonous calls
     */
    {
#if defined(TRACE_BY_SEQUENCE)
        PLASMA_sequence *sequence[3];
        PLASMA_request request[3] = { PLASMA_REQUEST_INITIALIZER, 
                                      PLASMA_REQUEST_INITIALIZER, 
                                      PLASMA_REQUEST_INITIALIZER };
        
        PLASMA_Sequence_Create(&sequence[0]);
        PLASMA_Sequence_Create(&sequence[1]);
        PLASMA_Sequence_Create(&sequence[2]);
        
        t = -cWtime();
#if defined(POTRI_SYNC)
        PLASMA_spotrf_Tile_Async(uplo, descA,                sequence[0], &request[0]);
        PLASMA_Sequence_Wait(sequence[0]);
        PLASMA_strtri_Tile_Async(uplo, PlasmaNonUnit, descA, sequence[1], &request[1]);
        PLASMA_Sequence_Wait(sequence[1]);
        PLASMA_slauum_Tile_Async(uplo, descA,                sequence[2], &request[2]);
        PLASMA_Sequence_Wait(sequence[2]);
#else
        PLASMA_spotrf_Tile_Async(uplo, descA,                sequence[0], &request[0]);
        PLASMA_strtri_Tile_Async(uplo, PlasmaNonUnit, descA, sequence[1], &request[1]);
        PLASMA_slauum_Tile_Async(uplo, descA,                sequence[2], &request[2]);
        PLASMA_Sequence_Wait(sequence[0]);
        PLASMA_Sequence_Wait(sequence[1]);
        PLASMA_Sequence_Wait(sequence[2]);
#endif
        t += cWtime();
        
        PLASMA_Sequence_Destroy(sequence[0]);
        PLASMA_Sequence_Destroy(sequence[1]);
        PLASMA_Sequence_Destroy(sequence[2]);
        
#else

#if defined(POTRI_SYNC)
        
        t = -cWtime();
        PLASMA_spotrf_Tile(uplo, descA);
        PLASMA_strtri_Tile(uplo, PlasmaNonUnit, descA);
        PLASMA_slauum_Tile(uplo, descA);
        t += cWtime();
        
#else
        /* Default: we use Asynchonous call with only one sequence */
        
        PLASMA_sequence *sequence;
        PLASMA_request request[2] = { PLASMA_REQUEST_INITIALIZER, 
                                      PLASMA_REQUEST_INITIALIZER };
        
        t = -cWtime();
        PLASMA_Sequence_Create(&sequence);
        PLASMA_spotrf_Tile_Async(uplo, descA, sequence, &request[0]);
        PLASMA_spotri_Tile_Async(uplo, descA, sequence, &request[1]);
        PLASMA_Sequence_Wait(sequence);
        t += cWtime();
        
        PLASMA_Sequence_Destroy(sequence);       
#endif
#endif
        *t_ = t;
    }
    
    /* Check the solution */
    if ( check )
    {
        dparam[TIMING_ANORM] = 0.0;
        dparam[TIMING_XNORM] = 0.0;
        dparam[TIMING_BNORM] = 0.0;
        dparam[TIMING_RES]   = 0.0;
    }

    PLASMA_Desc_Destroy(&descA);
    PLASMA_Finalize();
    free(AT);

    return 0;
}
