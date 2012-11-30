/**
 *
 * @generated s Thu Sep 15 12:09:48 2011
 *
 **/
#define _TYPE  float
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_ssyev_Tile"
/* See Lawn 41 page 120 */
#define _FMULS ((2. / 3.) * ((float)n * (float)n * (float)n))
#define _FADDS ((2. / 3.) * ((float)n * (float)n * (float)n))

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    float *AT, *U = NULL, *VT = NULL;
    float *S;
    PLASMA_desc *descA  = NULL;
    PLASMA_desc *descU  = NULL;
    PLASMA_desc *descVT = NULL;
    PLASMA_desc *descT  = NULL;
    real_Double_t t;
    int nb, nb2, nt;
    int n     = iparam[TIMING_N];
    int check = iparam[TIMING_CHECK];
    int lda   = n;
    int jobu  = PlasmaNoVec;
    int jobvt = PlasmaNoVec;

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
    S  = (float *)malloc(n*sizeof(float));
    if (jobu == PlasmaVec){
       U = (float *)malloc(nt*nt*nb2*sizeof(float));
       if ( (!U) ) {
          printf("Out of Memory -U-\n ");
          exit(0);
       }
    }
    if (jobvt == PlasmaVec){
       VT = (float *)malloc(nt*nt*nb2*sizeof(float));
       if ( (!VT) ) {
          printf("Out of Memory -VT-\n ");
          exit(0);
       }
    }
       
    /* Check if unable to allocate memory */
    if ( (!AT) || (!S) ) {
        printf("Out of Memory -\n ");
        exit(0);
    }

    /* Initialiaze Data */
    PLASMA_Desc_Create(&descA, AT, PlasmaRealFloat, nb, nb, nb*nb, lda, n, 0, 0, n, n);
    PLASMA_splrnt_Tile(descA, 51 );
    if (jobu == PlasmaVec)
       PLASMA_Desc_Create(&descU,  U , PlasmaRealFloat, nb, nb, nb*nb, lda, n, 0, 0, n, n);
    if (jobvt == PlasmaVec)
       PLASMA_Desc_Create(&descVT, VT, PlasmaRealFloat, nb, nb, nb*nb, lda, n, 0, 0, n, n);

    /* Save AT and bT in lapack layout for check */
    if ( check ) {
    }

    /* Allocate Workspace */
    PLASMA_Alloc_Workspace_sgesvd(n, n, &descT);

    t = -cWtime();
    PLASMA_sgesvd_Tile(jobu, jobvt, descA, S, descU, descVT, descT);
    t += cWtime();
    *t_ = t;
    
    /* Check the solution */
    if ( check )
      {
      }

    /* DeAllocate Workspace */
    PLASMA_Dealloc_Handle_Tile(&descT);

    PLASMA_Desc_Destroy(&descA);

    if (jobu == PlasmaVec)
       free( U );
    if (jobvt == PlasmaVec)
       free( VT );
    free( AT );
    free( S );
    PLASMA_Finalize();

    return 0;
}
