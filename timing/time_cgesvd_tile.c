/**
 *
 * @generated c Thu Sep 15 12:09:47 2011
 *
 **/
#define _TYPE  PLASMA_Complex32_t
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_cheev_Tile"
/* See Lawn 41 page 120 */
#define _FMULS ((2. / 3.) * ((float)n * (float)n * (float)n))
#define _FADDS ((2. / 3.) * ((float)n * (float)n * (float)n))

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    PLASMA_Complex32_t *AT, *U = NULL, *VT = NULL;
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
    AT  = (PLASMA_Complex32_t *)malloc(nt*nt*nb2*sizeof(PLASMA_Complex32_t));
    S  = (float *)malloc(n*sizeof(float));
    if (jobu == PlasmaVec){
       U = (PLASMA_Complex32_t *)malloc(nt*nt*nb2*sizeof(PLASMA_Complex32_t));
       if ( (!U) ) {
          printf("Out of Memory -U-\n ");
          exit(0);
       }
    }
    if (jobvt == PlasmaVec){
       VT = (PLASMA_Complex32_t *)malloc(nt*nt*nb2*sizeof(PLASMA_Complex32_t));
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
    PLASMA_Desc_Create(&descA, AT, PlasmaComplexFloat, nb, nb, nb*nb, lda, n, 0, 0, n, n);
    PLASMA_cplrnt_Tile(descA, 51 );
    if (jobu == PlasmaVec)
       PLASMA_Desc_Create(&descU,  U , PlasmaComplexFloat, nb, nb, nb*nb, lda, n, 0, 0, n, n);
    if (jobvt == PlasmaVec)
       PLASMA_Desc_Create(&descVT, VT, PlasmaComplexFloat, nb, nb, nb*nb, lda, n, 0, 0, n, n);

    /* Save AT and bT in lapack layout for check */
    if ( check ) {
    }

    /* Allocate Workspace */
    PLASMA_Alloc_Workspace_cgesvd(n, n, &descT);

    t = -cWtime();
    PLASMA_cgesvd_Tile(jobu, jobvt, descA, S, descU, descVT, descT);
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
