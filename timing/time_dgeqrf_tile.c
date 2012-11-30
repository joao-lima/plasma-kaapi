/**
 *
 * @generated d Thu Sep 15 12:09:47 2011
 *
 **/
#define _TYPE  double
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_dgeqrf_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GEQRF(n, n)
#define _FADDS FADDS_GEQRF(n, n)

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    double *A = NULL, *AT, *b = NULL, *bT, *x;
    PLASMA_desc        *descA, *descB, *descT;
    real_Double_t       t;
    int nb, nb2, nt;
    int n     = iparam[TIMING_N];
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
    nb  = iparam[TIMING_NB];
    nb2 = nb * nb;
    nt  = n / nb + ((n % nb == 0) ? 0 : 1);
    
    /* Allocate Data */
    AT  = (double *)malloc(nt*nt*nb2*sizeof(double));

    /* Check if unable to allocate memory */
    if ( !AT ){
        printf("Out of Memory \n ");
        exit(0);
    }

    /* Initialiaze Data */
    PLASMA_Desc_Create(&descA, AT, PlasmaRealDouble, nb, nb, nb*nb, n, n, 0, 0, n, n);
    LAPACKE_dlarnv_work(1, ISEED, nt*nt*nb2, AT);

    /* Allocate Workspace */
    PLASMA_Alloc_Workspace_dgels_Tile(n, n, &descT);

    /* Save AT in lapack layout for check */
    if ( check ) {
        A = (double *)malloc(lda*n    *sizeof(double));
        PLASMA_Tile_to_Lapack(descA, (void*)A, n);
    }

    t = -cWtime();
    PLASMA_dgeqrf_Tile( descA, descT );
    t += cWtime();
    *t_ = t;
    
    /* Check the solution */
    if ( check )
      {
        b  = (double *)malloc(ldb*nrhs *sizeof(double));
        bT = (double *)malloc(nt*nb2   *sizeof(double));
        x  = (double *)malloc(ldb*nrhs *sizeof(double));

        LAPACKE_dlarnv_work(1, ISEED, nt*nb2, bT);
        PLASMA_Desc_Create(&descB, bT, PlasmaRealDouble, nb, nb, nb*nb, n, nrhs, 0, 0, n, nrhs);
        PLASMA_Tile_to_Lapack(descB, (void*)b, n);

        PLASMA_dgeqrs_Tile( descA, descT, descB );

        PLASMA_Tile_to_Lapack(descB, (void*)x, n);

        dparam[TIMING_RES] = d_check_solution(n, n, nrhs, A, lda, b, x, ldb,
                                             &(dparam[TIMING_ANORM]), &(dparam[TIMING_BNORM]), 
                                             &(dparam[TIMING_XNORM]));

        PLASMA_Desc_Destroy(&descB);
        free( A ); 
        free( b ); 
        free( bT ); 
        free( x );
      }

    /* Allocate Workspace */
    PLASMA_Dealloc_Handle_Tile(&descT);

    PLASMA_Desc_Destroy(&descA);

    free( AT );
    PLASMA_Finalize();

    return 0;
}
