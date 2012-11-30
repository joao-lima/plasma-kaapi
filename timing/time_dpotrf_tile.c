/**
 *
 * @generated d Thu Sep 15 12:09:46 2011
 *
 **/
#define _TYPE  double
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_dpotrf_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_POTRF( n )
#define _FADDS FADDS_POTRF( n )

#include "./timing.c"

#include "core_cublas.h"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    double *A = NULL, *AT, *b = NULL, *bT, *x;
    real_Double_t       t;
    PLASMA_desc        *descA, *descB;
    int nb, nb2, nt;
    int n     = iparam[TIMING_N];
    int nrhs  = iparam[TIMING_NRHS];
    int check = iparam[TIMING_CHECK];
    int lda = n;
    int ldb = n;

    /* Initialize Plasma */ 
    PLASMA_Init( iparam[TIMING_THRDNBR] );
    PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_DYNAMIC_SCHEDULING );

#if defined(_CORE_CUBLAS_H_)
    core_cublas_init();
#endif
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
//    AT = (double *)malloc(nt*nt*nb2*sizeof(double));
    cudaHostAlloc((void**)&AT, nt*nt*nb2*sizeof(double), cudaHostAllocPortable);

    /* Check if unable to allocate memory */
    if ( !AT ){
        printf("Out of Memory \n ");
        exit(0);
    }

//    cudaHostRegister(AT, nt*nt*nb2*sizeof(double), cudaHostRegisterPortable);

    /* Initialiaze Data */
    PLASMA_Desc_Create(&descA, AT, PlasmaRealDouble, nb, nb, nb*nb, n, n, 0, 0, n, n);
    PLASMA_dplgsy_Tile((double)n, descA, 51 );

    /* Save AT in lapack layout for check */
    if ( check ) {
        A = (double *)malloc(lda*n    *sizeof(double));
        PLASMA_Tile_to_Lapack(descA, (void*)A, n);
    }

    /* PLASMA DPOSV */
    t = -cWtime();
    PLASMA_dpotrf_Tile(PlasmaUpper, descA);
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

        PLASMA_dpotrs_Tile( PlasmaUpper, descA, descB );

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

    PLASMA_Desc_Destroy(&descA);
    PLASMA_Finalize();

    cudaFreeHost(AT);
    kaapi_finalize();
//    cudaHostUnregister(AT);
 //   free(AT);

    return 0;
}
