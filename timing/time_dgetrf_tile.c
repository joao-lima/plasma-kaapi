/**
 *
 * @generated d Thu Sep 15 12:09:46 2011
 *
 **/
#define _TYPE  double
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_dgetrf_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GETRF(m, n)
#define _FADDS FADDS_GETRF(m, n)

#include "./timing.c"

#if defined(CONFIG_USE_CUDA)
#include "core_cublas.h"
#endif

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    double *A = NULL, *AT, *b = NULL, *bT, *x;
    PLASMA_desc        *descA, *descB;
    real_Double_t       t;
    int                *piv;
    int m     = iparam[TIMING_M];
    int n     = iparam[TIMING_N];
    int nb    = iparam[TIMING_NB];
    int nrhs  = iparam[TIMING_NRHS];
    int check = iparam[TIMING_CHECK];
    int lda   = m;
    int ldb   = m;

    /* Allocate Data */
    AT  = (double *)malloc(lda*n*sizeof(double));
    piv = (int *)malloc( min(m, n) * sizeof(int));

    /* Check if unable to allocate memory */
    if ( !AT || !piv ){
        printf("Out of Memory \n ");
        return -1;
    }
    
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
    
    /* Initialiaze Data */
    PLASMA_Desc_Create(&descA, AT, PlasmaRealDouble, nb, nb, nb*nb, lda, n, 0, 0, m, n);
    /*LAPACKE_dlarnv_work(1, ISEED, n*lda, AT);*/
    PLASMA_dplrnt_Tile(descA, 3456);

    /* Save AT in lapack layout for check */
    if ( check && (m == n) ) {
        A = (double *)malloc(lda*n    *sizeof(double));
        PLASMA_dTile_to_Lapack(descA, (void*)A, lda);
    }

#if defined(CONFIG_USE_CUDA)
    cudaHostRegister(AT, lda*n*sizeof(double), cudaHostRegisterPortable);
    cudaHostRegister(piv, min(m,n)*sizeof(int), cudaHostRegisterPortable);
#endif

#if defined(_CORE_CUBLAS_H_)
    core_cublas_init();
#endif

    t = -cWtime();
    PLASMA_dgetrf_Tile( descA, piv );
    t += cWtime();
    *t_ = t;
    
    /* Check the solution */
    if ( check && (m == n) )
      {
        b  = (double *)malloc(ldb*nrhs*sizeof(double));
        bT = (double *)malloc(ldb*nrhs*sizeof(double));
        x  = (double *)malloc(ldb*nrhs*sizeof(double));

        LAPACKE_dlarnv_work(1, ISEED, ldb*nrhs, b);

        PLASMA_Desc_Create(&descB, bT, PlasmaRealDouble, nb, nb, nb*nb, ldb, nrhs, 0, 0, n, nrhs);
        PLASMA_dLapack_to_Tile((void*)b, ldb, descB);

        PLASMA_dgetrs_Tile( PlasmaNoTrans, descA, piv, descB );

        PLASMA_dTile_to_Lapack(descB, (void*)x, ldb);
        PLASMA_Desc_Destroy(&descB);

        dparam[TIMING_RES] = d_check_solution(m, n, nrhs, A, lda, b, x, ldb,
                                              &(dparam[TIMING_ANORM]), &(dparam[TIMING_BNORM]), 
                                              &(dparam[TIMING_XNORM]));
        
        free( A ); free( b ); free( bT ); free( x );
      }

    PLASMA_Desc_Destroy(&descA);

    PLASMA_Finalize();
#if defined(CONFIG_USE_CUDA)
    cudaHostUnregister(AT);
    cudaHostUnregister(piv);
#endif
    free( AT );
    free( piv );

    return 0;
}
