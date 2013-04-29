/**
 *
 * @generated d Thu Sep 15 12:09:47 2011
 *
 **/
#define _TYPE  double
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_dgetrf_incpiv_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GETRF(n, n)
#define _FADDS FADDS_GETRF(n, n)

#include "./timing.c"

#if defined(CONFIG_USE_CUDA)
#include "core_cublas.h"
#endif

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    double *A = NULL, *AT, *b, *bT, *x;
    PLASMA_desc        *descA, *descB, *descL;
    real_Double_t       t;
    int                *piv;
    int nb, nb2, nt;
    int n     = iparam[TIMING_N];
    int nrhs  = iparam[TIMING_NRHS];
    int check = iparam[TIMING_CHECK];
    int lda      = n;
    int ldb      = n;

    /* Initialize Plasma */ 
    PLASMA_Init( iparam[TIMING_THRDNBR] );
    if ( iparam[TIMING_SCHEDULER] )
        PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_DYNAMIC_SCHEDULING );
    else
        PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_STATIC_SCHEDULING );

#if defined(_CORE_CUBLAS_H_)
    core_cublas_init();
#endif

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

#if defined(CONFIG_USE_CUDA)
    cudaHostRegister(AT, nt*nt*nb2*sizeof(double), cudaHostRegisterPortable);
#endif

    /* Initialiaze Data */
    PLASMA_Desc_Create(&descA, AT, PlasmaRealDouble, nb, nb, nb*nb, n, n, 0, 0, n, n);
    LAPACKE_dlarnv_work(1, ISEED, nt*nt*nb2, AT);

    /* Allocate Workspace */
    PLASMA_Alloc_Workspace_dgesv_incpiv_Tile(n, &descL, &piv);

    {
      int NB, MT, NT;
      size_t size;
      NB = nb;
      NT = (n%NB==0) ? (n/NB) : ((n/NB)+1);
      MT = (n%NB==0) ? (n/NB) : ((n/NB)+1);
      size = (size_t)MT*NT*NB * sizeof(int);
#if defined(CONFIG_USE_CUDA)
      cudaHostRegister((void*)piv, size, cudaHostRegisterPortable);
#endif
    }
#if defined(CONFIG_USE_CUDA)
    cudaHostRegister((void*)descL->mat, descL->lm*descL->ln*sizeof(double), cudaHostRegisterPortable);
#endif

    /* Save AT in lapack layout for check */
    if ( check ) {
        A = (double *)malloc(lda*n    *sizeof(double));
        PLASMA_Tile_to_Lapack(descA, (void*)A, n);
    }

    t = -cWtime();
    PLASMA_dgetrf_incpiv_Tile( descA, descL, piv );
    t += cWtime();
    *t_ = t;
    
    /* Check the solution */
    if ( check )
      {
        b  = (double *)malloc(ldb*nrhs *sizeof(double));
        bT = (double *)malloc(nt*nb2   *sizeof(double));
        x  = (double *)malloc(ldb*nrhs *sizeof(double));

        LAPACKE_dlarnv_work(1, ISEED, n*nrhs, b);
        PLASMA_Desc_Create(&descB, bT, PlasmaRealDouble, nb, nb, nb*nb, n, nrhs, 0, 0, n, nrhs);
        PLASMA_Lapack_to_Tile((void*)b, n, descB);

        PLASMA_dgetrs_incpiv_Tile( descA, descL, piv, descB );

        PLASMA_Tile_to_Lapack(descB, (void*)x, n);

        dparam[TIMING_RES] = d_check_solution(n, n, nrhs, A, lda, b, x, ldb,
                                             &(dparam[TIMING_ANORM]), &(dparam[TIMING_BNORM]), 
                                             &(dparam[TIMING_XNORM]));

        PLASMA_Desc_Destroy(&descB);
        free( A ); free( b ); free( bT ); free( x );
      }

    /* Deallocate Workspace */
    PLASMA_Dealloc_Handle_Tile(&descL);

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
