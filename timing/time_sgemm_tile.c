/**
 *
 * @generated s Thu Sep 15 12:09:46 2011
 *
 **/
#define _TYPE  float
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_sgemm_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GEMM(n, n, n)
#define _FADDS FADDS_GEMM(n, n, n)

#include "./timing.c"

#include "core_cublas.h"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    float *AT, *BT, *CT;
    float *A = NULL, *B = NULL, *C1 = NULL, *C2 = NULL;
    float alpha, beta;
    PLASMA_desc        *descA, *descB, *descC;
    real_Double_t       t;
    int nb, nb2, nt;
    int n       = iparam[TIMING_N];
    int check   = iparam[TIMING_CHECK];
    int lda     = n;
    
    /* Allocate Data */
    /* Initialize Plasma */ 
    PLASMA_Init( iparam[TIMING_THRDNBR] );
    if ( iparam[TIMING_SCHEDULER] )
        PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_DYNAMIC_SCHEDULING );
    else
        PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_STATIC_SCHEDULING );

    /*if ( !iparam[TIMING_AUTOTUNING] ) {*/
        PLASMA_Disable(PLASMA_AUTOTUNING);
        PLASMA_Set(PLASMA_TILE_SIZE,        iparam[TIMING_NB] );
    /* } */
    /* } else { */
    /*     PLASMA_Get(PLASMA_TILE_SIZE,        &iparam[TIMING_NB] ); */
    /* }  */
    nb  = iparam[TIMING_NB];
    nb2 = nb * nb;
    nt  = n / nb + ((n % nb == 0) ? 0 : 1);

    AT = (float *)malloc(nt*nt*nb2*sizeof(float));
    BT = (float *)malloc(nt*nt*nb2*sizeof(float));
    CT = (float *)malloc(nt*nt*nb2*sizeof(float));

    /* Check if unable to allocate memory */
    if ( (!AT) || (!BT) || (!CT) ) {
        printf("Out of Memory \n ");
        exit(0);
    }
    
    cudaHostRegister(AT, nt*nt*nb2*sizeof(float), cudaHostRegisterPortable);
    cudaHostRegister(BT, nt*nt*nb2*sizeof(float), cudaHostRegisterPortable);
    cudaHostRegister(CT, nt*nt*nb2*sizeof(float), cudaHostRegisterPortable);

     /* Initialiaze Data */
    LAPACKE_slarnv_work(1, ISEED, 1, &alpha);
    LAPACKE_slarnv_work(1, ISEED, 1, &beta);
    LAPACKE_slarnv_work(1, ISEED, nt*nt*nb2, AT);
    LAPACKE_slarnv_work(1, ISEED, nt*nt*nb2, BT);
    LAPACKE_slarnv_work(1, ISEED, nt*nt*nb2, CT);

    /* Initialize AT and bT for Symmetric Positif Matrix */
    PLASMA_Desc_Create(&descA, AT, PlasmaRealFloat, nb, nb, nb*nb, n, n, 0, 0, n, n);
    PLASMA_Desc_Create(&descB, BT, PlasmaRealFloat, nb, nb, nb*nb, n, n, 0, 0, n, n);
    PLASMA_Desc_Create(&descC, CT, PlasmaRealFloat, nb, nb, nb*nb, n, n, 0, 0, n, n);

    if (check)
      {
          C2 = (float *)malloc(n*lda*sizeof(float));
          PLASMA_Tile_to_Lapack(descC, (void*)C2, n);
      }

#if defined(_CORE_CUBLAS_H_)
    core_cublas_init();
#endif

    t = -cWtime();
    PLASMA_sgemm_Tile( PlasmaNoTrans, PlasmaNoTrans, alpha, descA, descB, beta, descC );
    t += cWtime();
    *t_ = t;
    
    /* Check the solution */
    if (check)
      {
          A = (float *)malloc(n*lda*sizeof(float));
          PLASMA_Tile_to_Lapack(descA, (void*)A, n);
          free(AT);

          B = (float *)malloc(n*lda*sizeof(float));
          PLASMA_Tile_to_Lapack(descB, (void*)B, n);
          free(BT);

          C1 = (float *)malloc(n*lda*sizeof(float));
          PLASMA_Tile_to_Lapack(descC, (void*)C1, n);
          free(CT);

          dparam[TIMING_RES] = s_check_gemm( PlasmaNoTrans, PlasmaNoTrans, n, n, n, 
                                            alpha, A, lda, B, lda, beta, C1, C2, lda,
                                            &(dparam[TIMING_ANORM]), &(dparam[TIMING_BNORM]), 
                                            &(dparam[TIMING_XNORM]));
          free(C2);
      }
    else {
        free( AT );
        free( BT );
        free( CT );
    }

    PLASMA_Desc_Destroy(&descA);
    PLASMA_Desc_Destroy(&descB);
    PLASMA_Desc_Destroy(&descC);
    PLASMA_Finalize();

    return 0;
}
