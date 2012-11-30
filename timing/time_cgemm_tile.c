/**
 *
 * @generated c Thu Sep 15 12:09:46 2011
 *
 **/
#define _TYPE  PLASMA_Complex32_t
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_cgemm_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GEMM(n, n, n)
#define _FADDS FADDS_GEMM(n, n, n)

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    PLASMA_Complex32_t *AT, *BT, *CT;
    PLASMA_Complex32_t *A = NULL, *B = NULL, *C1 = NULL, *C2 = NULL;
    PLASMA_Complex32_t alpha, beta;
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

    AT = (PLASMA_Complex32_t *)malloc(nt*nt*nb2*sizeof(PLASMA_Complex32_t));
    BT = (PLASMA_Complex32_t *)malloc(nt*nt*nb2*sizeof(PLASMA_Complex32_t));
    CT = (PLASMA_Complex32_t *)malloc(nt*nt*nb2*sizeof(PLASMA_Complex32_t));

    /* Check if unable to allocate memory */
    if ( (!AT) || (!BT) || (!CT) ) {
        printf("Out of Memory \n ");
        exit(0);
    }
    
     /* Initialiaze Data */
    LAPACKE_clarnv_work(1, ISEED, 1, &alpha);
    LAPACKE_clarnv_work(1, ISEED, 1, &beta);
    LAPACKE_clarnv_work(1, ISEED, nt*nt*nb2, AT);
    LAPACKE_clarnv_work(1, ISEED, nt*nt*nb2, BT);
    LAPACKE_clarnv_work(1, ISEED, nt*nt*nb2, CT);

    /* Initialize AT and bT for Symmetric Positif Matrix */
    PLASMA_Desc_Create(&descA, AT, PlasmaComplexFloat, nb, nb, nb*nb, n, n, 0, 0, n, n);
    PLASMA_Desc_Create(&descB, BT, PlasmaComplexFloat, nb, nb, nb*nb, n, n, 0, 0, n, n);
    PLASMA_Desc_Create(&descC, CT, PlasmaComplexFloat, nb, nb, nb*nb, n, n, 0, 0, n, n);

    if (check)
      {
          C2 = (PLASMA_Complex32_t *)malloc(n*lda*sizeof(PLASMA_Complex32_t));
          PLASMA_Tile_to_Lapack(descC, (void*)C2, n);
      }

    t = -cWtime();
    PLASMA_cgemm_Tile( PlasmaNoTrans, PlasmaNoTrans, alpha, descA, descB, beta, descC );
    t += cWtime();
    *t_ = t;
    
    /* Check the solution */
    if (check)
      {
          A = (PLASMA_Complex32_t *)malloc(n*lda*sizeof(PLASMA_Complex32_t));
          PLASMA_Tile_to_Lapack(descA, (void*)A, n);
          free(AT);

          B = (PLASMA_Complex32_t *)malloc(n*lda*sizeof(PLASMA_Complex32_t));
          PLASMA_Tile_to_Lapack(descB, (void*)B, n);
          free(BT);

          C1 = (PLASMA_Complex32_t *)malloc(n*lda*sizeof(PLASMA_Complex32_t));
          PLASMA_Tile_to_Lapack(descC, (void*)C1, n);
          free(CT);

          dparam[TIMING_RES] = c_check_gemm( PlasmaNoTrans, PlasmaNoTrans, n, n, n, 
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
