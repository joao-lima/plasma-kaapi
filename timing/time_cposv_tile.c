/**
 *
 * @generated c Thu Sep 15 12:09:46 2011
 *
 **/
#define _TYPE  PLASMA_Complex32_t
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_cposv_Tile"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_POTRF( n ) + FMULS_POTRS( n, nrhs ))
#define _FADDS (FADDS_POTRF( n ) + FADDS_POTRS( n, nrhs ))

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_)
{
    PLASMA_Complex32_t *A = NULL, *AT, *b = NULL, *bT, *x;
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
    if ( iparam[TIMING_SCHEDULER] )
        PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_DYNAMIC_SCHEDULING );
    else
        PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_STATIC_SCHEDULING );

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
    AT = (PLASMA_Complex32_t *)malloc(nt*nt*nb2*sizeof(PLASMA_Complex32_t));
    bT = (PLASMA_Complex32_t *)malloc(nt*nb2   *sizeof(PLASMA_Complex32_t));

    /* Check if unable to allocate memory */
    if ( (!AT) || (!bT) ) {
        printf("Out of Memory \n ");
        exit(0);
    }

    /* Initialize AT and bT for Symmetric Positif Matrix */
    PLASMA_Desc_Create(&descA, AT, PlasmaComplexFloat, nb, nb, nb*nb, n, n,    0, 0, n, n);
    PLASMA_Desc_Create(&descB, bT, PlasmaComplexFloat, nb, nb, nb*nb, n, nrhs, 0, 0, n, nrhs);
    PLASMA_cplghe_Tile((float)n, descA, 51 );
    LAPACKE_clarnv_work(1, ISEED, nt*nb2, bT);

    /* Save AT and bT in lapack layout for check */
    if ( check ) {
        A = (PLASMA_Complex32_t *)malloc(lda*n    *sizeof(PLASMA_Complex32_t));
        b = (PLASMA_Complex32_t *)malloc(ldb*nrhs *sizeof(PLASMA_Complex32_t));
        PLASMA_Tile_to_Lapack(descA, (void*)A, n);
        PLASMA_Tile_to_Lapack(descB, (void*)b, n);
    }

    /* PLASMA CPOSV */
    t = -cWtime();
    PLASMA_cposv_Tile(PlasmaUpper, descA, descB);
    t += cWtime();
    *t_ = t;

    /* Check the solution */
    if (check)
      {
        x = (PLASMA_Complex32_t *)malloc(ldb*nrhs *sizeof(PLASMA_Complex32_t));
        PLASMA_Tile_to_Lapack(descB, (void*)x, n);

        dparam[TIMING_RES] = c_check_solution(n, n, nrhs, A, lda, b, x, ldb,
                                             &(dparam[TIMING_ANORM]), &(dparam[TIMING_BNORM]), 
                                             &(dparam[TIMING_XNORM]));
        free(A); free(b); free(x);
      }

    PLASMA_Desc_Destroy(&descA);
    PLASMA_Desc_Destroy(&descB);

    free(AT); free(bT);

    PLASMA_Finalize();

    return 0;
}
