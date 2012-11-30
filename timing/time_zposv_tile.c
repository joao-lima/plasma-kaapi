/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_zposv_Tile"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_POTRF( n ) + FMULS_POTRS( n, nrhs ))
#define _FADDS (FADDS_POTRF( n ) + FADDS_POTRS( n, nrhs ))

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_)
{
    PLASMA_Complex64_t *A = NULL, *AT, *b = NULL, *bT, *x;
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
    AT = (PLASMA_Complex64_t *)malloc(nt*nt*nb2*sizeof(PLASMA_Complex64_t));
    bT = (PLASMA_Complex64_t *)malloc(nt*nb2   *sizeof(PLASMA_Complex64_t));

    /* Check if unable to allocate memory */
    if ( (!AT) || (!bT) ) {
        printf("Out of Memory \n ");
        exit(0);
    }

    /* Initialize AT and bT for Symmetric Positif Matrix */
    PLASMA_Desc_Create(&descA, AT, PlasmaComplexDouble, nb, nb, nb*nb, n, n,    0, 0, n, n);
    PLASMA_Desc_Create(&descB, bT, PlasmaComplexDouble, nb, nb, nb*nb, n, nrhs, 0, 0, n, nrhs);
    PLASMA_zplghe_Tile((double)n, descA, 51 );
    LAPACKE_zlarnv_work(1, ISEED, nt*nb2, bT);

    /* Save AT and bT in lapack layout for check */
    if ( check ) {
        A = (PLASMA_Complex64_t *)malloc(lda*n    *sizeof(PLASMA_Complex64_t));
        b = (PLASMA_Complex64_t *)malloc(ldb*nrhs *sizeof(PLASMA_Complex64_t));
        PLASMA_Tile_to_Lapack(descA, (void*)A, n);
        PLASMA_Tile_to_Lapack(descB, (void*)b, n);
    }

    /* PLASMA ZPOSV */
    t = -cWtime();
    PLASMA_zposv_Tile(PlasmaUpper, descA, descB);
    t += cWtime();
    *t_ = t;

    /* Check the solution */
    if (check)
      {
        x = (PLASMA_Complex64_t *)malloc(ldb*nrhs *sizeof(PLASMA_Complex64_t));
        PLASMA_Tile_to_Lapack(descB, (void*)x, n);

        dparam[TIMING_RES] = z_check_solution(n, n, nrhs, A, lda, b, x, ldb,
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
