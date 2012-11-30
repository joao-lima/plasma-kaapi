#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH lapack_dlamch

#define _NAME  "PLASMA_ztrsm"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_TRSM( PlasmaLeft, n, n )
#define _FADDS FADDS_TRSM( PlasmaLeft, n, n )

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_)
{
    PLASMA_Complex64_t *A, *B1, *B2 = NULL;
    PLASMA_Complex64_t alpha;
    real_Double_t       t;
    int n       = iparam[TIMING_N];
    int nrhs    = n;
    int check   = iparam[TIMING_CHECK];
    int lda     = n;

    /* Allocate Data */
    A  = (PLASMA_Complex64_t *)malloc(lda*n   *sizeof(PLASMA_Complex64_t));
    B1 = (PLASMA_Complex64_t *)malloc(lda*nrhs*sizeof(PLASMA_Complex64_t));

    /* Check if unable to allocate memory */
    if ( (!A) || (!B1) ) {
        printf("Out of Memory \n ");
        exit(0);
    }

    /* Initialize Plasma */
    PLASMA_Init( iparam[TIMING_THRDNBR] );
    if ( iparam[TIMING_SCHEDULER] )
        PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_DYNAMIC_SCHEDULING );
    else
        PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_STATIC_SCHEDULING );

    /*if ( !iparam[TIMING_AUTOTUNING] ) {*/
        PLASMA_Disable(PLASMA_AUTOTUNING);
        PLASMA_Set(PLASMA_TILE_SIZE, iparam[TIMING_NB] );
    /* } */

     /* Initialiaze Data */
    lapack_zlarnv(1, ISEED, n   *lda,  A );
    lapack_zlarnv(1, ISEED, nrhs*lda,  B1);
    lapack_zlarnv(1, ISEED, 1,  &alpha);
    int i;
    for(i=0;i<max(n, nrhs);i++)
      A[lda*i+i] = A[lda*i+i] + 2.0;

    if (check)
    {
        B2 = (PLASMA_Complex64_t *)malloc(lda*nrhs*sizeof(PLASMA_Complex64_t));
        memcpy(B2, B1, lda*nrhs*sizeof(PLASMA_Complex64_t));
    }

    t = -cWtime();
    PLASMA_ztrsm( PlasmaLeft, PlasmaUpper, PlasmaNoTrans, PlasmaUnit,
          n, nrhs, alpha, A, lda, B1, lda );
    t += cWtime();
    *t_ = t;

    /* Check the solution */
    if (check)
      {
    dparam[TIMING_RES] = z_check_trsm( PlasmaLeft, PlasmaUpper, PlasmaNoTrans, PlasmaUnit, n, nrhs,
                      alpha, A, lda, B1, B2, lda,
                      &(dparam[TIMING_ANORM]), &(dparam[TIMING_BNORM]),
                      &(dparam[TIMING_XNORM]));
    free(B2);
      }

    free( A );
    free( B1 );

    PLASMA_Finalize();

    return 0;
}
