/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_zgecfi"
/* See Lawn 41 page 120 */
#define _FMULS (0.0)
#define _FADDS (n * n * sizeof(_TYPE))

#include "./timing.c"

int z_check_conversion(int m, int n, int mba, int nba, int mbb, int nbb,
                      PLASMA_Complex64_t *A, PLASMA_Complex64_t *B, 
                      int (*mapA)(int, int, int, int, int, int), int (*mapB)(int, int, int, int, int, int)) {
    int i, j;

    for( j=0; j<n; j++) {
        for (i=0; i<m; i++) {
            if (A[ mapA(m, n, mba, nba, i, j) ] != B[ mapB(m, n, mbb, nbb, i, j) ] ) {
                return -1; 
            }
        }
    }
    return 0;
}

static int
RunTest(int *iparam, _PREC *dparam, real_Double_t *t_)
{
    PLASMA_Complex64_t *A, *Acpy = NULL;
    real_Double_t       t;
    int n       = iparam[TIMING_N]; 
    int nb      = iparam[TIMING_NB];
    int check   = iparam[TIMING_CHECK];

    n = ((n % nb) == 0) ? (n / nb) * nb : ((n / nb) + 1) * nb ;
    dparam[TIMING_ANORM] = (_PREC)n;
    dparam[TIMING_BNORM] = (_PREC)_FADDS;

    /* Allocate Data */
    A = (PLASMA_Complex64_t *)malloc(n*n*sizeof(PLASMA_Complex64_t));

    /* Check if unable to allocate memory */
    if ( (!A) ) {
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
    LAPACKE_zlarnv_work(1, ISEED, n*n, A);

    /* Save A and b  */
    if (check) {
        Acpy = (PLASMA_Complex64_t *)malloc(n*n*sizeof(PLASMA_Complex64_t));
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, lapack_const(PlasmaUpperLower), n, n, A, n, Acpy, n);
    }

    t = -cWtime();
    PLASMA_zgecfi( n, n, A, PlasmaCM, n, 1, PlasmaCCRB, nb, nb);
    t += cWtime();
    *t_ = t;
    
    /* Check the solution */
    if (check)
      {
        dparam[TIMING_RES] = (_PREC)z_check_conversion(n, n, n, 1, nb, nb, Acpy, A, map_CM, map_CCRB);
        free(Acpy);
      }

    free( A );

    PLASMA_Finalize();

    return 0;
}
