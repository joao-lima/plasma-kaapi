/**
 *
 * @generated c Thu Sep 15 12:09:48 2011
 *
 **/
#define _TYPE  PLASMA_Complex32_t
#define _PREC  float
#define _LAMCH LAPACKE_slamch

#define _NAME  "PLASMA_clapack_to_tile"
/* See Lawn 41 page 120 */
#define _FMULS (0.0)
#define _FADDS (n * n * sizeof(_TYPE))

#include "./timing.c"

int c_check_conversion(int m, int n, int mba, int nba, int mbb, int nbb,
                      PLASMA_Complex32_t *A, PLASMA_Complex32_t *B, 
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
    PLASMA_Complex32_t *A = NULL, *AT;
    PLASMA_desc        *descA;
    real_Double_t       t;
    int n       = iparam[TIMING_N];
    int nb      = iparam[TIMING_NB];
    int check   = iparam[TIMING_CHECK];

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

    n = ((n % nb) == 0) ? (n / nb) * nb : ((n / nb) + 1) * nb ;
    dparam[TIMING_ANORM] = (_PREC)n;

    /* Allocate Data */
    AT = (PLASMA_Complex32_t *)malloc(n*n*sizeof(PLASMA_Complex32_t));

    /* Check if unable to allocate memory */
    if ( (!AT) ) {
        printf("Out of Memory \n ");
        exit(0);
    }
    
     /* Initialiaze Data */
    PLASMA_Desc_Create(&descA, AT, PlasmaComplexFloat, nb, nb, nb*nb, n, n, 0, 0, n, n);
    LAPACKE_clarnv_work(1, ISEED, n*n, AT);

    /* Save A and b  */
    if (check) {
        A = (PLASMA_Complex32_t *)malloc(n*n*sizeof(PLASMA_Complex32_t));
        LAPACKE_clacpy_work(LAPACK_COL_MAJOR, lapack_const(PlasmaUpperLower), n, n, AT, n, A, n);
    }

    t = -cWtime();
    PLASMA_Lapack_to_Tile( (void *)A, n, descA);
    t += cWtime();
    *t_ = t;
    
    /* Check the solution */
    if (check)
      {
        dparam[TIMING_RES] = (_PREC)c_check_conversion(n, n, n, 1, nb, nb, A, AT, map_CM, map_CCRB);
        free(A);
      }

    PLASMA_Desc_Destroy(&descA);
    free( AT );

    PLASMA_Finalize();

    return 0;
}
