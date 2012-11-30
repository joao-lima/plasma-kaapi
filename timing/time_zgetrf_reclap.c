/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_zgetrf_reclap"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GETRF(n, nrhs)
#define _FADDS FADDS_GETRF(n, nrhs)

#include "../control/common.h"
#include "./timing.c"

void CORE_zgetrf_reclap_init(void);
extern plasma_context_t*  plasma_context_self(void);

/*
 * WARNING: the check is only working with LAPACK Netlib
 * which choose the same pivot than this code.
 * MKL has a different code and can pick a different pivot 
 * if two elments have the same absolute value but not the 
 * same sign for example.
 */

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    PLASMA_Complex64_t *A, *A2 = NULL;
    real_Double_t       t;
    int                *ipiv, *ipiv2 = NULL;
    int i;
    int m     = iparam[TIMING_N];
    int n     = iparam[TIMING_NRHS];
    int check = iparam[TIMING_CHECK];
    int lda   = m;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;

    /* Initialize Plasma */ 
    PLASMA_Init( iparam[TIMING_THRDNBR] );
    PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_DYNAMIC_SCHEDULING );

    PLASMA_Disable(PLASMA_AUTOTUNING);
    PLASMA_Set(PLASMA_TILE_SIZE,        iparam[TIMING_NB] );
    PLASMA_Set(PLASMA_INNER_BLOCK_SIZE, iparam[TIMING_IB] );

    /* Allocate Data */
    A  = (PLASMA_Complex64_t *)malloc(lda*n*sizeof(PLASMA_Complex64_t));

    /* Check if unable to allocate memory */
    if ( (! A) ) {
        printf("Out of Memory \n ");
        return -1;
    }

    /* Initialiaze Data */
    LAPACKE_zlarnv_work(1, ISEED, lda*n, A);

    /* Allocate Workspace */
    ipiv  = (int *)malloc( n*sizeof(int) );

    /* Save A in lapack layout for check */
    if ( check ) {
        A2 = (PLASMA_Complex64_t *)malloc(lda*n*sizeof(PLASMA_Complex64_t));
        ipiv2 = (int *)malloc( n*sizeof(int) );
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR,' ', m, n, A, lda, A2, lda);
    
        LAPACKE_zgetrf_work(LAPACK_COL_MAJOR, m, n, A2, lda, ipiv2 );
    }

    plasma = plasma_context_self();
    PLASMA_Sequence_Create(&sequence);
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    QUARK_Task_Flag_Set(&task_flags, TASK_THREAD_COUNT, iparam[TIMING_THRDNBR] );

    plasma_dynamic_spawn();
    CORE_zgetrf_reclap_init();

    t = -cWtime();
    QUARK_CORE_zgetrf_reclap(plasma->quark, &task_flags,
                             m, n, n,
                             A, lda, ipiv,
                             sequence, &request,
                             0, 0,
                             iparam[TIMING_THRDNBR]);
    PLASMA_Sequence_Wait(sequence);
    t += cWtime();
    *t_ = t;
    
    PLASMA_Sequence_Destroy(sequence);

    /* Check the solution */
    if ( check )
    {
        double *work = (double *)malloc(max(m,n)*sizeof(double));

        /* Check ipiv */
        for(i=0; i<n; i++)
        {
            if( ipiv[i] != ipiv2[i] ) {
                fprintf(stderr, "\nPLASMA (ipiv[%d] = %d, A[%d] = %e) / LAPACK (ipiv[%d] = %d, A[%d] = [%e])\n",
                        i, ipiv[i],  i, creal(A[  i * lda + i ]), 
                        i, ipiv2[i], i, creal(A2[ i * lda + i ])); 
                break;
            }
        }

        dparam[TIMING_ANORM] = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaMaxNorm), 
                                                   m, n, A, lda, work);
        dparam[TIMING_XNORM] = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaMaxNorm), 
                                                   m, n, A2, lda, work);
        dparam[TIMING_BNORM] = 0.0;

        CORE_zaxpy( m, n, -1.0, A, lda, A2, lda);

        dparam[TIMING_RES] = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaMaxNorm), 
                                                 m, n, A2, lda, work);

        free( A2 );
        free( ipiv2 );
        free( work );
    }
    
    free( A  );
    free( ipiv );
    PLASMA_Finalize();

    return 0;
}