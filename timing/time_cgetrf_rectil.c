/**
 *
 * @generated c Thu Sep 15 12:09:48 2011
 *
 **/
#define _TYPE  PLASMA_Complex32_t
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_cgetrf_rectil"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GETRF(n, nrhs)
#define _FADDS FADDS_GETRF(n, nrhs)

#include "../control/common.h"
#include "./timing.c"

void CORE_cgetrf_rectil_init(void);
extern plasma_context_t*  plasma_context_self(void);

/*
 * WARNING: the check is only working with LAPACK Netlib
 * which choose the same pivot than this code.
 * MKL has a different code and can pick a different pivot 
 * if two elments have the same absolute value but not the 
 * same sign for example.
 */

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    PLASMA_Complex32_t *A, *AT, *A2 = NULL;
    PLASMA_desc        *descA;
    real_Double_t       t;
    int                *ipiv, *ipiv2 = NULL;
    int i;
    int nb    = iparam[TIMING_NB];
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
    A  = (PLASMA_Complex32_t *)malloc(lda*n*sizeof(PLASMA_Complex32_t));
    AT = (PLASMA_Complex32_t *)malloc(lda*n*sizeof(PLASMA_Complex32_t));

    /* Check if unable to allocate memory */
    if ( ( !AT ) || (! A) ) {
        printf("Out of Memory \n ");
        return -1;
    }

    /* Initialiaze Data */
    LAPACKE_clarnv_work(1, ISEED, lda*n, A);
/*     for(i=0; i<n; i++) { */
/*       A[i*lda+i] += (float)m; */
/*     } */

    PLASMA_Desc_Create(&descA, AT, PlasmaComplexFloat, nb, nb, nb*nb, lda, n, 0, 0, m, n);
    PLASMA_cLapack_to_Tile((void*)A, lda, descA);

    /* Allocate Workspace */
    ipiv  = (int *)malloc( n*sizeof(int) );

    /* Save AT in lapack layout for check */
    if ( check ) {
        A2 = (PLASMA_Complex32_t *)malloc(lda*n*sizeof(PLASMA_Complex32_t));
        ipiv2 = (int *)malloc( n*sizeof(int) );
        LAPACKE_clacpy_work(LAPACK_COL_MAJOR,' ', m, n, A, lda, A2, lda);
    
        LAPACKE_cgetrf_work(LAPACK_COL_MAJOR, m, n, A2, lda, ipiv2 );
    }

    plasma = plasma_context_self();
    PLASMA_Sequence_Create(&sequence);
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    QUARK_Task_Flag_Set(&task_flags, TASK_THREAD_COUNT, iparam[TIMING_THRDNBR] );

    plasma_dynamic_spawn();
    CORE_cgetrf_rectil_init();

    t = -cWtime();
    QUARK_CORE_cgetrf_rectil(plasma->quark, &task_flags,
                             *descA, AT, descA->mb*descA->nb, ipiv,
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
        float *work = (float *)malloc(max(m,n)*sizeof(float));

        PLASMA_cTile_to_Lapack(descA, (void*)A, lda);

        /* Check ipiv */
        for(i=0; i<n; i++)
        {
            if( ipiv[i] != ipiv2[i] ) {
                fprintf(stderr, "\nPLASMA (ipiv[%d] = %d, A[%d] = %e) / LAPACK (ipiv[%d] = %d, A[%d] = [%e])\n",
                        i, ipiv[i],  i, crealf(A[  i * lda + i ]), 
                        i, ipiv2[i], i, crealf(A2[ i * lda + i ])); 
                break;
            }
        }

        dparam[TIMING_ANORM] = LAPACKE_clange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaMaxNorm), 
                                                   m, n, A, lda, work);
        dparam[TIMING_XNORM] = LAPACKE_clange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaMaxNorm), 
                                                   m, n, A2, lda, work);
        dparam[TIMING_BNORM] = 0.0;

        CORE_caxpy( m, n, -1.0, A, lda, A2, lda);

        dparam[TIMING_RES] = LAPACKE_clange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaMaxNorm), 
                                                 m, n, A2, lda, work);

        free( A2 );
        free( ipiv2 );
        free( work );
    }
    
    /* Deallocate Workspace */
    PLASMA_Desc_Destroy(&descA);

    free( A  );
    free( AT );
    free( ipiv );
    PLASMA_Finalize();

    return 0;
}
