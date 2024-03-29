/**
 *
 * @file core_dgetrf_rectil.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Piotr Luszczek
 * @date 2009-11-15
 *
 * @generated d Thu Sep 15 12:09:01 2011
 * 
 **/

#include <math.h>
#include <cblas.h>
#include <lapacke.h>
#include "common.h"

#define A(m, n) BLKADDR(A, double, m, n)

static void CORE_dbarrier_thread(const int thidx, const int thcnt);
static void CORE_damax1_thread(const double localamx, 
                               const int thidx, const int thcnt, 
                               int *thwinner, double *diagvalue,
                               double *globalamx,
                               const int pividx, int *ipiv);

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 *  CORE_dgetf2 computes an LU factorization of a general M-by-N matrix A
 *  using partial pivoting with row interchanges.
 *
 *  WARNING: You cannot call this kernel on different matrices at the same time
 *
 *  The factorization has the form
 *
 *    A = P * L * U
 *
 *  where P is a permutation matrix, L is lower triangular with unit
 *  diagonal elements (lower trapezoidal if m > n), and U is upper
 *  triangular (upper trapezoidal if m < n).
 *
 *  This is the right-looking LAPACK Level 2 BLAS version of the algorithm.
 *
 *******************************************************************************
 *
 *  @param[in] M
 *          The number of rows of the matrix A.  M >= 0.
 *
 *  @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 *  @param[in,out] A
 *          On entry, the m by n matrix to be factored.
 *          On exit, the factors L and U from the factorization
 *          A = P*L*U; the unit diagonal elements of L are not stored.
 *
 *  @param[in] LDA
 *          The leading dimension of the array A.  LDA >= max(1,M).
 *
 *  @param[out] IPIV
 *          The pivot indices; for 1 <= i <= min(M,N), row i of the
 *          matrix was interchanged with row IPIV(i).
 *
 *  @param[out] INFO
 *          = k if U(k,k) is exactly zero. The factorization
 *               has been completed, but the factor U is exactly
 *               singular, and division by zero will occur if it is used
 *               to solve a system of equations.
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval -k, the k-th argument had an illegal value
 *
 */

#define AMAX1BUF_SIZE (48 << 1)

/* 48 threads should be enough for everybody */
static volatile double CORE_damax1buf[AMAX1BUF_SIZE]; 
static double sfmin;

void 
CORE_dgetrf_rectil_init(void) { 
    int i;
    for (i = 0; i < AMAX1BUF_SIZE; ++i) CORE_damax1buf[i] = -1.0;
    sfmin =  LAPACKE_dlamch_work('S');
}

static void
CORE_damax1_thread(double localamx, int thidx, int thcnt, int *thwinner, 
                   double *diagvalue, double *globalamx, int pividx, int *ipiv) {
    if (thidx == 0) {
        int i, j = 0;
        double curval = localamx, tmp;
        double curamx = fabs(localamx);
        
        /* make sure everybody filled in their value */
        for (i = 1; i < thcnt; ++i) {
            while (CORE_damax1buf[i << 1] == -1.0) { /* wait for thread i to store its value */
            }
        }
        
        /* better not fuse the loop above and below to make sure data is sync'd */
       
        for (i = 1; i < thcnt; ++i) {
            tmp = CORE_damax1buf[ (i << 1) + 1];
            if (fabs(tmp) > curamx) {
                curamx = fabs(tmp);
                curval = tmp;
                j = i;
            }
        }
        
        if (0 == j)
            ipiv[0] = pividx;
        
        /* make sure everybody knows the amax value */
        for (i = 1; i < thcnt; ++i)
            CORE_damax1buf[ (i << 1) + 1] = curval;
        
        CORE_damax1buf[0] = -j - 2.0; /* set the index of the winning thread */
        CORE_damax1buf[1] = *diagvalue; /* set the index of the winning thread */
        
        *thwinner = j;
        *globalamx = curval;
        
        for (i = 1; i < thcnt; ++i)
            CORE_damax1buf[i << 1] = -3.0;
        
        /* make sure everybody read the max value */
        for (i = 1; i < thcnt; ++i) {
            while (CORE_damax1buf[i << 1] != -1.0) {
            }
        }
        
        CORE_damax1buf[0] = -1.0;
    } else {
        CORE_damax1buf[(thidx << 1) + 1] = localamx;
        CORE_damax1buf[thidx << 1] = -2.0;  /* announce to thread 0 that local amax was stored */
        while (CORE_damax1buf[0] == -1.0) { /* wait for thread 0 to finish calculating the global amax */
        }
        while (CORE_damax1buf[thidx << 1] != -3.0) { /* wait for thread 0 to store amax */
        }
        *thwinner  = -CORE_damax1buf[0] - 2.0;
        *diagvalue = CORE_damax1buf[1];
        *globalamx = CORE_damax1buf[(thidx << 1) + 1]; /* read the amax from the location adjacent to the one in the above loop */
        CORE_damax1buf[thidx << 1] = -1.0;  /* signal thread 0 that this thread is done reading */

        if (thidx == *thwinner)
            ipiv[0] = pividx;

        while (CORE_damax1buf[0] != -1.0) { /* wait for thread 0 to finish */
        }
    }
}

static void
CORE_dbarrier_thread(int thidx, int thcnt) {
    int idum1, idum2;
    double ddum1 = 0.;
    double ddum2 = 0.;
    /* it's probably faster to implement a dedicated barrier */
    CORE_damax1_thread( 1.0, thidx, thcnt, &idum1, &ddum1, &ddum2, 0, &idum2 );
}

static void
CORE_dgetrf_rectil_update(const PLASMA_desc A, int *IPIV, 
                          const int column, const int n1,     const int n2,
                          const int thidx,  const int thcnt, 
                          const int ft,     const int lt)
{
    int ld, lm, tmpM;
    int ip, j, it, i, ldft;
    double zone  = 1.0;
    double mzone = -1.0;
    double *Atop, *Atop2, *U, *L;
    int offset = A.i;

    ldft = BLKLDD(A, 0);
    Atop  = A(0, 0) + column * ldft;
    Atop2 = Atop    + n1     * ldft;

    if (thidx == 0)
    {
        /* Swap to the right */
        int *lipiv = IPIV+column;
        int idxMax = column+n1;
        for (j = column; j < idxMax; ++j, ++lipiv) {
            ip = (*lipiv) - offset - 1;
            if ( ip != j )
            {
                it = ip / A.mb;
                i  = ip % A.mb;
                ld = BLKLDD(A, it);
                cblas_dswap(n2, Atop2                     + j, ldft,
                            A(it, 0) + (column+n1)*ld + i, ld   );
            }
        }
        
        /* Trsm on the upper part */
        U = Atop2 + column;
        cblas_dtrsm( CblasColMajor, CblasLeft, CblasLower, 
                     CblasNoTrans, CblasUnit,
                     n1, n2, (zone), 
                     Atop  + column, ldft, 
                     U,              ldft );
        
        /* Signal to other threads that they can start update */
        CORE_dbarrier_thread( thidx, thcnt );

        /* First tile */
        L = Atop + column + n1;
        tmpM = ldft - column - n1;
        
        /* Apply the GEMM */
        cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, 
                     tmpM, n2, n1,
                     (mzone), L,      ldft, 
                     U,      ldft, 
                     (zone),  U + n1, ldft );
        
    }
    else 
    {
        ld = BLKLDD( A, ft );
        L  = A( ft, 0 ) + column * ld;
        lm = ft == A.mt-1 ? A.m - ft * A.mb : A.mb;
        
        U = Atop2 + column;
        
        /* Wait for pivoting and triangular solve to be finished
         * before to really start the update */
        CORE_dbarrier_thread( thidx, thcnt );

        /* First tile */
        /* Apply the GEMM */
        cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, 
                     lm, n2, n1,
                     (mzone), L,         ld, 
                                         U,         ldft, 
                     (zone),  L + n1*ld, ld );
    }
        
    /* Update the other blocks */
    for( it = ft+1; it < lt; it++)
    {
        ld = BLKLDD( A, it );
        L  = A( it, 0 ) + column * ld;
        lm = it == A.mt-1 ? A.m - it * A.mb : A.mb;
        
        /* Apply the GEMM */
        cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, 
                     lm, n2, n1,
                     (mzone), L,          ld, 
                                         U,          ldft, 
                     (zone),  L + n1*ld,  ld );
    }
}

static void 
CORE_dgetrf_rectil_rec(const PLASMA_desc A, int *IPIV, int *info, 
                       double *pivot,
                       const int thidx,  const int thcnt, 
                       const int column, const int width, 
                       const int ft,     const int lt)
{
    int ld, jp, n1, n2, lm, tmpM, piv_sf;
    int ip, j, it, i, ldft;
    int max_i, max_it, thwin;
    double zone  = 1.0;
    double mzone = -1.0;
    double tmp1;
    double tmp2;
    double pivval;
    double *Atop, *Atop2, *U, *L;
    double             abstmp1;
    int offset = A.i;

    ldft = BLKLDD(A, 0);
    Atop = A(0, 0) + column * ldft;

#if 0
    CORE_dbarrier_thread( thidx, thcnt );
    if (thidx == 0)
    {
        fprintf(stderr, "\n ------  column = %d / width = %d -------\n", column, width);
        int i, j;
        for(j=0;j<4;j++){
            for(i=0;i<4;i++){
                fprintf(stderr, "%e ", ((double*)A(0, 0))[j*4+i]);
            }
            for(i=0;i<4;i++){
                fprintf(stderr, "%e ", ((double*)A(1, 0))[j*4+i]);
            }
            fprintf(stderr, "\n");
        }
    }
    CORE_dbarrier_thread( thidx, thcnt );
#endif

    if ( width > 1 ) {
        /* Assumption: N = min( M, N ); */
        n1 = width / 2;
        n2 = width - n1;
        
        Atop2 = Atop + n1 * ldft;
        
        CORE_dgetrf_rectil_rec( A, IPIV, info, pivot,
                                thidx, thcnt, column, n1, ft, lt );

        if ( *info != 0 )
            return;

        if (thidx == 0)
        {
            /* Swap to the right */
            int *lipiv = IPIV+column;
            int idxMax = column+n1;
            for (j = column; j < idxMax; ++j, ++lipiv) {
                ip = (*lipiv) - offset - 1;
                if ( ip != j )
                {
                    it = ip / A.mb;
                    i  = ip % A.mb;
                    ld = BLKLDD(A, it);
                    cblas_dswap(n2, Atop2                     + j, ldft,
                                    A(it, 0) + (column+n1)*ld + i, ld   );
                }
            }   
            /* Trsm on the upper part */
            U = Atop2 + column;
            cblas_dtrsm( CblasColMajor, CblasLeft, CblasLower, 
                         CblasNoTrans, CblasUnit,
                         n1, n2, (zone), 
                         Atop  + column, ldft, 
                         U,              ldft );

            /* SIgnal to other threads that they can start update */
            CORE_dbarrier_thread( thidx, thcnt );
            pivval = *pivot;
            if ( pivval == 0.0 ) {
                *info = column+n1;
                return;
            } else {
                if ( fabs(pivval) >= sfmin ) {
                    piv_sf = 1;
                    pivval = 1.0 / pivval;
                } else {
                    piv_sf = 0;
                }
            }
                    
            /* First tile */
            { 
                L = Atop + column + n1;
                tmpM = ldft - column - n1;
                
                /* Scale last column of L */
                if ( piv_sf ) {
                    cblas_dscal( tmpM, (pivval), L+(n1-1)*ldft, 1 );
                } else {
                    int i;
                    Atop2 = L+(n1-1)*ldft;
                    for( i=0; i < tmpM; i++, Atop2++)
                        *Atop2 = *Atop2 / pivval;
                }

                /* Apply the GEMM */
                cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, 
                             tmpM, n2, n1,
                             (mzone), L,      ldft, 
                                                 U,      ldft, 
                             (zone),  U + n1, ldft );
                
                /* Search Max in first column of U+n1 */
                tmp2    = U[n1];
                max_it  = ft;
                max_i   = cblas_idamax( tmpM, U+n1, 1 ) + n1;
                tmp1    = U[max_i];
                abstmp1 = fabs(tmp1);
                max_i   += column;
            }
        }
        else 
        {
            pivval = *pivot;
            if ( pivval == 0.0 ) {
                *info = column+n1;
                return;
            } else {
                if ( fabs(pivval) >= sfmin ) {
                    piv_sf = 1;
                    pivval = 1.0 / pivval;
                } else {
                    piv_sf = 0;
                }
            }
                    
            ld = BLKLDD( A, ft );
            L  = A( ft, 0 ) + column * ld;
            lm = ft == A.mt-1 ? A.m - ft * A.mb : A.mb;
            
            U = Atop2 + column;

            /* First tile */
            /* Scale last column of L */
            if ( piv_sf ) {
                cblas_dscal( lm, (pivval), L+(n1-1)*ld, 1 );
            } else {
                int i;
                Atop2 = L+(n1-1)*ld;
                for( i=0; i < lm; i++, Atop2++)
                    *Atop2 = *Atop2 / pivval;
            }

            /* Wait for pivoting and triangular solve to be finished
             * before to really start the update */
            CORE_dbarrier_thread( thidx, thcnt );

            /* Apply the GEMM */
            cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, 
                         lm, n2, n1,
                         (mzone), L,         ld, 
                                             U,         ldft, 
                         (zone),  L + n1*ld, ld );
                
            /* Search Max in first column of L+n1*ld */
            max_it  = ft;
            max_i   = cblas_idamax( lm, L+n1*ld, 1 );
            tmp1    = L[n1*ld+max_i];
            abstmp1 = fabs(tmp1);
        }
        
        /* Update the other blocks */
        for( it = ft+1; it < lt; it++)
        {
            ld = BLKLDD( A, it );
            L  = A( it, 0 ) + column * ld;
            lm = it == A.mt-1 ? A.m - it * A.mb : A.mb;
            
            /* Scale last column of L */
            if ( piv_sf ) {
                cblas_dscal( lm, (pivval), L+(n1-1)*ld, 1 );
            } else {
                int i;
                Atop2 = L+(n1-1)*ld;
                for( i=0; i < lm; i++, Atop2++)
                    *Atop2 = *Atop2 / pivval;
            }
            
            /* Apply the GEMM */
            cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, 
                         lm, n2, n1,
                         (mzone), L,          ld, 
                                             U,          ldft, 
                         (zone),  L + n1*ld,  ld );
            
            /* Search the max on the first column of L+n1*ld */
            jp = cblas_idamax( lm, L+n1*ld, 1 );
            if ( fabs( L[n1*ld+jp] ) > abstmp1 ) {
                tmp1 = L[n1*ld+jp];
                abstmp1 = fabs(tmp1);
                max_i  = jp;
                max_it = it;
            }
        }

        jp = offset + max_it*A.mb + max_i;
        CORE_damax1_thread( tmp1, thidx, thcnt, &thwin, 
                            &tmp2, pivot, jp + 1, IPIV + column + n1 );
        
        if ( thidx == 0 ) {
            U[n1] = *pivot; /* all threads have the pivot element: no need for synchronization */
        } 
        if (thwin == thidx) { /* the thread that owns the best pivot */
            if ( jp-offset != column+n1 ) /* if there is a need to exchange the pivot */
            {
                ld = BLKLDD(A, max_it);
                Atop2 = A( max_it, 0 ) + (column + n1 )* ld + max_i;
                *Atop2 = tmp2;
            }
        }
        
        CORE_dgetrf_rectil_rec( A, IPIV, info, pivot,
                                thidx, thcnt, column+n1, n2, ft, lt );
        if ( *info != 0 )
            return;

        if ( thidx == 0 )
        {
            /* Swap to the left */
            int *lipiv = IPIV+column+n1;
            int idxMax = column+width;
            for (j = column+n1; j < idxMax; ++j, ++lipiv) {
                ip = (*lipiv) - offset - 1;
                if ( ip != j )
                {
                    it = ip / A.mb;
                    i  = ip % A.mb;
                    ld = BLKLDD(A, it);
                    cblas_dswap(n1, Atop + j,                 ldft,
                                    A(it, 0) + column*ld + i, ld  );
                }
            }   
        }

    } else if ( width == 1 ) {

        /* Search maximum for column 0 */
        if ( column == 0 )
        {
            if ( thidx == 0 )
              tmp2 = Atop[column]; 
            
            /* First tmp1 */
            ld = BLKLDD(A, ft);
            Atop2   = A( ft, 0 );
            lm      = ft == A.mt-1 ? A.m - ft * A.mb : A.mb;
            max_it  = ft;
            max_i   = cblas_idamax( lm, Atop2, 1 );
            tmp1    = Atop2[max_i];
            abstmp1 = fabs(tmp1);
            
            /* Update */
            for( it = ft+1; it < lt; it++)
            {
                Atop2= A( it, 0 );
                lm   = it == A.mt-1 ? A.m - it * A.mb : A.mb;
                jp   = cblas_idamax( lm, Atop2, 1 );
                if (  fabs(Atop2[jp]) > abstmp1 ) {
                    tmp1 = Atop2[jp];
                    abstmp1 = fabs(tmp1);
                    max_i  = jp;
                    max_it = it;
                }
            }
            
            jp = offset + max_it*A.mb + max_i;
            CORE_damax1_thread( tmp1, thidx, thcnt, &thwin, 
                                &tmp2, pivot, jp + 1, IPIV + column );
            
            if ( thidx == 0 ) {
                Atop[0] = *pivot; /* all threads have the pivot element: no need for synchronization */
            } 
            if (thwin == thidx) { /* the thread that owns the best pivot */
                if ( jp-offset != 0 ) /* if there is a need to exchange the pivot */
                {
                    Atop2 = A( max_it, 0 ) + max_i;
                    *Atop2 = tmp2;
                }
            }
        }

        CORE_dbarrier_thread( thidx, thcnt );

        /* If it is the last column, we just scale */
        if ( column == (min(A.m, A.n))-1 ) {
        
            pivval = *pivot;
            if ( pivval != 0.0 ) {
                if ( thidx == 0 ) {
                    if ( fabs(pivval) >= sfmin ) {
                        pivval = 1.0 / pivval;
                        
                        /*
                         * We guess than we never enter the function with m == A.mt-1 
                         * because it means that there is only one thread 
                         */
                        lm = ft == A.mt-1 ? A.m - ft * A.mb : A.mb;
                        cblas_dscal( lm - column - 1, (pivval), Atop+column+1, 1 );
                        
                        for( it = ft+1; it < lt; it++)
                        {
                            ld = BLKLDD(A, it);
                            Atop2 = A( it, 0 ) + column * ld;
                            lm = it == A.mt-1 ? A.m - it * A.mb : A.mb;
                            cblas_dscal( lm, (pivval), Atop2, 1 );
                        }
                    } else {
                        /*
                         * We guess than we never enter the function with m == A.mt-1 
                         * because it means that there is only one thread 
                         */
                        int i;
                        Atop2 = Atop + column + 1;
                        lm = ft == A.mt-1 ? A.m - ft * A.mb : A.mb;
                        
                        for( i=0; i < lm-column-1; i++, Atop2++)
                            *Atop2 = *Atop2 / pivval;
                        
                        for( it = ft+1; it < lt; it++)
                            {
                                ld = BLKLDD(A, it);
                                Atop2 = A( it, 0 ) + column * ld;
                                lm = it == A.mt-1 ? A.m - it * A.mb : A.mb;
                                
                                for( i=0; i < lm; i++, Atop2++)
                                    *Atop2 = *Atop2 / pivval;
                            }
                    }
                }
                else
                {
                    if ( fabs(pivval) >= sfmin ) {
                        pivval = 1.0 / pivval;
                        
                        for( it = ft; it < lt; it++)
                        {
                            ld = BLKLDD(A, it);
                            Atop2 = A( it, 0 ) + column * ld;
                            lm = it == A.mt-1 ? A.m - it * A.mb : A.mb;
                            cblas_dscal( lm, (pivval), Atop2, 1 );
                        }
                    } else {
                        /*
                         * We guess than we never enter the function with m == A.mt-1 
                         * because it means that there is only one thread 
                         */
                        int i;
                        for( it = ft; it < lt; it++)
                        {
                            ld = BLKLDD(A, it);
                            Atop2 = A( it, 0 ) + column * ld;
                            lm = it == A.mt-1 ? A.m - it * A.mb : A.mb;
                            
                            for( i=0; i < lm; i++, Atop2++)
                                *Atop2 = *Atop2 / pivval;
                        }
                    }
                }
            }
            else {
                *info = column + 1;
                return;
            }
        }
    }
}

#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dgetrf_rectil = PCORE_dgetrf_rectil
#define CORE_dgetrf_rectil PCORE_dgetrf_rectil
#endif
int CORE_dgetrf_rectil(const PLASMA_desc A, int *IPIV, int *info)
{
    int ft, lt; 
    int thidx = info[1];
    int thcnt = min( info[2], A.mt );
    int minMN = min( A.m, A.n );
    double pivot;

    if ( A.nt > 1 ) {
        coreblas_error(1, "Illegal value of A.nt");
        return -1;
    }

    if ( thidx >= thcnt )
      return 0;

    int q = A.mt / thcnt;
    int r = A.mt % thcnt;

    if (thidx < r) {
        q++;
        ft = thidx * q;
        lt = ft + q;
    } else {
        ft = r * (q + 1) + (thidx - r) * q;
        lt = ft + q;
        lt = min( lt, A.mt );
    }
    
    info[0] = 0;
    CORE_dgetrf_rectil_rec( A, IPIV, info, &pivot,
                            thidx, thcnt, 0, minMN, ft, lt);
   
    if ( A.n > minMN ) {
        CORE_dgetrf_rectil_update( A, IPIV, 
                                   0, minMN, A.n-minMN,
                                   thidx, thcnt, 
                                   ft, lt);
    }

    return info[0];
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dgetrf_rectil(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc A, double *Amn, int size,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo,
                              int nbthread)
{
    DAG_CORE_GETRF;
    QUARK_Insert_Task(quark, CORE_dgetrf_rectil_quark, task_flags,
        sizeof(PLASMA_desc),                &A,             VALUE,
        sizeof(double)*size,     Amn,               INOUT,
        sizeof(int)*A.n,                     IPIV,              OUTPUT,
        sizeof(PLASMA_sequence*),           &sequence,      VALUE,
        sizeof(PLASMA_request*),            &request,       VALUE,
        sizeof(PLASMA_bool),                &check_info,    VALUE,
        sizeof(int),                        &iinfo,         VALUE,
        sizeof(int),                        &nbthread,      VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dgetrf_rectil_quark = PCORE_dgetrf_rectil_quark
#define CORE_dgetrf_rectil_quark PCORE_dgetrf_rectil_quark
#endif
void CORE_dgetrf_rectil_quark(Quark* quark)
{
    PLASMA_desc A;
    double *Amn;
    int *IPIV;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    PLASMA_bool check_info;
    int iinfo;

    int info[3];
    int maxthreads;

    quark_unpack_args_8(quark, A, Amn, IPIV, sequence, request, 
                        check_info, iinfo, maxthreads );

    info[1] = QUARK_Get_RankInTask(quark);
    info[2] = maxthreads;

    CORE_dgetrf_rectil( A, IPIV, info );
    if (info[1] == 0 && info[0] != PLASMA_SUCCESS && check_info)
        plasma_sequence_flush(quark, sequence, request, iinfo + info[0] );
}
