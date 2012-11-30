/**
 *
 * @file core_sgetrf_reclap.c
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
 * @generated s Thu Sep 15 12:09:01 2011
 * 
 **/

#include <math.h>
#include <cblas.h>
#include <lapacke.h>
#include "common.h"

static void CORE_sbarrier_thread(const int thidx, const int thcnt);
static void CORE_samax1_thread(const float localamx, 
                               const int thidx, const int thcnt, 
                               int *thwinner, float *globalamx,
                               const int pividx, int *ipiv);

/* Laswp with inc = 1 */
static void CORE_slaswap1(const int ncol, float *a, const int lda, 
                          const int idxStart, const int idxMax, const int *piv);

/***************************************************************************//**
 *
 * @ingroup CORE_float
 *
 *  CORE_sgetrf_reclap computes an LU factorization of a general M-by-N matrix A
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
static void
psplit(int n, int pidx, int pcnt, int *poff_p, int *psiz_p) {
    int q = n / pcnt, r = n % pcnt;
    
    if (pidx < r) {
        q++;
        *psiz_p = q;
        *poff_p = pidx * q;
    } else {
        *psiz_p = q;
        *poff_p = r * (q + 1) + (pidx - r) * q;
    }
}

#define AMAX1BUF_SIZE (48 << 1)
/* 48 threads should be enough for everybody */
static volatile float CORE_samax1buf[AMAX1BUF_SIZE]; 
static float sfmin;

void
CORE_sgetrf_reclap_init(void) { 
    int i;

    for (i = 0; i < AMAX1BUF_SIZE; ++i) CORE_samax1buf[i] = -1.0;
    sfmin =  LAPACKE_slamch_work('S');
}

static void
CORE_samax1_thread(float localamx, int thidx, int thcnt, int *thwinner, 
                   float *globalamx, int pividx, int *ipiv) {
    if (thidx == 0) {
        int i, j = 0;
        float curval = localamx, tmp;
        float curamx = fabsf(localamx);
        
        /* make sure everybody filled in their value */
        for (i = 1; i < thcnt; ++i) {
            while (CORE_samax1buf[i << 1] == -1.0) { /* wait for thread i to store its value */
            }
        }
        
        /* better not fuse the loop above and below to make sure data is sync'd */
        
        for (i = 1; i < thcnt; ++i) {
            tmp = CORE_samax1buf[ (i << 1) + 1];
            if (fabsf(tmp) > curamx) {
                curamx = fabsf(tmp);
                curval = tmp;
                j = i;
            }
        }
        
        if (0 == j)
            ipiv[0] = pividx;
        
        /* make sure everybody knows the amax value */
        for (i = 1; i < thcnt; ++i)
            CORE_samax1buf[ (i << 1) + 1] = curval;

        CORE_samax1buf[0] = -j - 2.0; /* set the index of the winning thread */
        
        *thwinner = j;
        *globalamx = curval;
        
        for (i = 1; i < thcnt; ++i)
            CORE_samax1buf[i << 1] = -3.0;
        
        /* make sure everybody read the max value */
        for (i = 1; i < thcnt; ++i) {
            while (CORE_samax1buf[i << 1] != -1.0) {
            }
        }
        
        CORE_samax1buf[0] = -1.0;
    } else {
        CORE_samax1buf[(thidx << 1) + 1] = localamx;
        CORE_samax1buf[thidx << 1] = -2.0;  /* announce to thread 0 that local amax was stored */
        while (CORE_samax1buf[0] == -1.0) { /* wait for thread 0 to finish calculating the global amax */
        }
        while (CORE_samax1buf[thidx << 1] != -3.0) { /* wait for thread 0 to store amax */
        }
        *globalamx = CORE_samax1buf[(thidx << 1) + 1]; /* read the amax from the location adjacent to the one in the above loop */
        *thwinner = -CORE_samax1buf[0] - 2.0;
        CORE_samax1buf[thidx << 1] = -1.0;  /* signal thread 0 that this thread is done reading */
        
        if (thidx == *thwinner)
            ipiv[0] = pividx;
        
        while (CORE_samax1buf[0] != -1.0) { /* wait for thread 0 to finish */
        }
    }
}


static inline void CORE_sgetrf_reclap_update(const int M, const int column, const int n1, const int n2,
                                             float *A, const int LDA, int *IPIV, 
                                             const int thidx, const int thcnt)
{
    static float posone =  1.0;
    static float negone = -1.0;
    float *Atop  = A    + column*LDA;
    float *Atop2 = Atop + n1    *LDA;
    int coff, ccnt, lm, loff;

    CORE_sbarrier_thread( thidx, thcnt );
    
    psplit( n2, thidx, thcnt, &coff, &ccnt );

    if (ccnt > 0) {
        CORE_slaswap1( ccnt, Atop2 + coff*LDA, LDA, column, n1 + column, IPIV ); /* swap to the right */
        
        cblas_strsm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
                     n1, ccnt, (posone), Atop + column, LDA, Atop2 + coff*LDA + column, LDA );
    }
    
    /* __sync_synchronize(); */ /* hopefully we will not need memory fences */
    
    /* need to wait for pivoting and triangular solve to finish */
    CORE_sbarrier_thread( thidx, thcnt );
    
    psplit( M, thidx, thcnt, &loff, &lm );
    if (thidx == 0) {
        loff = column + n1;
        lm  -= column + n1;
    };
    
    cblas_sgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, lm, n2, n1,
                 (negone), Atop+loff, LDA, Atop2 + column, LDA, (posone), Atop2+loff, LDA );
}

static void 
CORE_sgetrf_reclap_rec(const int M, const int N, 
                       float *A, const int LDA, 
                       int *IPIV, int *info, 
                       const int thidx, const int thcnt, const int column)
{
    int jp, n1, n2, lm, loff;
    float tmp1, tmp2, tmp3;
    float *Atop = A + column*LDA;
    
    /* Assumption: N = min( M, N ); */
    if (N > 1) {
        int coff, ccnt;
        
        n1 = N / 2;
        n2 = N - n1;
        
        CORE_sgetrf_reclap_rec( M, n1, A, LDA, IPIV, info, 
                                thidx, thcnt, column );
        if ( *info != 0 )
            return;
        
        CORE_sgetrf_reclap_update(M, column, n1, n2,
                                  A, LDA, IPIV, 
                                  thidx, thcnt);
        
        CORE_sgetrf_reclap_rec( M, n2, A, LDA, IPIV, info, 
                                thidx, thcnt, column + n1 );
        if ( *info != 0 )
            return;
        
        psplit( n1, thidx, thcnt, &coff, &ccnt );
        
        if (ccnt > 0) {
            CORE_slaswap1( ccnt, Atop+coff*LDA, LDA, n1 + column, N + column, IPIV ); /* swap to the left */
        }
        
    } else {
        int thrd;
        
        CORE_sbarrier_thread( thidx, thcnt );
        
        psplit( M, thidx, thcnt, &loff, &lm );
        
        if (thidx == 0) {
            loff = column;
            lm -= column;
        }
        
        tmp2 = Atop[column]; /* all threads read the pivot element in case they need it */
        
        jp = cblas_isamax( lm, Atop + loff, 1 );
        tmp1 = Atop[loff + jp];
        
        CORE_samax1_thread( tmp1, thidx, thcnt, &thrd, 
                            &tmp3, loff + jp + 1, IPIV + column );
        
        Atop[column] = tmp3; /* all threads set the pivot element: no need for synchronization */
        
        if ( tmp3 != 0.0 ) {
            if ( fabsf(tmp3) >= sfmin ) {
                float tmp = (float)1.0 / tmp3;
                n1 = (thidx == 0) ? 1 : 0;
                cblas_sscal( lm - n1, (tmp), Atop + loff + n1, 1 );
            } else {
                int i;
                float *Atop2;
                n1 = (thidx == 0) ? 1 : 0;
                Atop2 = Atop + loff + n1;

                for( i=0; i < lm-n1; i++, Atop2++)
                    *Atop2 = *Atop2 / tmp3;
            }

            if (thrd == thidx) { /* the thread that owns the best pivot */
              if (loff + jp != column) /* if there is a need to exchange the pivot */
                Atop[loff + jp] = tmp2 / tmp3;
            }
        
        } else {
            *info = column + 1;
            return;
        }

        CORE_sbarrier_thread( thidx, thcnt );
    }
}

#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sgetrf_reclap = PCORE_sgetrf_reclap
#define CORE_sgetrf_reclap PCORE_sgetrf_reclap
#endif
int CORE_sgetrf_reclap(int M, int N, 
                       float *A, int LDA, 
                       int *IPIV, int *info)
{
    int thidx = info[1];
    int thcnt = min( info[2], M / N );
    int minMN = min(M, N);

    if( M < 0 ) {
        coreblas_error(1, "illegal value of M");
        return -1;
    }
    if( N < 0 ) {
        coreblas_error(2, "illegal value of N");
        return -2;
    }
    if( LDA < max(1, M) ) {
        coreblas_error(5, "illegal value of LDA");
        return -5;
    }
    
    /*
     * Quick return
     */
    if ( (M == 0) || (N == 0) || (thidx >= thcnt) ){
      return PLASMA_SUCCESS;
    }

    *info = 0;
    CORE_sgetrf_reclap_rec( M, minMN, A, LDA, IPIV, info, 
                            thidx, thcnt, 0 );

    if ( N > minMN ) {
        CORE_sgetrf_reclap_update(M, 0, minMN, N-minMN,
                                  A, LDA, IPIV, 
                                  thidx, thcnt);
    }

    return info[0];
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sgetrf_reclap(Quark *quark, Quark_Task_Flags *task_flags,
                              int m, int n, int nb,
                              float *A, int lda,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo,
                              int nbthread)
{
    DAG_CORE_GETRF;
    QUARK_Insert_Task(quark, CORE_sgetrf_reclap_quark, task_flags,
        sizeof(int),                        &m,             VALUE,
        sizeof(int),                        &n,             VALUE,
        sizeof(float)*nb*nb,    A,                     INOUT,
        sizeof(int),                        &lda,           VALUE,
        sizeof(int)*nb,                      IPIV,                  OUTPUT,
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
#pragma weak CORE_sgetrf_reclap_quark = PCORE_sgetrf_reclap_quark
#define CORE_sgetrf_reclap_quark PCORE_sgetrf_reclap_quark
#endif
void CORE_sgetrf_reclap_quark(Quark* quark)
{
    int M;
    int N;
    float *A;
    int LDA;
    int *IPIV;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    PLASMA_bool check_info;
    int iinfo;

    int info[3];
    int maxthreads;

    quark_unpack_args_10(quark, M, N, A, LDA, IPIV, sequence, request, 
                         check_info, iinfo, maxthreads );

    info[1] = QUARK_Get_RankInTask(quark);
    info[2] = maxthreads;

    CORE_sgetrf_reclap( M, N, A, LDA, IPIV, info );
    if (info[1] == 0 && info[0] != PLASMA_SUCCESS && check_info)
        plasma_sequence_flush(quark, sequence, request, iinfo + info[0] );
}


/*******************************************************************
 *   Additional routines
 */

static void
CORE_sbarrier_thread(int thidx, int thcnt) 
{
    int idum1, idum2;
    float ddum2;
    /* it's probably faster to implement a dedicated barrier */
    CORE_samax1_thread( 1.0, thidx, thcnt, &idum1, &ddum2, 0, &idum2 );
}

static void
CORE_slaswap1(const int ncol, float *a, const int lda, 
              const int idxStart, const int idxMax, const int *piv) 
{
    int i, j;
    float tmp;
    
    for (j = 0; j < ncol; ++j) {
        for (i = idxStart; i < idxMax; ++i) {
            tmp = a[j*lda + piv[i] - 1];
            a[j*lda + piv[i] - 1] = a[i + j*lda];
            a[i + j*lda] = tmp;
        }
    }
}

