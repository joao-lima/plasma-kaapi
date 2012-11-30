/**
 *
 * @file pclansy.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Emmanuel Agullo
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated c Thu Sep 15 12:09:21 2011
 *
 **/
#include <stdlib.h>
#include <math.h>
#include "common.h"

#define A(m, n, i, j, ldt)  (BLKADDR(A, PLASMA_Complex32_t, m, n)+((j)*(ldt)+(i)))
/***************************************************************************//**
 *
 **/
void plasma_pclansy(plasma_context_t *plasma)
{
    PLASMA_enum norm;
    PLASMA_enum uplo;
    PLASMA_desc A;
    float *work;
    float *result;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int m, n;
    int next_m;
    int next_n;
    int ldam, ldan;
    int step, lrank;
    int X, X1, X2, Y, Y1, Y2;

    float* lwork;
    float normtmp, normtmp2;

    plasma_unpack_args_7(norm, uplo, A, work, result, sequence, request);
    *result = 0.0;

    if (PLASMA_RANK == 0)
        memset(work, 0, PLASMA_SIZE*sizeof(float));
    ss_init(PLASMA_SIZE, 1, 0);

    switch (norm) {
    /*
     *  PlasmaMaxNorm
     */
    case PlasmaMaxNorm:
        n = 0;
        m = PLASMA_RANK;
        while (m >= A.mt && n < A.nt) {
            n++;
            m = m-A.mt+n;
        }

        while (n < A.nt) {
            next_m = m;
            next_n = n;

            next_m += PLASMA_SIZE;
            while (next_m >= A.mt && next_n < A.nt) {
                next_n++;
                next_m = next_m-A.mt+next_n;
            }

            if (m == n) {
                X1 = m == 0      ?  A.i       %A.mb   : 0;
                X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
                X = X2 - X1;

                ldam = BLKLDD(A, m);
                CORE_clansy(PlasmaMaxNorm, uplo, X, A(m, n, X1, X1, ldam), ldam, NULL, &normtmp);
            }
            else {
                /*
                 *  PlasmaLower
                 */
                if (uplo == PlasmaLower) {
                    X1 = m == 0      ?  A.i       %A.mb   : 0;
                    X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
                    X = X2 - X1;

                    Y1 = n == 0      ?  A.j       %A.nb   : 0;
                    Y2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                    Y = Y2 - Y1;

                    ldam = BLKLDD(A, m);
                    CORE_clange(PlasmaMaxNorm, X, Y, A(m, n, X1, Y1, ldam), ldam, NULL, &normtmp);
                }
                /*
                 *  PlasmaUpper
                 */
                else {
                    X1 = n == 0      ?  A.i       %A.mb   : 0;
                    X2 = n == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
                    X = X2 - X1;

                    Y1 = m == 0      ?  A.j       %A.nb   : 0;
                    Y2 = m == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                    Y = Y2 - Y1;

                    ldan = BLKLDD(A, n);
                    CORE_clange(PlasmaMaxNorm, X, Y, A(n, m, X1, Y1, ldan), ldan, NULL, &normtmp);
                }
            }
            if (normtmp > work[PLASMA_RANK])
                work[PLASMA_RANK] = normtmp;

            m = next_m;
            n = next_n;
        }
        ss_cond_set(PLASMA_RANK, 0, 1);
        break;
    /*
     *  PlasmaOneNorm / PlasmaInfNorm
     */
    case PlasmaOneNorm:
    case PlasmaInfNorm:
        m = PLASMA_RANK;
        normtmp2 = 0.0;
        lwork = (float*)plasma_private_alloc(plasma, A.mb, PlasmaRealDouble);

        while (m < A.mt) {
            X1 = m == 0      ?  A.i       %A.mb   : 0;
            X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
            X = X2 - X1;

            ldam = BLKLDD(A, m);
            memset(lwork, 0, A.mb*sizeof(float));
            /*
             *  PlasmaLower
             */
            if (uplo == PlasmaLower) {
                for (n = 0; n < m; n++) {
                    Y1 = n == 0 ? A.j%A.nb : 0;
                    Y = A.nb - Y1;
                    CORE_scasum(PlasmaRowwise, PlasmaUpperLower, X, Y, A(m, n, X1, Y1, ldam), ldam, lwork);
                }
                CORE_scasum(PlasmaRowwise, uplo, X, X, A(m, m, X1, X1, ldam), ldam, lwork);

                for (n = m+1; n < A.mt; n++) {
                    Y = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                    ldan = BLKLDD(A, n);
                    CORE_scasum(PlasmaColumnwise, PlasmaUpperLower, Y, X, A(n, m, 0, X1, ldan), ldan, lwork);
                }
            }
            /*
             *  PlasmaUpper
             */
            else {
                for (n = 0; n < m; n++) {
                    Y1 = n == 0 ?  A.j%A.nb : 0;
                    Y = A.nb - Y1;
                    CORE_scasum(PlasmaColumnwise, PlasmaUpperLower, Y, X, A(n, m, Y1, X1, A.nb), A.nb, lwork);
                }
                CORE_scasum(PlasmaRowwise, uplo, X, X, A(m, m, X1, X1, ldam), ldam, lwork);

                for ( n =m+1; n < A.mt; n++) {
                    Y = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                    CORE_scasum(PlasmaRowwise, PlasmaUpperLower, X, Y, A(m, n, X1, 0, ldam), ldam, lwork);
                }
            }
            CORE_slange(PlasmaMaxNorm, X, 1, lwork, 1, NULL, &normtmp);
            if (normtmp > normtmp2)
                normtmp2 = normtmp;

            m += PLASMA_SIZE;
        }
        work[PLASMA_RANK] = normtmp2;
        ss_cond_set(PLASMA_RANK, 0, 1);
        plasma_private_free(plasma, lwork);
        break;
    /*
     *  PlasmaFrobeniusNorm
     */
    case PlasmaFrobeniusNorm:
    default:;
    }

    if (norm != PlasmaFrobeniusNorm) {
        step = 1;
        lrank = PLASMA_RANK;
        while ( (lrank%2 == 0) && (PLASMA_RANK+step < PLASMA_SIZE) ) {
            ss_cond_wait(PLASMA_RANK+step, 0, step);
            work[PLASMA_RANK] = max(work[PLASMA_RANK], work[PLASMA_RANK+step]);
            lrank = lrank >> 1;
            step  = step << 1;
            ss_cond_set(PLASMA_RANK, 0, step);
        }
        if (PLASMA_RANK > 0) {
            while( lrank != 0 ) {
                if (lrank%2 == 1) {
                    ss_cond_set(PLASMA_RANK, 0, step);
                    lrank = 0;
                } else {
                    lrank = lrank >> 1;
                    step  = step << 1;
                    ss_cond_set(PLASMA_RANK, 0, step);
                }
            }
        }

        if (PLASMA_RANK == 0)
            *result = work[0];
    }
    ss_finalize();
}

/***************************************************************************//**
 *
 **/
void plasma_pclansy_quark(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc A, float *work, float *result,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int X, X1, X2, Y, Y1;
    int ldam;
    int m, n;
    int szeW, pos;
    float* lwork;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    *result = 0.0;
    switch ( norm ) {
    /*
     *  PlasmaMaxNorm
     */
    case PlasmaMaxNorm:
        szeW = A.mt*(A.mt+1)/2;
        pos = 0;
        lwork = (float*)plasma_shared_alloc(plasma, szeW, PlasmaRealDouble);
        memset(lwork, 0, szeW*sizeof(float));
        for(m = 0; m < A.mt; m++) {
            X1 = m == 0      ?  A.i       %A.mb   : 0;
            X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
            X = X2 - X1;
            ldam = BLKLDD(A, m);
            QUARK_CORE_clansy_f1(
                plasma->quark, &task_flags,
                PlasmaMaxNorm, uplo, X,
                A(m, m, X1, X1, ldam), ldam, ldam*X,
                0, &(lwork[pos]),
                lwork, szeW);
            pos++;
            /*
             *  PlasmaLower
             */
            if (uplo == PlasmaLower) {
                for(n=0; n<m; n++) {
                    Y1 = n == 0 ? A.j%A.nb : 0;
                    Y  = A.nb - Y1;
                    QUARK_CORE_clange_f1(
                        plasma->quark, &task_flags,
                        PlasmaMaxNorm, X, Y,
                        A(m, n, X1, Y1, ldam), ldam, ldam*Y,
                        0, &(lwork[pos]),
                        lwork, szeW);
                    pos++;
                }
            }
            /*
             *  PlasmaUpper
             */
            else {
                for(n = m+1; n < A.mt; n++) {
                    Y = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                    QUARK_CORE_clange_f1(
                        plasma->quark, &task_flags,
                        PlasmaMaxNorm, X, Y,
                        A(m, n, X1, 0, ldam), ldam, ldam*Y,
                        0, &(lwork[pos]),
                        lwork, szeW);
                    pos++;
                }
            }
        }
        QUARK_CORE_slange(
            plasma->quark, &task_flags,
            PlasmaMaxNorm, szeW, 1,
            lwork, 1, szeW,
            0, result);

        QUARK_CORE_free(plasma->quark, &task_flags, lwork, szeW*sizeof(PLASMA_Complex32_t));
        break;
    /*
     *  PlasmaOneNorm / PlasmaInfNorm
     */
    case PlasmaOneNorm:
    case PlasmaInfNorm:
        lwork = (float *)plasma_shared_alloc(plasma, A.m+1, PlasmaRealDouble);
        memset(lwork, 0, (A.m+1)*sizeof(float));
        for(m = 0; m < A.mt; m++) {
            X1 = m == 0      ?  A.i       %A.mb   : 0;
            X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
            X = X2 - X1;
            ldam = BLKLDD(A, m);
            QUARK_CORE_scasum_f1(
                plasma->quark, &task_flags,
                PlasmaRowwise, uplo, X, X,
                A(m, m, X1, X1, ldam), ldam, ldam*X,
                &(lwork[m*A.mb+1]), A.mb,
                lwork, A.m);
            /*
             *  PlasmaLower
             */
            if (uplo == PlasmaLower) {
                for(n = 0; n < m; n++) {
                    Y1 = n == 0 ? A.j%A.nb : 0;
                    Y = A.nb - Y1;
                    QUARK_CORE_scasum_f1(
                        plasma->quark, &task_flags,
                        PlasmaRowwise, PlasmaUpperLower, X, Y,
                        A(m, n, X1, Y1, ldam), ldam, ldam*Y,
                        &(lwork[m*A.mb+1]), A.mb,
                        lwork, A.m);
                    QUARK_CORE_scasum_f1(
                        plasma->quark, &task_flags,
                        PlasmaColumnwise, PlasmaUpperLower, X, Y,
                        A(m, n, X1, Y1, ldam), ldam, ldam*Y,
                        &(lwork[n*A.mb+1]), A.mb,
                        lwork, A.m);
                }
            }
            /*
             *  PlasmaUpper
             */
            else {
                for(n = m+1; n < A.mt; n++) {
                    Y = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                    QUARK_CORE_scasum_f1(
                        plasma->quark, &task_flags,
                        PlasmaRowwise, PlasmaUpperLower, X, Y,
                        A(m, n, X1, 0, ldam), ldam, ldam*Y,
                        &(lwork[m*A.mb+1]), A.mb,
                        lwork, A.m);
                    QUARK_CORE_scasum_f1(
                        plasma->quark, &task_flags,
                        PlasmaColumnwise, PlasmaUpperLower, X, Y,
                        A(m, n, X1, 0, ldam), ldam, ldam*Y,
                        &(lwork[n*A.mb+1]), A.mb,
                        lwork, A.m);
                }
            }
        }
        QUARK_CORE_slange(
            plasma->quark, &task_flags,
            PlasmaMaxNorm, A.m+1, 1,
            lwork, 1, A.m+1,
            0, result);

        QUARK_CORE_free(plasma->quark, &task_flags, lwork, (A.m+1)*sizeof(PLASMA_Complex32_t));
        break;
    /*
     *  PlasmaFrobeniusNorm - not implemented
     */
    case PlasmaFrobeniusNorm:
    default:;
    }
}
