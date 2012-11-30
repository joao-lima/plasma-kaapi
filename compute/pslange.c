/**
 *
 * @file pslange.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Emmanuel Agullo
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Thu Sep 15 12:09:21 2011
 *
 **/
#include <stdlib.h>
#include <math.h>
#include "common.h"

#define A(m, n, i, j, ldt)  (BLKADDR(A, float, m, n)+((j)*(ldt)+(i)))
/***************************************************************************//**
 *
 **/
void plasma_pslange(plasma_context_t *plasma)
{
    PLASMA_enum norm;
    PLASMA_desc A;
    float *work;
    float *result;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int m, n;
    int next_m;
    int next_n;
    int ldam;
    int step, lrank;
    int X, X1, X2, Y, Y1, Y2;

    float* lwork;
    float normtmp, normtmp2;

    plasma_unpack_args_6(norm, A, work, result, sequence, request);
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
            m = m-A.mt;
        }

        while (n < A.nt) {
            next_m = m;
            next_n = n;

            next_m += PLASMA_SIZE;
            while (next_m >= A.mt && next_n < A.nt) {
                next_n++;
                next_m = next_m-A.mt;
            }

            X1 = m == 0      ?  A.i       %A.mb   : 0;
            X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
            X = X2 - X1;

            Y1 = n == 0      ?  A.j       %A.nb   : 0;
            Y2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
            Y = Y2 - Y1;

            ldam = BLKLDD(A, m);
            CORE_slange(PlasmaMaxNorm, X, Y, A(m, n, X1, Y1, ldam), ldam, NULL, &normtmp);

            if (normtmp > work[PLASMA_RANK])
                work[PLASMA_RANK] = normtmp;

            m = next_m;
            n = next_n;
        }
        ss_cond_set(PLASMA_RANK, 0, 1);
        break;
    /*
     *  PlasmaOneNorm
     */
    case PlasmaOneNorm:

        n = PLASMA_RANK;
        normtmp2 = 0.0;
        lwork = (float*)plasma_private_alloc(plasma, A.nb, PlasmaRealDouble);

        while (n < A.nt) {
            Y1 = n == 0      ?  A.j       %A.nb   : 0;
            Y2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
            Y = Y2 - Y1;
            memset(lwork, 0, A.nb*sizeof(float));
            for (m = 0; m < A.mt; m++) {
                X1 = m == 0      ?  A.i       %A.mb   : 0;
                X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
                X = X2 - X1;

                ldam = BLKLDD(A, m);
                CORE_sasum(
                    PlasmaColumnwise, PlasmaUpperLower,
                    X, Y,
                    A(m, n, X1, Y1, ldam), ldam,
                    lwork);
            }
            CORE_slange(PlasmaMaxNorm, Y, 1, lwork, 1, NULL, &normtmp);
            if (normtmp > normtmp2)
                normtmp2 = normtmp;

            n += PLASMA_SIZE;
        }
        work[PLASMA_RANK] = normtmp2;
        ss_cond_set(PLASMA_RANK, 0, 1);
        plasma_private_free(plasma, lwork);
        break;
    /*
     *  PlasmaInfNorm
     */
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

            for (n = 0; n < A.nt; n++) {
                Y1 = n == 0      ?  A.j       %A.nb   : 0;
                Y2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                Y = Y2 - Y1;
                CORE_sasum(
                    PlasmaRowwise, PlasmaUpperLower,
                    X, Y,
                    A(m, n, X1, Y1, ldam), ldam,
                    lwork);
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
     *  PlasmaFrobeniusNorm - not implemented
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
void plasma_pslange_quark(PLASMA_enum norm, PLASMA_desc A, float *work, float *result,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int X, X1, X2, Y, Y1, Y2;
    int ldam;
    int m, n;
    int szeW;
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
        szeW = A.mt*A.nt;
        lwork = (float*)plasma_shared_alloc(plasma, szeW, PlasmaRealDouble);
        memset(lwork, 0, szeW*sizeof(float));
        for(m = 0; m < A.mt; m++) {
            X1 = m == 0      ?  A.i       %A.mb   : 0;
            X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
            X = X2 - X1;
            ldam = BLKLDD(A, m);
            for(n = 0; n < A.nt; n++) {
                Y1 = n == 0      ?  A.j       %A.nb   : 0;
                Y2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                Y = Y2 - Y1;
                QUARK_CORE_slange_f1(
                    plasma->quark, &task_flags,
                    PlasmaMaxNorm, X, Y,
                    A(m, n, X1, Y1, ldam), ldam, ldam*Y,
                    0, &(lwork[A.mt*n+m]),
                    lwork, szeW);
            }
        }
        QUARK_CORE_slange(
            plasma->quark, &task_flags,
            PlasmaMaxNorm, A.mt, A.nt,
            lwork, A.mt, szeW,
            0, result);

        QUARK_CORE_free(plasma->quark, &task_flags, lwork, szeW*sizeof(float));
        break;

    /*
     *  PlasmaOneNorm
     */
    case PlasmaOneNorm:
        lwork = (float*)plasma_shared_alloc(plasma, (A.n+1), PlasmaRealDouble);
        memset(lwork, 0, (A.n+1)*sizeof(float));
        for(m = 0; m < A.mt; m++) {
            X1 = m == 0      ?  A.i       %A.mb   : 0;
            X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
            X = X2 - X1;
            ldam = BLKLDD(A, m);
            for(n = 0; n < A.nt; n++) {
                Y1 = n == 0      ?  A.j       %A.nb   : 0;
                Y2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                Y = Y2 - Y1;
                QUARK_CORE_sasum_f1(
                    plasma->quark, &task_flags,
                    PlasmaColumnwise, PlasmaUpperLower, X, Y,
                    A(m, n, X1, Y1, ldam), ldam, ldam*Y,
                    &(lwork[n*A.nb+1]), A.nb,
                    lwork, A.n);
            }
        }
        QUARK_CORE_slange(
            plasma->quark, &task_flags,
            PlasmaMaxNorm, A.n+1, 1,
            lwork, 1, A.n+1,
            0, result);

        QUARK_CORE_free(plasma->quark, &task_flags, lwork, (A.n+1)*sizeof(float));
        break;
    /*
     *  PlasmaInfNorm
     */
    case PlasmaInfNorm:
        lwork = (float*)plasma_shared_alloc(plasma, (A.m+1), PlasmaRealDouble);
        memset(lwork, 0, (A.m+1)*sizeof(float));
        for(m = 0; m < A.mt; m++) {
            X1 = m == 0      ?  A.i       %A.mb   : 0;
            X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
            X = X2 - X1;
            ldam = BLKLDD(A, m);
            for(n = 0; n < A.nt; n++) {
                Y1 = n == 0      ?  A.j       %A.nb   : 0;
                Y2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                Y = Y2 - Y1;
                QUARK_CORE_sasum_f1(
                    plasma->quark, &task_flags,
                    PlasmaRowwise, PlasmaUpperLower, X, Y,
                    A(m, n, X1, Y1, ldam), ldam, ldam*Y,
                    &(lwork[m*A.mb+1]), A.mb,
                    lwork, A.m);
            }
        }
        QUARK_CORE_slange(
            plasma->quark, &task_flags,
            PlasmaMaxNorm, A.m+1, 1,
            lwork, 1, A.m+1,
            0, result);

        QUARK_CORE_free(plasma->quark, &task_flags, lwork, (A.m+1)*sizeof(float));
        break;
    /*
     *  PlasmaFrobeniusNorm - not implemented
     */
    case PlasmaFrobeniusNorm:
    default:;
    }
}