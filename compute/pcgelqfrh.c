/**
 *
 * @file pcgelqfrh.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Dulceneia Becker
 * @date 2010-11-15
 * @generated c Thu Sep 15 12:09:20 2011
 *
 **/
#include "common.h"

#define A(m,n)  BLKADDR(A, PLASMA_Complex32_t, (m), (n))
#define T(m,n)  BLKADDR(T, PLASMA_Complex32_t, (m), (n))
#define T2(m,n) BLKADDR(T, PLASMA_Complex32_t, (m), (n)+A.nt)
/***************************************************************************//**
 *  Parallel tile LQ factorization (reduction Householder) - static / sequential
 **/
void plasma_pcgelqfrh(plasma_context_t *plasma)
{
    PLASMA_desc A;
    PLASMA_desc T;
    int BS;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int k, m, n;
    int K, N, RD;
    int ldak, ldam;
    int tempkm, tempNn, tempmm, tempnn, tempNRDn;
    int ib;

    if (PLASMA_RANK != 0) return;

    plasma_unpack_args_5(A, T, BS, sequence, request);
    ib = PLASMA_IB;

    PLASMA_Complex32_t *work, *tau;
    work = (PLASMA_Complex32_t*)plasma_private_alloc(plasma, ib*T.nb, T.dtyp);
    tau  = (PLASMA_Complex32_t*)plasma_private_alloc(plasma, A.nb, A.dtyp);

    if (sequence->status != PLASMA_SUCCESS)
        return;

    K = min(A.mt, A.nt);
    for (k = 0; k < K; k++) {
        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
        ldak = BLKLDD(A, k);
        for (N = k;
             N < A.nt-1 || N == k;  // No rightmost single-column subdomain
             N += BS) {
            tempNn = N == A.nt-1 ? A.n-N*A.nb : A.nb;
            CORE_cgelqt(
                tempkm, tempNn, ib,
                A(k, N), ldak,
                T(k, N), T.mb,
                tau, work);

            for (m = k+1; m < A.mt; m++) {
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                ldam = BLKLDD(A, m);
                CORE_cunmlq(
                    PlasmaRight, PlasmaConjTrans,
                    tempmm, tempNn, tempNn, ib,
                    A(k, N), ldak,
                    T(k, N), T.mb,
                    A(m, N), ldam,
                    work , A.nb);
            }
            for (n = N+1;
                 (n < N+BS && n < A.nt) || n == A.nt-1; // Suck in rightmost single-column domain
                 n++) {
                tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                CORE_ctslqt(
                    tempkm, tempnn, ib,
                    A(k, N), ldak,
                    A(k, n), ldak,
                    T(k, n), T.mb,
                    tau, work);

                for (m = k+1; m < A.mt; m++) {
                    tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                    ldam = BLKLDD(A, m);
                    CORE_ctsmlq(
                        PlasmaRight, PlasmaConjTrans,
                        tempmm, A.mb, tempmm, tempnn, A.mb, ib,
                        A(m, N), ldam,
                        A(m, n), ldam,
                        A(k, n), ldak,
                        T(k, n), T.mb,
                        work , A.nb);
                }
            }
        }
        for (RD = BS; RD < A.nt-k; RD *= 2) {
            for (N = k;
                 N+RD < A.nt-1; // No reduction with rightmost single-column subdomain
                 N += 2*RD) {
                tempNRDn = N+RD == A.nt-1 ? A.n-(N+RD)*A.nb : A.nb;
                CORE_cttlqt(
                    tempkm, tempNRDn, ib,
                    A (k, N   ), ldak,
                    A (k, N+RD), ldak,
                    T2(k, N+RD), T.mb,
                    tau, work);

                for (m = k+1; m < A.mt; m++) {
                    tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                    ldam = BLKLDD(A, m);
                    CORE_cttmlq(
                        PlasmaRight, PlasmaConjTrans,
                        tempmm, A.nb, tempmm, tempNRDn, tempkm, ib,
                        A (m, N   ), ldam,
                        A (m, N+RD), ldam,
                        A (k, N+RD), ldak,
                        T2(k, N+RD), T.mb,
                        work , A.nb);
                }
            }
        }
    }
}

/***************************************************************************//**
 *  Parallel tile LQ factorization (reduction Householder) - dynamic scheduling
 **/
void plasma_pcgelqfrh_quark(PLASMA_desc A, PLASMA_desc T, int BS,
                            PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int K, N, RD;
    int ldak, ldam;
    int tempkm, tempNn, tempmm, tempnn, tempNRDn;
    int ib;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    ib = PLASMA_IB;
    K = min(A.mt, A.nt);
    for (k = 0; k < K; k++) {

        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
        ldak = BLKLDD(A, k);
//      for (N = k; N < A.nt; N += BS) {
        for (N = k;
             N < A.nt-1 || N == k;  // No rightmost single-column subdomain
             N += BS) {
            tempNn = N == A.nt-1 ? A.n-N*A.nb : A.nb;
            QUARK_CORE_cgelqt(
                plasma->quark, &task_flags,
                tempkm, tempNn, ib, T.nb,
                A(k, N), ldak,
                T(k, N), T.mb);

            for (m = k+1; m < A.mt; m++) {
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                ldam = BLKLDD(A, m);
                QUARK_CORE_cunmlq(
                    plasma->quark, &task_flags,
                    PlasmaRight, PlasmaConjTrans,
                    tempmm, tempNn, tempNn, ib, T.nb,
                    A(k, N), ldak,
                    T(k, N), T.mb,
                    A(m, N), ldam);
            }
//          for (n = N+1; n < N+BS && n < A.nt; n++) {
            for (n = N+1;
                 (n < N+BS && n < A.nt) || n == A.nt-1; // Suck in rightmost single-column domain
                 n++) {
                tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                QUARK_CORE_ctslqt(
                    plasma->quark, &task_flags,
                    tempkm, tempnn, ib, T.nb,
                    A(k, N), ldak,
                    A(k, n), ldak,
                    T(k, n), T.mb);

                for (m = k+1; m < A.mt; m++) {
                    tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                    ldam = BLKLDD(A, m);
                    QUARK_CORE_ctsmlq(
                        plasma->quark, &task_flags,
                        PlasmaRight, PlasmaConjTrans,
                        tempmm, A.nb, tempmm, tempnn, A.mb, ib, T.nb,
                        A(m, N), ldam,
                        A(m, n), ldam,
                        A(k, n), ldak,
                        T(k, n), T.mb);
                }
            }
        }
        for (RD = BS; RD < A.nt-k; RD *= 2) {
//          for (N = k; N+RD < A.nt; N += 2*RD) {
            for (N = k;
                 N+RD < A.nt-1; // No reduction with rightmost single-column subdomain
                 N += 2*RD) {
                tempNRDn = N+RD == A.nt-1 ? A.n-(N+RD)*A.nb : A.nb;
                QUARK_CORE_cttlqt(
                    plasma->quark, &task_flags,
                    tempkm, tempNRDn, ib, T.nb,
                    A (k, N   ), ldak,
                    A (k, N+RD), ldak,
                    T2(k, N+RD), T.mb);

                for (m = k+1; m < A.mt; m++) {
                    tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                    ldam = BLKLDD(A, m);
                    QUARK_CORE_cttmlq(
                        plasma->quark, &task_flags,
                        PlasmaRight, PlasmaConjTrans,
                        tempmm, A.nb, tempmm, tempNRDn, A.mb, ib, T.nb,
                        A (m, N   ), ldam,
                        A (m, N+RD), ldam,
                        A (k, N+RD), ldak,
                        T2(k, N+RD), T.mb);
                }
            }
        }
    }
}
