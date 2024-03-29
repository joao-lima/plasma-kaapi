/**
 *
 * @file pdgerbb.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Hatem Ltaief
 * @author Azzam Haidar
 * @date 2010-11-15
 * @generated d Thu Sep 15 12:09:24 2011
 *
 **/
#include "common.h"

/***************************************************************************//**
 *  Parallel tile BAND Bidiagonal Reduction - dynamic scheduler
 *  Could be optimized by using the algorithms from Trefethen book
 *
 * WARNING: do never call this function because ormqr and unmlq are
 * not implementing all the cases required in static.
 *
 **/
void plasma_pdgerbb(plasma_context_t *plasma)
{
    PLASMA_desc A;
    PLASMA_desc T;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int k;
    int tempkm, tempkn;

    plasma_unpack_args_4(A, T, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;

    if (A.m >= A.n){
       for (k = 0; k < A.nt; k++) {
           tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
           tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
  
           plasma_static_call_4(plasma_pdgeqrf,
               PLASMA_desc, plasma_desc_submatrix(A, k*A.mb, k*A.nb,  A.m-k*A.mb, tempkn),
               PLASMA_desc, plasma_desc_submatrix(T, k*T.mb, k*T.nb,  T.m-k*T.mb, tempkn),
               PLASMA_sequence*, sequence,
               PLASMA_request*, request);
  
           plasma_static_call_7(plasma_pdormqr,
               PLASMA_enum, PlasmaLeft,
               PLASMA_enum, PlasmaTrans,
               PLASMA_desc, plasma_desc_submatrix(A, k*A.mb, k*A.nb,  A.m-k*A.mb, tempkn),
               PLASMA_desc, plasma_desc_submatrix(A, k*A.mb, (k+1)*A.nb,  A.m-k*A.mb, A.n-(k+1)*A.nb),
               PLASMA_desc, plasma_desc_submatrix(T, k*T.mb, k*T.nb,  T.m-k*T.mb, tempkn),
               PLASMA_sequence*, sequence,
               PLASMA_request*, request);
  
           if (k+1 < A.nt){
              tempkn = k+1 == A.nt-1 ? A.n-(k+1)*A.nb : A.nb;
  
              plasma_static_call_4(plasma_pdgelqf,
                  PLASMA_desc, plasma_desc_submatrix(A, k*A.mb, (k+1)*A.nb, tempkm, A.n-(k+1)*A.nb),
                  PLASMA_desc, plasma_desc_submatrix(T, k*T.mb, (k+1)*T.nb, tempkm, T.n-(k+1)*T.nb),
                  PLASMA_sequence*, sequence,
                  PLASMA_request*, request);
  
              plasma_static_call_7(plasma_pdormlq,
                  PLASMA_enum, PlasmaRight,
                  PLASMA_enum, PlasmaTrans,
                  PLASMA_desc, plasma_desc_submatrix(A, k*A.mb, (k+1)*A.nb, tempkm, A.n-(k+1)*A.nb),
                  PLASMA_desc, plasma_desc_submatrix(A, (k+1)*A.mb, (k+1)*A.nb, A.m-(k+1)*A.mb, A.n-(k+1)*A.nb),
                  PLASMA_desc, plasma_desc_submatrix(T, k*T.mb, (k+1)*T.nb, tempkm, T.n-(k+1)*T.nb),
                  PLASMA_sequence*, sequence,
                  PLASMA_request*, request);
           }
       }
    }
    else{
       for (k = 0; k < A.mt; k++) {
           tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
           tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
   
           plasma_static_call_4(plasma_pdgelqf,
               PLASMA_desc, plasma_desc_submatrix(A, k*A.mb, k*A.nb, tempkm, A.n-k*A.nb),
               PLASMA_desc, plasma_desc_submatrix(T, k*T.mb, k*T.nb, tempkm, T.n-k*T.nb),
               PLASMA_sequence*, sequence,
               PLASMA_request*, request);
   
           plasma_static_call_7(plasma_pdormlq,
               PLASMA_enum, PlasmaRight,
               PLASMA_enum, PlasmaTrans,
               PLASMA_desc, plasma_desc_submatrix(A, k*A.mb, k*A.nb, tempkm, A.n-k*A.nb),
               PLASMA_desc, plasma_desc_submatrix(A, (k+1)*A.mb, k*A.nb, A.m-(k+1)*A.mb, A.n-k*A.nb),
               PLASMA_desc, plasma_desc_submatrix(T, k*T.mb, k*T.nb, tempkm, T.n-k*T.nb),
               PLASMA_sequence*, sequence,
               PLASMA_request*, request);
   
           if (k+1 < A.mt){
              tempkm = k+1 == A.mt-1 ? A.m-(k+1)*A.mb : A.mb;
              tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
   
              plasma_static_call_4(plasma_pdgeqrf,
                  PLASMA_desc, plasma_desc_submatrix(A, (k+1)*A.mb, k*A.nb,  A.m-(k+1)*A.mb, tempkn),
                  PLASMA_desc, plasma_desc_submatrix(T, (k+1)*T.mb, k*T.nb,  T.m-(k+1)*T.mb, tempkn),
                  PLASMA_sequence*, sequence,
                  PLASMA_request*, request);
       
              plasma_static_call_7(plasma_pdormqr,
                  PLASMA_enum, PlasmaLeft,
                  PLASMA_enum, PlasmaTrans,
                  PLASMA_desc, plasma_desc_submatrix(A, (k+1)*A.mb, k*A.nb,  A.m-(k+1)*A.mb, tempkn),
                  PLASMA_desc, plasma_desc_submatrix(A, (k+1)*A.mb, (k+1)*A.nb,  A.m-(k+1)*A.mb, A.n-(k+1)*A.nb),
                  PLASMA_desc, plasma_desc_submatrix(T, (k+1)*T.mb, k*T.nb,  T.m-(k+1)*T.mb, tempkn),
                  PLASMA_sequence*, sequence,
                  PLASMA_request*, request);
           }
       }
    }
}

/***************************************************************************//**
 *  Parallel tile BAND Bidiagonal Reduction - dynamic scheduler
 *  Could be optimized by using the algorithms from Trefethen book
 **/
void plasma_pdgerbb_quark(PLASMA_desc A, PLASMA_desc T,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k;
    int tempkm, tempkn;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    if (A.m >= A.n){
       for (k = 0; k < A.nt; k++) {
           tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
           tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
  
           plasma_pdgeqrf_quark(
               plasma_desc_submatrix(A, k*A.mb, k*A.nb,  A.m-k*A.mb, tempkn),
               plasma_desc_submatrix(T, k*T.mb, k*T.nb,  T.m-k*T.mb, tempkn),
               sequence, request);
  
           plasma_pdormqr_quark(
               PlasmaLeft,
               PlasmaTrans,
               plasma_desc_submatrix(A, k*A.mb, k*A.nb,  A.m-k*A.mb, tempkn),
               plasma_desc_submatrix(A, k*A.mb, (k+1)*A.nb,  A.m-k*A.mb, A.n-(k+1)*A.nb),
               plasma_desc_submatrix(T, k*T.mb, k*T.nb,  T.m-k*T.mb, tempkn),
               sequence, request);
  
           if (k+1 < A.nt){
              tempkn = k+1 == A.nt-1 ? A.n-(k+1)*A.nb : A.nb;
  
              plasma_pdgelqf_quark(
                  plasma_desc_submatrix(A, k*A.mb, (k+1)*A.nb, tempkm, A.n-(k+1)*A.nb),
                  plasma_desc_submatrix(T, k*T.mb, (k+1)*T.nb, tempkm, T.n-(k+1)*T.nb),
                  sequence, request);
  
              plasma_pdormlq_quark(
                  PlasmaRight, PlasmaTrans,
                  plasma_desc_submatrix(A, k*A.mb, (k+1)*A.nb, tempkm, A.n-(k+1)*A.nb),
                  plasma_desc_submatrix(A, (k+1)*A.mb, (k+1)*A.nb, A.m-(k+1)*A.mb, A.n-(k+1)*A.nb),
                  plasma_desc_submatrix(T, k*T.mb, (k+1)*T.nb, tempkm, T.n-(k+1)*T.nb),
                  sequence, request);
           }
       }
    }
    else{
       for (k = 0; k < A.mt; k++) {
           tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
           tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
   
           plasma_pdgelqf_quark(
               plasma_desc_submatrix(A, k*A.mb, k*A.nb, tempkm, A.n-k*A.nb),
               plasma_desc_submatrix(T, k*T.mb, k*T.nb, tempkm, T.n-k*T.nb),
               sequence, request);
   
           plasma_pdormlq_quark(
               PlasmaRight, PlasmaTrans,
               plasma_desc_submatrix(A, k*A.mb, k*A.nb, tempkm, A.n-k*A.nb),
               plasma_desc_submatrix(A, (k+1)*A.mb, k*A.nb, A.m-(k+1)*A.mb, A.n-k*A.nb),
               plasma_desc_submatrix(T, k*T.mb, k*T.nb, tempkm, T.n-k*T.nb),
               sequence, request);
           
           if (k+1 < A.mt){
              tempkm = k+1 == A.mt-1 ? A.m-(k+1)*A.mb : A.mb;
              tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
   
              plasma_pdgeqrf_quark(
                   plasma_desc_submatrix(A, (k+1)*A.mb, k*A.nb,  A.m-(k+1)*A.mb, tempkn),
                   plasma_desc_submatrix(T, (k+1)*T.mb, k*T.nb,  T.m-(k+1)*T.mb, tempkn),
                   sequence, request);
       
              plasma_pdormqr_quark(
                  PlasmaLeft, PlasmaTrans,
                  plasma_desc_submatrix(A, (k+1)*A.mb, k*A.nb,  A.m-(k+1)*A.mb, tempkn),
                  plasma_desc_submatrix(A, (k+1)*A.mb, (k+1)*A.nb,  A.m-(k+1)*A.mb, A.n-(k+1)*A.nb),
                  plasma_desc_submatrix(T, (k+1)*T.mb, k*T.nb,  T.m-(k+1)*T.mb, tempkn),
                  sequence, request);
           }
       }
    }
}
