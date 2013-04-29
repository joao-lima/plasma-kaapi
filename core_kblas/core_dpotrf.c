/*
 -- Innovative Computing Laboratory
 -- Electrical Engineering and Computer Science Department
 -- University of Tennessee
 -- (C) Copyright 2008-2010
 
 Redistribution  and  use  in  source and binary forms, with or without
 modification,  are  permitted  provided  that the following conditions
 are met:
 
 * Redistributions  of  source  code  must  retain  the above copyright
 notice,  this  list  of  conditions  and  the  following  disclaimer.
 * Redistributions  in  binary  form must reproduce the above copyright
 notice,  this list of conditions and the following disclaimer in the
 documentation  and/or other materials provided with the distribution.
 * Neither  the  name of the University of Tennessee, Knoxville nor the
 names of its contributors may be used to endorse or promote products
 derived from this software without specific prior written permission.
 
 THIS  SOFTWARE  IS  PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 ``AS IS''  AND  ANY  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A  PARTICULAR  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL,  EXEMPLARY,  OR  CONSEQUENTIAL  DAMAGES  (INCLUDING,  BUT NOT
 LIMITED  TO,  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA,  OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 THEORY  OF  LIABILITY,  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF  THIS  SOFTWARE,  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <lapacke.h>

#include "core_kblas.h"

#include "common.h"

void CORE_dpotrf_cpu(int uplo, int N, double *A, int LDA, int *INFO)
{
#if 0
  fprintf(stdout, "%s: uplo=%d N=%d A=%p LDA=%d\n", __FUNCTION__,
          uplo, N, A, LDA);
  fflush(stdout);
#endif
  *INFO = LAPACKE_dpotrf_work(
                              LAPACK_COL_MAJOR,
                              lapack_const(uplo),
                              N, A, LDA );
}

#define   A(i,j)    ((double*)&A[(i*blocksize + (j*blocksize*LDA))])
int CORE_dpotrf_parallel_cpu(int uplo, int N, double *A, int LDA, int *INFO)
{
  int k, m, n;
  int ldak, ldam;
  int tempkm, tempmm;
  
  const int blocksize = 64;
  const int nb = (N+blocksize-1) / 64;
  
  double zone  = (double) 1.0;
  double mzone = (double)-1.0;
  
  if(n > 500) /* from control/context.c */
  {
    if (uplo == PlasmaLower) {
      for (k = 0; k < nb; k++) {
        tempkm = k == nb-1 ? N - k*blocksize : blocksize;
//        ldak = BLKLDD(A, k);
        ldak = LDA;
        CORE_dpotrf(PlasmaLower, tempkm, A(k, k), ldak, INFO);
        
        for (m = k+1; m < nb; m++) {
//          tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
          tempmm = m == nb-1 ? m - m*blocksize : blocksize;
//          ldam = BLKLDD(A, m);
          ldam = LDA;
          CORE_dtrsm(PlasmaRight, PlasmaLower, PlasmaTrans, PlasmaNonUnit,
                           tempmm, blocksize,
                           zone, A(k, k), ldak,
                           A(m, k), ldam);
        }
        for (m = k+1; m < nb; m++) {
//          tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
          tempmm = m == nb-1 ? m - m*blocksize : blocksize;
//          ldam = BLKLDD(A, m);
          ldam = LDA;
          CORE_dsyrk(PlasmaLower, PlasmaNoTrans,
                           tempmm, blocksize,
                           -1.0, A(m, k), ldam,
                           1.0, A(m, m), ldam);
          
          for (n = k+1; n < m; n++) {
            CORE_dgemm( PlasmaNoTrans, PlasmaTrans,
                             tempmm, blocksize, blocksize,
                             mzone, A(m, k), ldam,
                             A(n, k), blocksize,
                             zone,  A(m, n), ldam);
          }
        }
      }
    }
  }
  else
  {
    CORE_dpotrf_cpu(uplo, N, A, LDA, INFO);
  }
  
  return 0;
}

void CORE_dpotrf_quark_cpu(Quark *quark)
{
  int uplo;
  int n;
  double *A;
  int lda;
  PLASMA_sequence *sequence;
  PLASMA_request *request;
  int iinfo;
  
  int info;
  
  quark_unpack_args_7(quark, uplo, n, A, lda, sequence, request, iinfo);

  CORE_dpotrf_parallel_cpu(uplo, n, A, lda, &info);
  
  if (sequence->status == PLASMA_SUCCESS && info != 0)
    plasma_sequence_flush(quark, sequence, request, iinfo+info);
}
