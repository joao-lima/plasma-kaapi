/**
 *
 * @file pzhbrdt.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.0
 * @author Azzam Haidar
 * @date 2011-05-15
 * @precisions normal z -> c d s
 *
 **/
#include <lapacke.h>
#include "common.h"

/*
 * struct timeval {time_t tv_sec; suseconds_t tv_usec;};
 */
#include <sys/time.h>
double Wtimming(void);
double Wtimming(void)
{
    struct timeval tp;
    gettimeofday( &tp, NULL );
    return tp.tv_sec + 1e-6 * tp.tv_usec;
}
    double tblg;



#undef REAL
#define COMPLEX
#define A(_m, _n) (PLASMA_Complex64_t *)plasma_geteltaddr(&A, (_m), (_n), eltsize)
/***************************************************************************//**
 *  Parallel Reduction from BAND tridiagonal to the final condensed form - dynamic scheduler
 **/
void plasma_pzhbrdt_quark(PLASMA_enum uplo,
                          PLASMA_desc A, double *D, double *E, PLASMA_desc T,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

#ifdef COMPLEX
    static PLASMA_Complex64_t zone  = (PLASMA_Complex64_t) 1.0;
    static double             dzero = (double) 0.0;
    PLASMA_Complex64_t ztmp;
    double absztmp;
#endif

    PLASMA_Complex64_t *C, *S;
    int blksweep, lcsweep, blkid, lcNB;
    int N, NB, NT, grsiz, lcgrsiz;
    int i;
    size_t eltsize = plasma_element_size(A.dtyp);

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;

    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    NT = A.nt;
    N  = A.m;
    NB = A.mb;

    /* Quick return */
    if (N == 0){
        return;
    }

    if (NB == 0) {
        memset(D, 0,     N*sizeof(double));
        memset(E, 0, (N-1)*sizeof(double));
#ifdef COMPLEX
        for (i=0; i<N; i++)
            D[i]  = cabs(*A(i,i));
#else
        for (i=0; i<N; i++)
          D[i]  = *A(i,i);
#endif
        return;
    }

    /*
     * Barrier is used because the bulge have to wait until
     * the reduction to band has been finish.
     * otherwise, I can remove this BARRIER when I integrate
     * the function dependencies link inside the reduction to
     * band. Keep in min the case when NB=1, where no bulge-chasing.
     */
    /***************************************************************/
    QUARK_Barrier(plasma->quark);
    tblg   = -Wtimming();
    /***************************************************************/

    /*
     * Case NB=1 ==> matrix is already Bidiagonal. no need to bulge.
     * Make diagonal and superdiagonal elements real, storing them in
     * D and E. if PlasmaLower, first transform lower bidiagonal form
     * to upper bidiagonal by applying plane rotations/ Householder
     * from the left, overwriting superdiagonal elements then make
     * elements real of the resulting upper Bidiagonal. if PlasmaUpper
     * then make its elements real.  For Q, PT: ZSCAL should be done
     * in case of WANTQ.
     */
    if (NB == 1){
        memset(D, 0, N    *sizeof(double));
        memset(E, 0, (N-1)*sizeof(double));
#ifdef COMPLEX
        if(uplo==PlasmaLower){
            for (i=0; i<N; i++)
            {
                D[i] = creal( *A(i, i) );               /* diag value */
                if( i < (N-1)) {                            /* lower off-diag value */
                    ztmp        = *A((i+1),i);
                    absztmp     = cabs(ztmp);
                    *A((i+1),i) = absztmp;
                    E[i]        = absztmp;
                    if(absztmp != dzero)
                        ztmp = (PLASMA_Complex64_t) (ztmp / absztmp);
                    else
                        ztmp = zone;
                    if(i<(N-2)) *A((i+2),(i+1)) = *A((i+2),(i+1)) * ztmp;
                    /* for Q: ZSCAL should be done in case of WANTQ */
                }
            }
        } else { /* PlasmaUpper */
            for (i=0; i<N; i++)
            {
                D[i]  =  creal( *A(i,i) );               /* diag value*/
                if(i<(N-1)) {                            /* lower off-diag value */
                    ztmp        = *A(i, (i+1));
                    absztmp     = cabs(ztmp);
                    *A(i,(i+1)) = absztmp;
                    E[i]        = absztmp;
                    if(absztmp != dzero)
                        ztmp = (PLASMA_Complex64_t) (ztmp / absztmp);
                    else
                        ztmp = zone;
                    if(i<(N-2)) *A((i+1),(i+2)) = *A((i+1),(i+2)) * ztmp;
                    /* for Q: ZSCAL should be done in case of WANTQ. HERE NEED THE multiply by CONJ(T) */
                }
            }
        } /* end PlasmaUpper*/
#else
       if( uplo == PlasmaLower ){
           for (i=0; i < N-1; i++) {
               D[i] = *A(i,   i);
               E[i] = *A(i+1, i);
           }
           D[i] = *A(i, i);
       } else {
           for (i=0; i < N-1; i++) {
               D[i] = *A(i, i  );
               E[i] = *A(i, i+1);
           }
           D[i] = *A(i, i);
       }
#endif
       return;
    }

    /* Case N<NB ==> matrix is very small and better to call lapack XHETRD. */
    if( N <= 0 ) /* this will be removed we don t need it. */
    {
        PLASMA_Complex64_t *work, *TTau;
        int info, ldwork = N*N;
        work = (PLASMA_Complex64_t *) plasma_shared_alloc(plasma, ldwork, PlasmaComplexDouble);
        TTau = (PLASMA_Complex64_t *) plasma_shared_alloc(plasma, N,   PlasmaComplexDouble);

        info = LAPACKE_zhetrd_work(LAPACK_COL_MAJOR, lapack_const(uplo), N,
                                   A(0,0), A.lm, D, E, TTau, work, ldwork);
        plasma_shared_free(plasma, (void*) work);
        plasma_shared_free(plasma, (void*) TTau);

        if( info == 0 )
            sequence->status = PLASMA_SUCCESS;
        else
            plasma_sequence_flush(plasma->quark, sequence, request, info);
        return;
    }

    /* General case NB > 1 && N > NB */
    C   = (PLASMA_Complex64_t *) plasma_shared_alloc(plasma, N,   PlasmaComplexDouble);
    S   = (PLASMA_Complex64_t *) plasma_shared_alloc(plasma, N,   PlasmaComplexDouble);

    /***************************************************************************
     *                       START BULGE CHASING CODE
     **************************************************************************/
    /*
     * Initialisation of local parameter. those parameter should be
     * input or tuned parameter.
     */
    grsiz = 1;
    if( NB > 160 ) {
        grsiz = 1;
    }
    else if( NB > 100 ) {
        grsiz = 1;
        /*
        if( N < 5000 )
            grsiz = 1;
        else
            grsiz = 2;
        */
    } else {
        grsiz = 2;
    }
    grsiz = max(1, grsiz);

    /*grsiz=1;*/
    /*printf("  Version -dp- N %5d   NB %5d    lcNB %5d  grsiz %5d    A.ln %5d   A.nb %5d \n",N,NB,lcNB,grsiz,A.ln,A.nb);*/
    for (blksweep = 0; blksweep<NT; blksweep++){
       lcNB =  blksweep == NT-1 ? A.n-blksweep*A.nb : A.nb;
       /*printf("  Version -dp- N %5d   NB %5d    lcNB %5d  grsiz %5d blksweep%5d     NT %5d \n",N,NB,lcNB,grsiz,blksweep,NT);*/
       for (lcsweep = 0; lcsweep<lcNB; lcsweep++){
          for (blkid = blksweep; blkid<NT; blkid=blkid+grsiz){
              lcgrsiz = (blkid+1) < NT ? grsiz : NT-blkid;
              /*printf("  Version -dp- N %5d   NB %5d    lcNB %5d  grsiz %5d lcgrsiz %5d  blkid %5d \n",N,NB,lcNB,grsiz,lcgrsiz,blkid);*/
              QUARK_CORE_ztrdalg_v2(
                        plasma->quark, &task_flags,
                        uplo, &A, C, S,
                        lcgrsiz, lcsweep, blkid, blksweep);
             }
       }
    }
    /*
     * Barrier used only for now, to be sure that everything
     * is done before copying the D and E and free workspace.
     * this will be removed later when D and E are directly filled
     * during the bulge process.
     */
    QUARK_Barrier(plasma->quark);
    tblg   += Wtimming();
    printf("   done with bulge %lf \n\n\n",tblg);

    plasma_shared_free(plasma, (void*) C);
    plasma_shared_free(plasma, (void*) S);

    /*
     * STORE THE RESULTING diagonal/off-diagonal in D AND E
     */
    memset(D, 0,  N   *sizeof(double));
    memset(E, 0, (N-1)*sizeof(double));
    /* Make diagonal and superdiagonal elements real,
     * storing them in D and E
     */
    /* In complex case, the off diagonal element are
     * not necessary real. we have to make off-diagonal
     * elements real and copy them to E.
     * When using HouseHolder elimination,
     * the ZLARFG give us a real as output so, all the
     * diagonal/off-diagonal element except the last one are already
     * real and thus we need only to take the abs of the last
     * one.
     *  */
#ifdef COMPLEX
    if(uplo==PlasmaLower){
       for (i=0; i < N-1 ; i++)
       {
          D[i] = creal( *A(i,i) );
          /*
           * Alternative for Householder case, all off-diag
           * are real except the last off-diag, where we
           * have to take the abs
           */
          if(i<(N-2))
              E[i] = creal(*A(i+1, i));
          else
              E[i] = cabs( *A(i+1, i));
       }
        D[i] = creal( *A(i, i) );
    } else { /* PlasmaUpper */
        for (i=0; i<N-1; i++)
        {
            D[i]  =  creal( *A(i,i) );
            /*
             * Alternative for Householder case, all off-diag
             * are real except the last off-diag, where we
             * have to take the abs
             */
            if( i < (N-2) )
                E[i] = creal(*A(i, (i+1)));
            else
                E[i] = cabs(*A(i, (i+1)));
        }
        D[i] = creal( *A(i, i) );
    } /* end PlasmaUpper */
#else
    if( uplo == PlasmaLower ){
        for (i=0; i < N-1; i++) {
            D[i] = *A(i,   i);
            E[i] = *A(i+1, i);
        }
        D[i] = *A(i, i);
    } else {
        for (i=0; i < N-1; i++) {
            D[i] = *A(i, i  );
            E[i] = *A(i, i+1);
        }
        D[i] = *A(i, i);
    }
#endif

} /* END FUNCTION */
