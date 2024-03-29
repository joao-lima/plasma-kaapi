        PROGRAM EXAMPLE_SGESV_F
*
*********************************************************************
*     PLASMA example routine (version 2.4.2)                        
*     Author: Bilel Hadri                                           
*     Release Date: November, 15th 2010                             
*     PLASMA is a software package provided by Univ. of Tennessee,  
*     Univ. of California Berkeley and Univ. of Colorado Denver.    
*     @generated s Thu Sep 15 12:09:38 2011
*********************************************************************
*
        IMPLICIT NONE
*
        INCLUDE "plasmaf.h"
*
*  Purpose
*  =======
*
*   FORTRAN EXAMPLE FOR PLASMA_SGESV
*   Example for solving a system of linear equations using LU 
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER           CORES, N, NRHS
      PARAMETER         ( CORES = 2 )
      PARAMETER         ( N = 10 )
      PARAMETER         ( NRHS = 5 )
      COMPLEX*16        ZONE
      PARAMETER         ( ZONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      COMPLEX*16        A1( N, N ), B1( N, NRHS )
      COMPLEX*16        A2( N, N ), B2( N, NRHS )
      DOUBLE PRECISION  RWORK( N )
      INTEGER*4         HL( 2 ), HPIV( 2 )
      INTEGER           I, INFO
      INTEGER           ISEED( 4 )
      DOUBLE PRECISION  XNORM, ANORM, BNORM, RNORM, EPS
      DOUBLE PRECISION  DLAMCH, SLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL          ZLARNV, DLAMCH, SLANGE
      EXTERNAL          PLASMA_INIT, PLASMA_ALLOC_WORKSPACE_SGESV_INCPIV
      EXTERNAL          PLASMA_SGESV, PLASMA_FINALIZE
      EXTERNAL          PLASMA_DEALLOC_HANDLE
      EXTERNAL          SGEMM
*     ..
*     .. Executable Statements ..
*
      DO  I = 1, 4
          ISEED( I ) = 1
      ENDDO
*
*     Initialize Plasma
*
      CALL PLASMA_INIT( CORES, INFO )
      WRITE(*,*) "-- PLASMA is initialized on", CORES, "cores."
*
*     Initialization of the matrix A1
*
      CALL ZLARNV( 1, ISEED, N*N, A1 )
      A2(:,:)=A1(:,:)
*
*     Initialization of the RHS
*
      CALL ZLARNV( 1, ISEED, N*NRHS, B1 )
      B2(:,:)=B1(:,:)
*
*     Allocate L and IPIV
*
      CALL PLASMA_ALLOC_WORKSPACE_SGESV_INCPIV( N, HL, HPIV, INFO )
*
*     PLASMA SGESV
*
      CALL PLASMA_SGESV_INCPIV( N, NRHS, A2, N, HL, HPIV,B2, N, INFO )
*
*     Check the solution
*
      XNORM = SLANGE('I',N, NRHS, B2, N, RWORK)
      ANORM = SLANGE('I',N, N, A1, N, RWORK)
      BNORM = SLANGE('I',N, NRHS, B1, N, RWORK)

      CALL SGEMM('No transpose','No transpose', N, NRHS, N, ZONE,
     $              A1, N, B2, N, -ZONE, B1, N)

      RNORM = SLANGE('I',N, NRHS, B1, N, RWORK)

      EPS= DLAMCH('Epsilon')

      WRITE(*,*) '============'
      WRITE(*,*) 'Checking the Residual of the solution '
      WRITE(*,*) '-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps)=',
     $        RNORM / ((ANORM * XNORM + BNORM) * N * EPS)

      IF ((RNORM > 60.0).AND.( INFO < 0 )) THEN
          WRITE(*,*) "-- Error in SGESV example !"
      ELSE
          WRITE(*,*) "-- Run of SGESV example successful !"
      ENDIF

*
*     Deallocate L and IPIV
*
      CALL PLASMA_DEALLOC_HANDLE( HL, INFO )
      CALL PLASMA_DEALLOC_HANDLE( HPIV, INFO )
*
*     Finalize Plasma
*
      CALL PLASMA_FINALIZE( INFO )
*
*     End of EXAMPLE_SGESV.
*
      END PROGRAM EXAMPLE_SGESV_F
