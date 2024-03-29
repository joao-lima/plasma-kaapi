      PROGRAM EXAMPLE_DPOSV_F
*     
*********************************************************************
*     PLASMA example routine (version 2.4.2)                        
*     Author: Bilel Hadri                                           
*     Release Date: November, 15th 2010                             
*     PLASMA is a software package provided by Univ. of Tennessee,  
*     Univ. of California Berkeley and Univ. of Colorado Denver.    
*     @generated d Thu Sep 15 12:09:38 2011
*********************************************************************
*     
      IMPLICIT NONE
*     
      INCLUDE "plasmaf.h"
*     
*     Purpose
*     =======
*     
*     FORTRAN EXAMPLE FOR PLASMA_DPOSV
*     Example for solving a system of linear equations using Cholesky 
*     factorization
*     
*     =====================================================================
*     
*     .. Parameters ..
      INTEGER           CORES, N, NRHS
      PARAMETER         ( CORES = 2 )
      PARAMETER         ( N = 15 )
      PARAMETER         ( NRHS = 5 )
      COMPLEX*16        ZONE
      PARAMETER         ( ZONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      COMPLEX*16  A1( N, N ), B1( N, NRHS ), WORK( 2*N )
      COMPLEX*16  A2( N, N ), B2( N, NRHS )
      DOUBLE PRECISION  D( N )
      DOUBLE PRECISION  RWORK ( N )
      INTEGER           I, INFO
      INTEGER           ISEED( 4 )
      DOUBLE PRECISION  XNORM, ANORM, BNORM, RNORM, EPS
      DOUBLE PRECISION  DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL          DLARNV, ZLAGHE, DLAMCH, DLANGE
      EXTERNAL          PLASMA_INIT, PLASMA_DPOSV, PLASMA_FINALIZE
      EXTERNAL          DGEMM
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
      CALL DLARNV( 1, ISEED, N, D )
      CALL ZLAGHE( N, N-1, D, A1, N, ISEED, WORK, INFO )
*     
*     Make it definite positive
*     
      DO I = 1, N
         A1( I, I ) = A1( I, I ) + N
      ENDDO
      A2(:,:)=A1(:,:)
*     
*     Initialization of the RHS
*     
      CALL ZLARNV( 1, ISEED, N*NRHS, B1 )
      B2(:,:)=B1(:,:)
*     
*     Perform the Cholesky solve
*     
      CALL PLASMA_DPOSV( PlasmaUpper, N, NRHS, A2, N, B2, N, INFO )
*     
*     Check the solution
*     
      XNORM = DLANGE('I',N, NRHS, B2, N, RWORK)
      ANORM = DLANGE('I',N, N, A1, N, RWORK)
      BNORM = DLANGE('I',N, NRHS, B1, N, RWORK)

      CALL DGEMM('No transpose','No transpose', N, NRHS, N, ZONE,
     $     A1, N, B2, N, -ZONE, B1, N)

      RNORM = DLANGE('I',N, NRHS, B1, N, RWORK)

      EPS= DLAMCH('Epsilon')

      WRITE(*,*) '============'
      WRITE(*,*) 'Checking the Residual of the solution '
      WRITE(*,*) '-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps)=',
     $     RNORM / ((ANORM * XNORM + BNORM) * N * EPS)

      IF ((RNORM > 60.0).AND.( INFO < 0 )) THEN
         WRITE(*,*) "-- Error in DPOSV example !"
      ELSE
         WRITE(*,*) "-- Run of DPOSV example successful !"
      ENDIF
*     
*     Finalize Plasma
*     
      CALL PLASMA_FINALIZE( INFO )
*     
*     End of EXAMPLE_DPOSV.
*     
      END PROGRAM EXAMPLE_DPOSV_F
