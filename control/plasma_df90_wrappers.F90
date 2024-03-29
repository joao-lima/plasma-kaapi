!
!     Copyright © 2011 The Numerical Algorithms Group Ltd. All rights reserved.
!   
!     Redistribution and use in source and binary forms, with or without
!     modification, are permitted provided that the following conditions are
!     met:
!     - Redistributions of source code must retain the above copyright notice,
!       this list of conditions, and the following disclaimer.
!     - Redistributions in binary form must reproduce the above copyright
!       notice, this list of conditions and the following disclaimer listed in
!       this license in the documentation and/or other materials provided with
!       the distribution.
!     - Neither the name of the copyright holders nor the names of its
!       contributors may be used to endorse or promote products derived from
!       this software without specific prior written permission.
!     
!     This software is provided by the copyright holders and contributors "as
!     is" and any express or implied warranties, including, but not limited
!     to, the implied warranties of merchantability and fitness for a
!     particular purpose are disclaimed. in no event shall the copyright owner
!     or contributors be liable for any direct, indirect, incidental, special,
!     exemplary, or consequential damages (including, but not limited to,
!     procurement of substitute goods or services; loss of use, data, or
!     profits; or business interruption) however caused and on any theory of
!     liability, whether in contract, strict liability, or tort (including
!     negligence or otherwise) arising in any way out of the use of this
!     software, even if advised of the possibility of such damage.
!
!
!
! @file plasma_df90_wrappers.F90
!
!  PLASMA fortran wrapper for BLAS and LAPACK subroutines.
!  PLASMA is a software package provided by Univ. of Tennessee,
!  Univ. of California Berkeley and Univ. of Colorado Denver
!
! @version 2.4.2
! @author Numerical Algorithm Group
! @author Mathieu Faverge
! @date 2011-09-15
! @generated d Thu Sep 15 12:09:27 2011
!
!
! Wrappers to PLASMA functions are provided for the following BLAS
! subroutines since the PLASMA and BLAS interfaces match exactly:
! DGEMM  PLASMA_dgemm
! DSYMM  PLASMA_dsymm
! DSYR2K PLASMA_dsyr2k
! DSYRK  PLASMA_dsyrk
! DSYMM  PLASMA_dsymm
! DSYR2K PLASMA_dsyr2k
! DSYRK  PLASMA_dsyrk
! DTRMM  PLASMA_dtrmm
! DTRSM  PLASMA_dtrsm
!
! Wrappers to PLASMA functions are provided for the following LAPACK
! subroutines since the PLASMA and LAPACK interfaces match exactly:
! DGESV  PLASMA_dgesv
! DGETRF PLASMA_dgetrf
! DGETRS PLASMA_dgetrs
! DSYGST PLASMA_dsygst
! DLASWP PLASMA_dlaswp
! DLAUUM PLASMA_dlauum
! DPOSV  PLASMA_dposv
! DPOTRF PLASMA_dpotrf
! DPOTRI PLASMA_dpotri
! DPOTRS PLASMA_dpotrs
! DTRTRI PLASMA_dtrtri
#define PRECISION_d

      subroutine plasma_wrap_DGESV(N,NRHS,A,LDA,IPIV,B,LDB,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: N
            integer, intent(in) :: NRHS
            integer, intent(out) :: INFO
            integer, intent(out), target :: IPIV(*)
            double precision, intent(inout), target :: A(LDA,*)
            double precision, intent(inout), target :: B(LDB,*)
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_DGESV"
            call PLASMA_DGESV(N,NRHS,A,LDA,IPIV,B,LDB,INFO)
      end subroutine plasma_wrap_DGESV

      subroutine plasma_wrap_DGETRF(M,N,A,LDA,IPIV,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: M
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            integer, intent(out), target :: IPIV(*)
            double precision, intent(inout), target :: A(LDA,*)
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_DGETRF"
            call PLASMA_DGETRF(M,N,A,LDA,IPIV,INFO)
      end subroutine plasma_wrap_DGETRF

      subroutine plasma_wrap_DGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: N
            integer, intent(in) :: NRHS
            integer, intent(in), target :: IPIV(*)
            integer, intent(out) :: INFO
            character, intent(in) :: TRANS
            double precision, intent(in), target :: A(LDA,*)
            double precision, intent(inout), target :: B(LDB,*)
            integer :: local_TRANS
            if(TRANS=='N' .or. TRANS=='n')then
               local_TRANS = PlasmaNoTrans
            else if(TRANS=='T' .or. TRANS=='t')then
               local_TRANS = PlasmaTrans
            else if(TRANS=='C' .or. TRANS=='c')then
               local_TRANS = PlasmaTrans
            else
               local_TRANS = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_DGETRS"
            call PLASMA_DGETRS(local_TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
      end subroutine plasma_wrap_DGETRS

      subroutine plasma_wrap_DSYGST(ITYPE,UPLO,N,A,LDA,B,LDB,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: ITYPE
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            double precision, intent(in), target :: B(LDB,*)
            double precision, intent(inout), target :: A(LDA,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_DSYGST"
            call PLASMA_DSYGST(ITYPE,local_UPLO,N,A,LDA,B,LDB,INFO)
      end subroutine plasma_wrap_DSYGST

      subroutine plasma_wrap_DLASWP(N,A,LDA,K1,K2,IPIV,INCX)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: INCX
            integer, intent(in) :: K1
            integer, intent(in) :: K2
            integer, intent(in) :: LDA
            integer, intent(in) :: N
            integer, intent(in), target :: IPIV(*)
            double precision, intent(inout), target :: A(LDA,*)
            integer :: local_ret
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_DLASWP"
            call PLASMA_DLASWP(N,A,LDA,K1,K2,IPIV,INCX,local_ret)
      end subroutine plasma_wrap_DLASWP

      subroutine plasma_wrap_DLAUUM(UPLO,N,A,LDA,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            double precision, intent(inout), target :: A(LDA,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_DLAUUM"
            call PLASMA_DLAUUM(local_UPLO,N,A,LDA,INFO)
      end subroutine plasma_wrap_DLAUUM

      subroutine plasma_wrap_DPOSV(UPLO,N,NRHS,A,LDA,B,LDB,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: N
            integer, intent(in) :: NRHS
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            double precision, intent(inout), target :: A(LDA,*)
            double precision, intent(inout), target :: B(LDB,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_DPOSV"
            call PLASMA_DPOSV(local_UPLO,N,NRHS,A,LDA,B,LDB,INFO)
      end subroutine plasma_wrap_DPOSV

      subroutine plasma_wrap_DPOTRF(UPLO,N,A,LDA,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            double precision, intent(inout), target :: A(LDA,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_DPOTRF"
            call PLASMA_DPOTRF(local_UPLO,N,A,LDA,INFO)
      end subroutine plasma_wrap_DPOTRF

      subroutine plasma_wrap_DPOTRI(UPLO,N,A,LDA,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            double precision, intent(inout), target :: A(LDA,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_DPOTRI"
            call PLASMA_DPOTRI(local_UPLO,N,A,LDA,INFO)
      end subroutine plasma_wrap_DPOTRI

      subroutine plasma_wrap_DPOTRS(UPLO,N,NRHS,A,LDA,B,LDB,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: N
            integer, intent(in) :: NRHS
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            double precision, intent(in), target :: A(LDA,*)
            double precision, intent(inout), target :: B(LDB,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_DPOTRS"
            call PLASMA_DPOTRS(local_UPLO,N,NRHS,A,LDA,B,LDB,INFO)
      end subroutine plasma_wrap_DPOTRS

      subroutine plasma_wrap_DTRTRI(UPLO,DIAG,N,A,LDA,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            character, intent(in) :: DIAG
            character, intent(in) :: UPLO
            double precision, intent(inout), target :: A(LDA,*)
            integer :: local_DIAG
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if(DIAG=='U' .or. DIAG=='u')then
               local_DIAG = PlasmaUnit
            else if(DIAG=='N' .or. DIAG=='n')then
               local_DIAG = PlasmaNonUnit
            else
               local_DIAG = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_DTRTRI"
            call PLASMA_DTRTRI(local_UPLO,local_DIAG,N,A,LDA,INFO)
      end subroutine plasma_wrap_DTRTRI

      subroutine plasma_wrap_DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: K
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: LDC
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: TRANSA
            character, intent(in) :: TRANSB
            double precision, intent(in) :: ALPHA
            double precision, intent(in) :: BETA
            double precision, intent(in), target :: A(LDA,*)
            double precision, intent(in), target :: B(LDB,*)
            double precision, intent(inout), target :: C(LDC,*)
            integer :: local_TRANSA
            integer :: local_TRANSB
            integer :: local_ret
            if(TRANSA=='N' .or. TRANSA=='n')then
               local_TRANSA = PlasmaNoTrans
            else if(TRANSA=='T' .or. TRANSA=='t')then
               local_TRANSA = PlasmaTrans
            else if(TRANSA=='C' .or. TRANSA=='c')then
               local_TRANSA = PlasmaTrans
            else
               local_TRANSA = -1
            end if
            if(TRANSB=='N' .or. TRANSB=='n')then
               local_TRANSB = PlasmaNoTrans
            else if(TRANSB=='T' .or. TRANSB=='t')then
               local_TRANSB = PlasmaTrans
            else if(TRANSB=='C' .or. TRANSB=='c')then
               local_TRANSB = PlasmaTrans
            else
               local_TRANSB = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_DGEMM"
            call PLASMA_DGEMM(local_TRANSA,local_TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_DGEMM

#if defined(PRECISION_z) || defined(PRECISION_c)
      subroutine plasma_wrap_DSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: LDC
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: SIDE
            character, intent(in) :: UPLO
            double precision, intent(in) :: ALPHA
            double precision, intent(in) :: BETA
            double precision, intent(in), target :: A(LDA,*)
            double precision, intent(in), target :: B(LDB,*)
            double precision, intent(inout), target :: C(LDC,*)
            integer :: local_SIDE
            integer :: local_UPLO
            integer :: local_ret
            if(SIDE=='L' .or. SIDE=='l')then
               local_SIDE = PlasmaLeft
            else if(SIDE=='R' .or. SIDE=='r')then
               local_SIDE = PlasmaRight
            else
               local_SIDE = -1
            end if
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_DSYMM"
            call PLASMA_DSYMM(local_SIDE,local_UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_DSYMM

      subroutine plasma_wrap_DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: K
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: LDC
            integer, intent(in) :: N
            character, intent(in) :: TRANS
            character, intent(in) :: UPLO
            double precision, intent(in) :: ALPHA
            double precision, intent(in), target :: A(LDA,*)
            double precision, intent(in), target :: B(LDB,*)
            double precision, intent(inout), target :: C(LDC,*)
            double precision, intent(in) :: BETA
            integer :: local_TRANS
            integer :: local_UPLO
            integer :: local_ret
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if(TRANS=='N' .or. TRANS=='n')then
               local_TRANS = PlasmaNoTrans
            else if(TRANS=='T' .or. TRANS=='t')then
               local_TRANS = PlasmaTrans
            else if(TRANS=='C' .or. TRANS=='c')then
               local_TRANS = PlasmaTrans
            else
               local_TRANS = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_DSYR2K"
            call PLASMA_DSYR2K(local_UPLO,local_TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_DSYR2K

      subroutine plasma_wrap_DSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: K
            integer, intent(in) :: LDA
            integer, intent(in) :: LDC
            integer, intent(in) :: N
            character, intent(in) :: TRANS
            character, intent(in) :: UPLO
            double precision, intent(in), target :: A(LDA,*)
            double precision, intent(inout), target :: C(LDC,*)
            double precision, intent(in) :: ALPHA
            double precision, intent(in) :: BETA
            integer :: local_TRANS
            integer :: local_UPLO
            integer :: local_ret
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if(TRANS=='N' .or. TRANS=='n')then
               local_TRANS = PlasmaNoTrans
            else if(TRANS=='T' .or. TRANS=='t')then
               local_TRANS = PlasmaTrans
            else if(TRANS=='C' .or. TRANS=='c')then
               local_TRANS = PlasmaTrans
            else
               local_TRANS = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_DSYRK"
            call PLASMA_DSYRK(local_UPLO,local_TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_DSYRK
#endif

      subroutine plasma_wrap_DSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: LDC
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: SIDE
            character, intent(in) :: UPLO
            double precision, intent(in) :: ALPHA
            double precision, intent(in) :: BETA
            double precision, intent(in), target :: A(LDA,*)
            double precision, intent(in), target :: B(LDB,*)
            double precision, intent(inout), target :: C(LDC,*)
            integer :: local_SIDE
            integer :: local_UPLO
            integer :: local_ret
            if(SIDE=='L' .or. SIDE=='l')then
               local_SIDE = PlasmaLeft
            else if(SIDE=='R' .or. SIDE=='r')then
               local_SIDE = PlasmaRight
            else
               local_SIDE = -1
            end if
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_DSYMM"
            call PLASMA_DSYMM(local_SIDE,local_UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_DSYMM

      subroutine plasma_wrap_DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: K
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: LDC
            integer, intent(in) :: N
            character, intent(in) :: TRANS
            character, intent(in) :: UPLO
            double precision, intent(in) :: ALPHA
            double precision, intent(in) :: BETA
            double precision, intent(in), target :: A(LDA,*)
            double precision, intent(in), target :: B(LDB,*)
            double precision, intent(inout), target :: C(LDC,*)
            integer :: local_TRANS
            integer :: local_UPLO
            integer :: local_ret
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if(TRANS=='N' .or. TRANS=='n')then
               local_TRANS = PlasmaNoTrans
            else if(TRANS=='T' .or. TRANS=='t')then
               local_TRANS = PlasmaTrans
            else if(TRANS=='C' .or. TRANS=='c')then
               local_TRANS = PlasmaTrans
            else
               local_TRANS = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_DSYR2K"
            call PLASMA_DSYR2K(local_UPLO,local_TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_DSYR2K

      subroutine plasma_wrap_DSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: K
            integer, intent(in) :: LDA
            integer, intent(in) :: LDC
            integer, intent(in) :: N
            character, intent(in) :: TRANS
            character, intent(in) :: UPLO
            double precision, intent(in) :: ALPHA
            double precision, intent(in) :: BETA
            double precision, intent(in), target :: A(LDA,*)
            double precision, intent(inout), target :: C(LDC,*)
            integer :: local_TRANS
            integer :: local_UPLO
            integer :: local_ret
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if(TRANS=='N' .or. TRANS=='n')then
               local_TRANS = PlasmaNoTrans
            else if(TRANS=='T' .or. TRANS=='t')then
               local_TRANS = PlasmaTrans
            else if(TRANS=='C' .or. TRANS=='c')then
               local_TRANS = PlasmaTrans
            else
               local_TRANS = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_DSYRK"
            call PLASMA_DSYRK(local_UPLO,local_TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_DSYRK

      subroutine plasma_wrap_DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: DIAG
            character, intent(in) :: SIDE
            character, intent(in) :: TRANSA
            character, intent(in) :: UPLO
            double precision, intent(in) :: ALPHA
            double precision, intent(in), target :: A(LDA,*)
            double precision, intent(inout), target :: B(LDB,*)
            integer :: local_DIAG
            integer :: local_SIDE
            integer :: local_TRANSA
            integer :: local_UPLO
            integer :: local_ret
            if(SIDE=='L' .or. SIDE=='l')then
               local_SIDE = PlasmaLeft
            else if(SIDE=='R' .or. SIDE=='r')then
               local_SIDE = PlasmaRight
            else
               local_SIDE = -1
            end if
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if(TRANSA=='N' .or. TRANSA=='n')then
               local_TRANSA = PlasmaNoTrans
            else if(TRANSA=='T' .or. TRANSA=='t')then
               local_TRANSA = PlasmaTrans
            else if(TRANSA=='C' .or. TRANSA=='c')then
               local_TRANSA = PlasmaTrans
            else
               local_TRANSA = -1
            end if
            if(DIAG=='U' .or. DIAG=='u')then
               local_DIAG = PlasmaUnit
            else if(DIAG=='N' .or. DIAG=='n')then
               local_DIAG = PlasmaNonUnit
            else
               local_DIAG = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_DTRMM"
            call PLASMA_DTRMM(local_SIDE,local_UPLO,local_TRANSA,local_DIAG,M,N,ALPHA,A,LDA,B,LDB,local_ret)
      end subroutine plasma_wrap_DTRMM

      subroutine plasma_wrap_DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: DIAG
            character, intent(in) :: SIDE
            character, intent(in) :: TRANSA
            character, intent(in) :: UPLO
            double precision, intent(in) :: ALPHA
            double precision, intent(in), target :: A(LDA,*)
            double precision, intent(inout), target :: B(LDB,*)
            integer :: local_DIAG
            integer :: local_SIDE
            integer :: local_TRANSA
            integer :: local_UPLO
            integer :: local_ret
            if(SIDE=='L' .or. SIDE=='l')then
               local_SIDE = PlasmaLeft
            else if(SIDE=='R' .or. SIDE=='r')then
               local_SIDE = PlasmaRight
            else
               local_SIDE = -1
            end if
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if(TRANSA=='N' .or. TRANSA=='n')then
               local_TRANSA = PlasmaNoTrans
            else if(TRANSA=='T' .or. TRANSA=='t')then
               local_TRANSA = PlasmaTrans
            else if(TRANSA=='C' .or. TRANSA=='c')then
               local_TRANSA = PlasmaTrans
            else
               local_TRANSA = -1
            end if
            if(DIAG=='U' .or. DIAG=='u')then
               local_DIAG = PlasmaUnit
            else if(DIAG=='N' .or. DIAG=='n')then
               local_DIAG = PlasmaNonUnit
            else
               local_DIAG = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_DTRSM"
            call PLASMA_DTRSM(local_SIDE,local_UPLO,local_TRANSA,local_DIAG,M,N,ALPHA,A,LDA,B,LDB,local_ret)
      end subroutine plasma_wrap_DTRSM
