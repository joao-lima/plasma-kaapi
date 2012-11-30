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
! @file plasma_df90.F90
!
!  PLASMA fortran 90 interface
!  PLASMA is a software package provided by Univ. of Tennessee,
!  Univ. of California Berkeley and Univ. of Colorado Denver
!
! @version 2.4.2
! @author Numerical Algorithm Group
! @author Mathieu Faverge
! @date 2011-09-15
! @generated d Thu Sep 15 12:09:27 2011
!
#define PRECISION_d

module plasma_d

!      private PLASMA_dgebrd_c
      interface
         function PLASMA_dgebrd_c(jobu,jobvt,M,N,A,LDA,D,E,U,LDU,VT,LDVT,T) &
          & bind(c, name='PLASMA_dgebrd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgebrd_c
            integer(kind=c_int), value :: jobu
            integer(kind=c_int), value :: jobvt
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: D
            type(c_ptr), value :: E
            type(c_ptr), value :: U
            integer(kind=c_int), value :: LDU
            type(c_ptr), value :: VT
            integer(kind=c_int), value :: LDVT
            type(c_ptr), value :: T
         end function PLASMA_dgebrd_c
      end interface

!      private PLASMA_dgelqf_c
      interface
         function PLASMA_dgelqf_c(M,N,A,LDA,T) &
          & bind(c, name='PLASMA_dgelqf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgelqf_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: T
         end function PLASMA_dgelqf_c
      end interface

!      private PLASMA_dgelqs_c
      interface
         function PLASMA_dgelqs_c(M,N,NRHS,A,LDA,T,B,LDB) &
          & bind(c, name='PLASMA_dgelqs')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgelqs_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_dgelqs_c
      end interface

!      private PLASMA_dgels_c
      interface
         function PLASMA_dgels_c(trans,M,N,NRHS,A,LDA,T,B,LDB) &
          & bind(c, name='PLASMA_dgels')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgels_c
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_dgels_c
      end interface

!      private PLASMA_dgemm_c
      interface
         function PLASMA_dgemm_c(transA,transB,M,N,K,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='PLASMA_dgemm')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgemm_c
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: transB
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            real(kind=c_double), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
         end function PLASMA_dgemm_c
      end interface

!      private PLASMA_dgeqrf_c
      interface
         function PLASMA_dgeqrf_c(M,N,A,LDA,T) &
          & bind(c, name='PLASMA_dgeqrf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgeqrf_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: T
         end function PLASMA_dgeqrf_c
      end interface

!      private PLASMA_dgeqrs_c
      interface
         function PLASMA_dgeqrs_c(M,N,NRHS,A,LDA,T,B,LDB) &
          & bind(c, name='PLASMA_dgeqrs')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgeqrs_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_dgeqrs_c
      end interface

!      private PLASMA_dgesv_c
      interface
         function PLASMA_dgesv_c(N,NRHS,A,LDA,IPIV,B,LDB) &
          & bind(c, name='PLASMA_dgesv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgesv_c
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_dgesv_c
      end interface

!      private PLASMA_dgesv_incpiv_c
      interface
         function PLASMA_dgesv_incpiv_c(N,NRHS,A,LDA,L,IPIV,B,LDB) &
          & bind(c, name='PLASMA_dgesv_incpiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgesv_incpiv_c
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_dgesv_incpiv_c
      end interface

!      private PLASMA_dgesvd_c
      interface
         function PLASMA_dgesvd_c(jobu,jobvt,M,N,A,LDA,S,U,LDU,VT,LDVT,T) &
          & bind(c, name='PLASMA_dgesvd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgesvd_c
            integer(kind=c_int), value :: jobu
            integer(kind=c_int), value :: jobvt
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: S
            type(c_ptr), value :: U
            integer(kind=c_int), value :: LDU
            type(c_ptr), value :: VT
            integer(kind=c_int), value :: LDVT
            type(c_ptr), value :: T
         end function PLASMA_dgesvd_c
      end interface

!      private PLASMA_dgetrf_c
      interface
         function PLASMA_dgetrf_c(M,N,A,LDA,IPIV) &
          & bind(c, name='PLASMA_dgetrf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgetrf_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: IPIV
         end function PLASMA_dgetrf_c
      end interface

!      private PLASMA_dgetrf_incpiv_c
      interface
         function PLASMA_dgetrf_incpiv_c(M,N,A,LDA,L,IPIV) &
          & bind(c, name='PLASMA_dgetrf_incpiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgetrf_incpiv_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
         end function PLASMA_dgetrf_incpiv_c
      end interface

!      private PLASMA_dgetrs_c
      interface
         function PLASMA_dgetrs_c(trans,N,NRHS,A,LDA,IPIV,B,LDB) &
          & bind(c, name='PLASMA_dgetrs')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgetrs_c
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_dgetrs_c
      end interface

!      private PLASMA_dgetrs_incpiv_c
      interface
         function PLASMA_dgetrs_incpiv_c(trans,N,NRHS,A,LDA,L,IPIV,B,LDB) &
          & bind(c, name='PLASMA_dgetrs_incpiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgetrs_incpiv_c
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_dgetrs_incpiv_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
!      private PLASMA_dsymm_c
      interface
         function PLASMA_dsymm_c(side,uplo,M,N,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='PLASMA_dsymm')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsymm_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            real(kind=c_double), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
         end function PLASMA_dsymm_c
      end interface

!      private PLASMA_dsyrk_c
      interface
         function PLASMA_dsyrk_c(uplo,trans,N,K,alpha,A,LDA,beta,C,LDC) &
          & bind(c, name='PLASMA_dsyrk')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsyrk_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            real(kind=c_double), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
         end function PLASMA_dsyrk_c
      end interface

!      private PLASMA_dsyr2k_c
      interface
         function PLASMA_dsyr2k_c(uplo,trans,N,K,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='PLASMA_dsyr2k')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsyr2k_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            real(kind=c_double), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
         end function PLASMA_dsyr2k_c
      end interface
#endif

!      private PLASMA_dsyev_c
      interface
         function PLASMA_dsyev_c(jobz,uplo,N,A,LDA,W,T,Q,LDQ) &
          & bind(c, name='PLASMA_dsyev')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsyev_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
         end function PLASMA_dsyev_c
      end interface

!      private PLASMA_dsygv_c
      interface
         function PLASMA_dsygv_c(itype,jobz,uplo,N,A,LDA,B,LDB,W,T,Q,LDQ) &
          & bind(c, name='PLASMA_dsygv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsygv_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
         end function PLASMA_dsygv_c
      end interface

!      private PLASMA_dsygst_c
      interface
         function PLASMA_dsygst_c(itype,uplo,N,A,LDA,B,LDB) &
          & bind(c, name='PLASMA_dsygst')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsygst_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_dsygst_c
      end interface

!      private PLASMA_dsytrd_c
      interface
         function PLASMA_dsytrd_c(jobz,uplo,N,A,LDA,D,E,T,Q,LDQ) &
          & bind(c, name='PLASMA_dsytrd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsytrd_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: D
            type(c_ptr), value :: E
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
         end function PLASMA_dsytrd_c
      end interface

!      private PLASMA_dlange_c
      interface
         function PLASMA_dlange_c(norm,M,N,A,LDA,work) &
          & bind(c, name='PLASMA_dlange')
            use iso_c_binding
            implicit none
            real(kind=c_double) :: PLASMA_dlange_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: work
         end function PLASMA_dlange_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
!      private PLASMA_dlansy_c
      interface
         function PLASMA_dlansy_c(norm,uplo,N,A,LDA,work) &
          & bind(c, name='PLASMA_dlansy')
            use iso_c_binding
            implicit none
            real(kind=c_double) :: PLASMA_dlansy_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: work
         end function PLASMA_dlansy_c
      end interface
#endif

!      private PLASMA_dlansy_c
      interface
         function PLASMA_dlansy_c(norm,uplo,N,A,LDA,work) &
          & bind(c, name='PLASMA_dlansy')
            use iso_c_binding
            implicit none
            real(kind=c_double) :: PLASMA_dlansy_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: work
         end function PLASMA_dlansy_c
      end interface

!      private PLASMA_dlaswp_c
      interface
         function PLASMA_dlaswp_c(N,A,LDA,K1,K2,IPIV,INCX) &
          & bind(c, name='PLASMA_dlaswp')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dlaswp_c
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            integer(kind=c_int), value :: K1
            integer(kind=c_int), value :: K2
            type(c_ptr), value :: IPIV
            integer(kind=c_int), value :: INCX
         end function PLASMA_dlaswp_c
      end interface

!      private PLASMA_dlauum_c
      interface
         function PLASMA_dlauum_c(uplo,N,A,LDA) &
          & bind(c, name='PLASMA_dlauum')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dlauum_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
         end function PLASMA_dlauum_c
      end interface

! PLASMA_dplgsy uses unsigned ints, not generating Fortran interface

! PLASMA_dplgsy uses unsigned ints, not generating Fortran interface

! PLASMA_dplrnt uses unsigned ints, not generating Fortran interface

!      private PLASMA_dposv_c
      interface
         function PLASMA_dposv_c(uplo,N,NRHS,A,LDA,B,LDB) &
          & bind(c, name='PLASMA_dposv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dposv_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_dposv_c
      end interface

!      private PLASMA_dpotrf_c
      interface
         function PLASMA_dpotrf_c(uplo,N,A,LDA) &
          & bind(c, name='PLASMA_dpotrf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dpotrf_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
         end function PLASMA_dpotrf_c
      end interface

!      private PLASMA_dpotri_c
      interface
         function PLASMA_dpotri_c(uplo,N,A,LDA) &
          & bind(c, name='PLASMA_dpotri')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dpotri_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
         end function PLASMA_dpotri_c
      end interface

!      private PLASMA_dpotrs_c
      interface
         function PLASMA_dpotrs_c(uplo,N,NRHS,A,LDA,B,LDB) &
          & bind(c, name='PLASMA_dpotrs')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dpotrs_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_dpotrs_c
      end interface

!      private PLASMA_dsymm_c
      interface
         function PLASMA_dsymm_c(side,uplo,M,N,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='PLASMA_dsymm')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsymm_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            real(kind=c_double), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
         end function PLASMA_dsymm_c
      end interface

!      private PLASMA_dsyrk_c
      interface
         function PLASMA_dsyrk_c(uplo,trans,N,K,alpha,A,LDA,beta,C,LDC) &
          & bind(c, name='PLASMA_dsyrk')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsyrk_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            real(kind=c_double), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
         end function PLASMA_dsyrk_c
      end interface

!      private PLASMA_dsyr2k_c
      interface
         function PLASMA_dsyr2k_c(uplo,trans,N,K,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='PLASMA_dsyr2k')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsyr2k_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            real(kind=c_double), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
         end function PLASMA_dsyr2k_c
      end interface

!      private PLASMA_dtrmm_c
      interface
         function PLASMA_dtrmm_c(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB) &
          & bind(c, name='PLASMA_dtrmm')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dtrmm_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            real(kind=c_double), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_dtrmm_c
      end interface

!      private PLASMA_dtrsm_c
      interface
         function PLASMA_dtrsm_c(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB) &
          & bind(c, name='PLASMA_dtrsm')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dtrsm_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            real(kind=c_double), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_dtrsm_c
      end interface

!      private PLASMA_dtrsmpl_c
      interface
         function PLASMA_dtrsmpl_c(N,NRHS,A,LDA,L,IPIV,B,LDB) &
          & bind(c, name='PLASMA_dtrsmpl')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dtrsmpl_c
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_dtrsmpl_c
      end interface

!      private PLASMA_dtrtri_c
      interface
         function PLASMA_dtrtri_c(uplo,diag,N,A,LDA) &
          & bind(c, name='PLASMA_dtrtri')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dtrtri_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: diag
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
         end function PLASMA_dtrtri_c
      end interface

!      private PLASMA_dorglq_c
      interface
         function PLASMA_dorglq_c(M,N,K,A,LDA,T,B,LDB) &
          & bind(c, name='PLASMA_dorglq')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dorglq_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_dorglq_c
      end interface

!      private PLASMA_dorgqr_c
      interface
         function PLASMA_dorgqr_c(M,N,K,A,LDA,T,B,LDB) &
          & bind(c, name='PLASMA_dorgqr')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dorgqr_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_dorgqr_c
      end interface

!      private PLASMA_dormlq_c
      interface
         function PLASMA_dormlq_c(side,trans,M,N,K,A,LDA,T,B,LDB) &
          & bind(c, name='PLASMA_dormlq')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dormlq_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_dormlq_c
      end interface

!      private PLASMA_dormqr_c
      interface
         function PLASMA_dormqr_c(side,trans,M,N,K,A,LDA,T,B,LDB) &
          & bind(c, name='PLASMA_dormqr')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dormqr_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_dormqr_c
      end interface

!      private PLASMA_dgecfi_c
      interface
         function PLASMA_dgecfi_c(m,n,A,fin,imb,inb,fout,omb,onb) &
          & bind(c, name='PLASMA_dgecfi')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgecfi_c
            integer(kind=c_int), value :: m
            integer(kind=c_int), value :: n
            type(c_ptr), value :: A
            integer(kind=c_int), value :: fin
            integer(kind=c_int), value :: imb
            integer(kind=c_int), value :: inb
            integer(kind=c_int), value :: fout
            integer(kind=c_int), value :: omb
            integer(kind=c_int), value :: onb
         end function PLASMA_dgecfi_c
      end interface

!      private PLASMA_dgetmi_c
      interface
         function PLASMA_dgetmi_c(m,n,A,fin,mb,nb) &
          & bind(c, name='PLASMA_dgetmi')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgetmi_c
            integer(kind=c_int), value :: m
            integer(kind=c_int), value :: n
            type(c_ptr), value :: A
            integer(kind=c_int), value :: fin
            integer(kind=c_int), value :: mb
            integer(kind=c_int), value :: nb
         end function PLASMA_dgetmi_c
      end interface

!      private PLASMA_dgebrd_Tile_c
      interface
         function PLASMA_dgebrd_Tile_c(jobu,jobvt,A,D,E,U,VT,T) &
          & bind(c, name='PLASMA_dgebrd_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgebrd_Tile_c
            integer(kind=c_int), value :: jobu
            integer(kind=c_int), value :: jobvt
            type(c_ptr), value :: A
            type(c_ptr), value :: D
            type(c_ptr), value :: E
            type(c_ptr), value :: U
            type(c_ptr), value :: VT
            type(c_ptr), value :: T
         end function PLASMA_dgebrd_Tile_c
      end interface

!      private PLASMA_dgelqf_Tile_c
      interface
         function PLASMA_dgelqf_Tile_c(A,T) &
          & bind(c, name='PLASMA_dgelqf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgelqf_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
         end function PLASMA_dgelqf_Tile_c
      end interface

!      private PLASMA_dgelqs_Tile_c
      interface
         function PLASMA_dgelqs_Tile_c(A,B,T) &
          & bind(c, name='PLASMA_dgelqs_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgelqs_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: T
         end function PLASMA_dgelqs_Tile_c
      end interface

!      private PLASMA_dgels_Tile_c
      interface
         function PLASMA_dgels_Tile_c(trans,A,T,B) &
          & bind(c, name='PLASMA_dgels_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgels_Tile_c
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: T
         end function PLASMA_dgels_Tile_c
      end interface

!      private PLASMA_dgemm_Tile_c
      interface
         function PLASMA_dgemm_Tile_c(transA,transB,alpha,A,B,beta,C) &
          & bind(c, name='PLASMA_dgemm_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgemm_Tile_c
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: transB
            real(kind=c_double), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
         end function PLASMA_dgemm_Tile_c
      end interface

!      private PLASMA_dgeqrf_Tile_c
      interface
         function PLASMA_dgeqrf_Tile_c(A,T) &
          & bind(c, name='PLASMA_dgeqrf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgeqrf_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
         end function PLASMA_dgeqrf_Tile_c
      end interface

!      private PLASMA_dgeqrs_Tile_c
      interface
         function PLASMA_dgeqrs_Tile_c(A,T,B) &
          & bind(c, name='PLASMA_dgeqrs_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgeqrs_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: T
         end function PLASMA_dgeqrs_Tile_c
      end interface

!      private PLASMA_dgesv_Tile_c
      interface
         function PLASMA_dgesv_Tile_c(A,IPIV,B) &
          & bind(c, name='PLASMA_dgesv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgesv_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
         end function PLASMA_dgesv_Tile_c
      end interface

!      private PLASMA_dgesv_incpiv_Tile_c
      interface
         function PLASMA_dgesv_incpiv_Tile_c(A,L,IPIV,B) &
          & bind(c, name='PLASMA_dgesv_incpiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgesv_incpiv_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
         end function PLASMA_dgesv_incpiv_Tile_c
      end interface

!      private PLASMA_dgesvd_Tile_c
      interface
         function PLASMA_dgesvd_Tile_c(jobu,jobvt,A,S,U,VT,T) &
          & bind(c, name='PLASMA_dgesvd_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgesvd_Tile_c
            integer(kind=c_int), value :: jobu
            integer(kind=c_int), value :: jobvt
            type(c_ptr), value :: A
            type(c_ptr), value :: S
            type(c_ptr), value :: U
            type(c_ptr), value :: VT
            type(c_ptr), value :: T
         end function PLASMA_dgesvd_Tile_c
      end interface

!      private PLASMA_dgetrf_Tile_c
      interface
         function PLASMA_dgetrf_Tile_c(A,IPIV) &
          & bind(c, name='PLASMA_dgetrf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgetrf_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
         end function PLASMA_dgetrf_Tile_c
      end interface

!      private PLASMA_dgetrf_incpiv_Tile_c
      interface
         function PLASMA_dgetrf_incpiv_Tile_c(A,L,IPIV) &
          & bind(c, name='PLASMA_dgetrf_incpiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgetrf_incpiv_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
         end function PLASMA_dgetrf_incpiv_Tile_c
      end interface

!      private PLASMA_dgetrs_Tile_c
      interface
         function PLASMA_dgetrs_Tile_c(trans,A,IPIV,B) &
          & bind(c, name='PLASMA_dgetrs_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgetrs_Tile_c
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
         end function PLASMA_dgetrs_Tile_c
      end interface

!      private PLASMA_dgetrs_incpiv_Tile_c
      interface
         function PLASMA_dgetrs_incpiv_Tile_c(A,L,IPIV,B) &
          & bind(c, name='PLASMA_dgetrs_incpiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dgetrs_incpiv_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
         end function PLASMA_dgetrs_incpiv_Tile_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
!      private PLASMA_dsymm_Tile_c
      interface
         function PLASMA_dsymm_Tile_c(side,uplo,alpha,A,B,beta,C) &
          & bind(c, name='PLASMA_dsymm_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsymm_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            real(kind=c_double), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
         end function PLASMA_dsymm_Tile_c
      end interface

!      private PLASMA_dsyrk_Tile_c
      interface
         function PLASMA_dsyrk_Tile_c(uplo,trans,alpha,A,beta,C) &
          & bind(c, name='PLASMA_dsyrk_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsyrk_Tile_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            real(kind=c_double), value :: alpha
            type(c_ptr), value :: A
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
         end function PLASMA_dsyrk_Tile_c
      end interface

!      private PLASMA_dsyr2k_Tile_c
      interface
         function PLASMA_dsyr2k_Tile_c(uplo,trans,alpha,A,B,beta,C) &
          & bind(c, name='PLASMA_dsyr2k_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsyr2k_Tile_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            real(kind=c_double), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
         end function PLASMA_dsyr2k_Tile_c
      end interface
#endif

!      private PLASMA_dsyev_Tile_c
      interface
         function PLASMA_dsyev_Tile_c(jobz,uplo,A,W,T,Q) &
          & bind(c, name='PLASMA_dsyev_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsyev_Tile_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
         end function PLASMA_dsyev_Tile_c
      end interface

!      private PLASMA_dsygv_Tile_c
      interface
         function PLASMA_dsygv_Tile_c(itype,jobz,uplo,A,B,W,T,Q) &
          & bind(c, name='PLASMA_dsygv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsygv_Tile_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
         end function PLASMA_dsygv_Tile_c
      end interface

!      private PLASMA_dsygst_Tile_c
      interface
         function PLASMA_dsygst_Tile_c(itype,uplo,A,B) &
          & bind(c, name='PLASMA_dsygst_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsygst_Tile_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
         end function PLASMA_dsygst_Tile_c
      end interface

!      private PLASMA_dsytrd_Tile_c
      interface
         function PLASMA_dsytrd_Tile_c(jobz,uplo,A,D,E,T,Q) &
          & bind(c, name='PLASMA_dsytrd_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsytrd_Tile_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: D
            type(c_ptr), value :: E
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
         end function PLASMA_dsytrd_Tile_c
      end interface

!      private PLASMA_dlange_Tile_c
      interface
         function PLASMA_dlange_Tile_c(norm,A,work) &
          & bind(c, name='PLASMA_dlange_Tile')
            use iso_c_binding
            implicit none
            real(kind=c_double) :: PLASMA_dlange_Tile_c
            integer(kind=c_int), value :: norm
            type(c_ptr), value :: A
            type(c_ptr), value :: work
         end function PLASMA_dlange_Tile_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
!      private PLASMA_dlansy_Tile_c
      interface
         function PLASMA_dlansy_Tile_c(norm,uplo,A,work) &
          & bind(c, name='PLASMA_dlansy_Tile')
            use iso_c_binding
            implicit none
            real(kind=c_double) :: PLASMA_dlansy_Tile_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: work
         end function PLASMA_dlansy_Tile_c
      end interface
#endif

!      private PLASMA_dlansy_Tile_c
      interface
         function PLASMA_dlansy_Tile_c(norm,uplo,A,work) &
          & bind(c, name='PLASMA_dlansy_Tile')
            use iso_c_binding
            implicit none
            real(kind=c_double) :: PLASMA_dlansy_Tile_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: work
         end function PLASMA_dlansy_Tile_c
      end interface

!      private PLASMA_dlaswp_Tile_c
      interface
         function PLASMA_dlaswp_Tile_c(A,K1,K2,IPIV,INCX) &
          & bind(c, name='PLASMA_dlaswp_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dlaswp_Tile_c
            type(c_ptr), value :: A
            integer(kind=c_int), value :: K1
            integer(kind=c_int), value :: K2
            type(c_ptr), value :: IPIV
            integer(kind=c_int), value :: INCX
         end function PLASMA_dlaswp_Tile_c
      end interface

!      private PLASMA_dlauum_Tile_c
      interface
         function PLASMA_dlauum_Tile_c(uplo,A) &
          & bind(c, name='PLASMA_dlauum_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dlauum_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
         end function PLASMA_dlauum_Tile_c
      end interface

! PLASMA_dplgsy_Tile uses unsigned ints, not generating Fortran interface

! PLASMA_dplgsy_Tile uses unsigned ints, not generating Fortran interface

! PLASMA_dplrnt_Tile uses unsigned ints, not generating Fortran interface

!      private PLASMA_dposv_Tile_c
      interface
         function PLASMA_dposv_Tile_c(uplo,A,B) &
          & bind(c, name='PLASMA_dposv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dposv_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
         end function PLASMA_dposv_Tile_c
      end interface

!      private PLASMA_dpotrf_Tile_c
      interface
         function PLASMA_dpotrf_Tile_c(uplo,A) &
          & bind(c, name='PLASMA_dpotrf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dpotrf_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
         end function PLASMA_dpotrf_Tile_c
      end interface

!      private PLASMA_dpotri_Tile_c
      interface
         function PLASMA_dpotri_Tile_c(uplo,A) &
          & bind(c, name='PLASMA_dpotri_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dpotri_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
         end function PLASMA_dpotri_Tile_c
      end interface

!      private PLASMA_dpotrs_Tile_c
      interface
         function PLASMA_dpotrs_Tile_c(uplo,A,B) &
          & bind(c, name='PLASMA_dpotrs_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dpotrs_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
         end function PLASMA_dpotrs_Tile_c
      end interface

!      private PLASMA_dsymm_Tile_c
      interface
         function PLASMA_dsymm_Tile_c(side,uplo,alpha,A,B,beta,C) &
          & bind(c, name='PLASMA_dsymm_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsymm_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            real(kind=c_double), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
         end function PLASMA_dsymm_Tile_c
      end interface

!      private PLASMA_dsyrk_Tile_c
      interface
         function PLASMA_dsyrk_Tile_c(uplo,trans,alpha,A,beta,C) &
          & bind(c, name='PLASMA_dsyrk_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsyrk_Tile_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            real(kind=c_double), value :: alpha
            type(c_ptr), value :: A
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
         end function PLASMA_dsyrk_Tile_c
      end interface

!      private PLASMA_dsyr2k_Tile_c
      interface
         function PLASMA_dsyr2k_Tile_c(uplo,trans,alpha,A,B,beta,C) &
          & bind(c, name='PLASMA_dsyr2k_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsyr2k_Tile_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            real(kind=c_double), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
         end function PLASMA_dsyr2k_Tile_c
      end interface

!      private PLASMA_dtrmm_Tile_c
      interface
         function PLASMA_dtrmm_Tile_c(side,uplo,transA,diag,alpha,A,B) &
          & bind(c, name='PLASMA_dtrmm_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dtrmm_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            real(kind=c_double), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
         end function PLASMA_dtrmm_Tile_c
      end interface

!      private PLASMA_dtrsm_Tile_c
      interface
         function PLASMA_dtrsm_Tile_c(side,uplo,transA,diag,alpha,A,B) &
          & bind(c, name='PLASMA_dtrsm_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dtrsm_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            real(kind=c_double), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
         end function PLASMA_dtrsm_Tile_c
      end interface

!      private PLASMA_dtrsmpl_Tile_c
      interface
         function PLASMA_dtrsmpl_Tile_c(A,L,IPIV,B) &
          & bind(c, name='PLASMA_dtrsmpl_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dtrsmpl_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
         end function PLASMA_dtrsmpl_Tile_c
      end interface

!      private PLASMA_dtrtri_Tile_c
      interface
         function PLASMA_dtrtri_Tile_c(uplo,diag,A) &
          & bind(c, name='PLASMA_dtrtri_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dtrtri_Tile_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: diag
            type(c_ptr), value :: A
         end function PLASMA_dtrtri_Tile_c
      end interface

!      private PLASMA_dorglq_Tile_c
      interface
         function PLASMA_dorglq_Tile_c(A,T,B) &
          & bind(c, name='PLASMA_dorglq_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dorglq_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
         end function PLASMA_dorglq_Tile_c
      end interface

!      private PLASMA_dorgqr_Tile_c
      interface
         function PLASMA_dorgqr_Tile_c(A,T,B) &
          & bind(c, name='PLASMA_dorgqr_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dorgqr_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
         end function PLASMA_dorgqr_Tile_c
      end interface

!      private PLASMA_dormlq_Tile_c
      interface
         function PLASMA_dormlq_Tile_c(side,trans,A,T,B) &
          & bind(c, name='PLASMA_dormlq_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dormlq_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
         end function PLASMA_dormlq_Tile_c
      end interface

!      private PLASMA_dormqr_Tile_c
      interface
         function PLASMA_dormqr_Tile_c(side,trans,A,T,B) &
          & bind(c, name='PLASMA_dormqr_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dormqr_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
         end function PLASMA_dormqr_Tile_c
      end interface

!      private PLASMA_Alloc_Workspace_dgelqf_c
      interface
         function PLASMA_Alloc_Workspace_dgelqf_c(M,N,T) &
          & bind(c, name='PLASMA_Alloc_Workspace_dgelqf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_dgelqf_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: T ! T is **, so pass by reference
         end function PLASMA_Alloc_Workspace_dgelqf_c
      end interface

!      private PLASMA_Alloc_Workspace_dgels_c
      interface
         function PLASMA_Alloc_Workspace_dgels_c(M,N,T) &
          & bind(c, name='PLASMA_Alloc_Workspace_dgels')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_dgels_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: T ! T is **, so pass by reference
         end function PLASMA_Alloc_Workspace_dgels_c
      end interface

!      private PLASMA_Alloc_Workspace_dgeqrf_c
      interface
         function PLASMA_Alloc_Workspace_dgeqrf_c(M,N,T) &
          & bind(c, name='PLASMA_Alloc_Workspace_dgeqrf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_dgeqrf_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: T ! T is **, so pass by reference
         end function PLASMA_Alloc_Workspace_dgeqrf_c
      end interface

!      private PLASMA_Alloc_Workspace_dgesv_incpiv_c
      interface
         function PLASMA_Alloc_Workspace_dgesv_incpiv_c(N,L,IPIV) &
          & bind(c, name='PLASMA_Alloc_Workspace_dgesv_incpiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_dgesv_incpiv_c
            integer(kind=c_int), value :: N
            type(c_ptr) :: L ! L is **, so pass by reference
            type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         end function PLASMA_Alloc_Workspace_dgesv_incpiv_c
      end interface

!      private PLASMA_Alloc_Workspace_dgetrf_incpiv_c
      interface
         function PLASMA_Alloc_Workspace_dgetrf_incpiv_c(M,N,L,IPIV) &
          & bind(c, name='PLASMA_Alloc_Workspace_dgetrf_incpiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_dgetrf_incpiv_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: L ! L is **, so pass by reference
            type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         end function PLASMA_Alloc_Workspace_dgetrf_incpiv_c
      end interface

!      private PLASMA_Alloc_Workspace_dgeev_c
      interface
         function PLASMA_Alloc_Workspace_dgeev_c(N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_dgeev')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_dgeev_c
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_dgeev_c
      end interface

!      private PLASMA_Alloc_Workspace_dgebrd_c
      interface
         function PLASMA_Alloc_Workspace_dgebrd_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_dgebrd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_dgebrd_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_dgebrd_c
      end interface

!      private PLASMA_Alloc_Workspace_dgesvd_c
      interface
         function PLASMA_Alloc_Workspace_dgesvd_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_dgesvd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_dgesvd_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_dgesvd_c
      end interface

!      private PLASMA_Alloc_Workspace_dsyev_c
      interface
         function PLASMA_Alloc_Workspace_dsyev_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_dsyev')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_dsyev_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_dsyev_c
      end interface

!      private PLASMA_Alloc_Workspace_dsygv_c
      interface
         function PLASMA_Alloc_Workspace_dsygv_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_dsygv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_dsygv_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_dsygv_c
      end interface

!      private PLASMA_Alloc_Workspace_dsytrd_c
      interface
         function PLASMA_Alloc_Workspace_dsytrd_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_dsytrd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_dsytrd_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_dsytrd_c
      end interface

!      private PLASMA_Alloc_Workspace_dgelqf_Tile_c
      interface
         function PLASMA_Alloc_Workspace_dgelqf_Tile_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_dgelqf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_dgelqf_Tile_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_dgelqf_Tile_c
      end interface

!      private PLASMA_Alloc_Workspace_dgels_Tile_c
      interface
         function PLASMA_Alloc_Workspace_dgels_Tile_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_dgels_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_dgels_Tile_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_dgels_Tile_c
      end interface

!      private PLASMA_Alloc_Workspace_dgeqrf_Tile_c
      interface
         function PLASMA_Alloc_Workspace_dgeqrf_Tile_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_dgeqrf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_dgeqrf_Tile_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_dgeqrf_Tile_c
      end interface

!      private PLASMA_Alloc_Workspace_dgesv_incpiv_Tile_c
      interface
         function PLASMA_Alloc_Workspace_dgesv_incpiv_Tile_c(N,descL,IPIV) &
          & bind(c, name='PLASMA_Alloc_Workspace_dgesv_incpiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_dgesv_incpiv_Tile_c
            integer(kind=c_int), value :: N
            type(c_ptr) :: descL ! descL is **, so pass by reference
            type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         end function PLASMA_Alloc_Workspace_dgesv_incpiv_Tile_c
      end interface

!      private PLASMA_Alloc_Workspace_dgetrf_incpiv_Tile_c
      interface
         function PLASMA_Alloc_Workspace_dgetrf_incpiv_Tile_c(N,descL,IPIV) &
          & bind(c, name='PLASMA_Alloc_Workspace_dgetrf_incpiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_dgetrf_incpiv_Tile_c
            integer(kind=c_int), value :: N
            type(c_ptr) :: descL ! descL is **, so pass by reference
            type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         end function PLASMA_Alloc_Workspace_dgetrf_incpiv_Tile_c
      end interface

!      private PLASMA_dLapack_to_Tile_c
      interface
         function PLASMA_dLapack_to_Tile_c(Af77,LDA,A) &
          & bind(c, name='PLASMA_dLapack_to_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dLapack_to_Tile_c
            type(c_ptr), value :: Af77
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: A
         end function PLASMA_dLapack_to_Tile_c
      end interface

!      private PLASMA_dTile_to_Lapack_c
      interface
         function PLASMA_dTile_to_Lapack_c(A,Af77,LDA) &
          & bind(c, name='PLASMA_dTile_to_Lapack')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dTile_to_Lapack_c
            type(c_ptr), value :: A
            type(c_ptr), value :: Af77
            integer(kind=c_int), value :: LDA
         end function PLASMA_dTile_to_Lapack_c
      end interface

  contains

      subroutine PLASMA_dgebrd(jobu,jobvt,M,N,A,LDA,D,E,U,LDU,VT,LDVT,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDU
         integer(kind=c_int), intent(in) :: LDVT
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: jobu
         integer(kind=c_int), intent(in) :: jobvt
         real(kind=c_double), intent(inout), target :: A(LDA,*)
         real(kind=c_double), intent(out), target :: D(*)
         real(kind=c_double), intent(out), target :: E(*)
         real(kind=c_double), intent(out), target :: U(LDU,*)
         real(kind=c_double), intent(out), target :: VT(LDVT,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dgebrd_c(jobu,jobvt,M,N,c_loc(A),LDA,c_loc(D),c_loc(E),c_loc(U),LDU,c_loc(VT),LDVT,T)
      end subroutine PLASMA_dgebrd

      subroutine PLASMA_dgelqf(M,N,A,LDA,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         real(kind=c_double), intent(inout), target :: A(LDA,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dgelqf_c(M,N,c_loc(A),LDA,T)
      end subroutine PLASMA_dgelqf

      subroutine PLASMA_dgelqs(M,N,NRHS,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         real(kind=c_double), intent(in), target :: A(LDA,*)
         real(kind=c_double), intent(inout), target :: B(LDB,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dgelqs_c(M,N,NRHS,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_dgelqs

      subroutine PLASMA_dgels(trans,M,N,NRHS,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: trans
         real(kind=c_double), intent(inout), target :: A(LDA,*)
         real(kind=c_double), intent(inout), target :: B(LDB,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dgels_c(trans,M,N,NRHS,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_dgels

      subroutine PLASMA_dgemm(transA,transB,M,N,K,alpha,A,LDA,B,LDB,beta,C,LDC,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDC
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: transB
         real(kind=c_double), intent(in) :: alpha
         real(kind=c_double), intent(in) :: beta
         real(kind=c_double), intent(in), target :: A(LDA,*)
         real(kind=c_double), intent(in), target :: B(LDB,*)
         real(kind=c_double), intent(inout), target :: C(LDC,*)
         info = PLASMA_dgemm_c(transA,transB,M,N,K,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
      end subroutine PLASMA_dgemm

      subroutine PLASMA_dgeqrf(M,N,A,LDA,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         real(kind=c_double), intent(inout), target :: A(LDA,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dgeqrf_c(M,N,c_loc(A),LDA,T)
      end subroutine PLASMA_dgeqrf

      subroutine PLASMA_dgeqrs(M,N,NRHS,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         real(kind=c_double), intent(inout), target :: A(LDA,*)
         real(kind=c_double), intent(inout), target :: B(LDB,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dgeqrs_c(M,N,NRHS,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_dgeqrs

      subroutine PLASMA_dgesv(N,NRHS,A,LDA,IPIV,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(out), target :: IPIV(*)
         real(kind=c_double), intent(inout), target :: A(LDA,*)
         real(kind=c_double), intent(inout), target :: B(LDB,*)
         info = PLASMA_dgesv_c(N,NRHS,c_loc(A),LDA,c_loc(IPIV),c_loc(B),LDB)
      end subroutine PLASMA_dgesv

      subroutine PLASMA_dgesv_incpiv(N,NRHS,A,LDA,L,IPIV,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         real(kind=c_double), intent(inout), target :: A(LDA,*)
         real(kind=c_double), intent(inout), target :: B(LDB,*)
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dgesv_incpiv_c(N,NRHS,c_loc(A),LDA,L,IPIV,c_loc(B),LDB)
      end subroutine PLASMA_dgesv_incpiv

      subroutine PLASMA_dgesvd(jobu,jobvt,M,N,A,LDA,S,U,LDU,VT,LDVT,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDU
         integer(kind=c_int), intent(in) :: LDVT
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: jobu
         integer(kind=c_int), intent(in) :: jobvt
         real(kind=c_double), intent(inout), target :: A(LDA,*)
         real(kind=c_double), intent(out), target :: S(*)
         real(kind=c_double), intent(out), target :: U(LDU,*)
         real(kind=c_double), intent(out), target :: VT(LDVT,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dgesvd_c(jobu,jobvt,M,N,c_loc(A),LDA,c_loc(S),c_loc(U),LDU,c_loc(VT),LDVT,T)
      end subroutine PLASMA_dgesvd

      subroutine PLASMA_dgetrf(M,N,A,LDA,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(out), target :: IPIV(*)
         real(kind=c_double), intent(inout), target :: A(LDA,*)
         info = PLASMA_dgetrf_c(M,N,c_loc(A),LDA,c_loc(IPIV))
      end subroutine PLASMA_dgetrf

      subroutine PLASMA_dgetrf_incpiv(M,N,A,LDA,L,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         real(kind=c_double), intent(inout), target :: A(LDA,*)
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L    ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dgetrf_incpiv_c(M,N,c_loc(A),LDA,L,IPIV)
      end subroutine PLASMA_dgetrf_incpiv

      subroutine PLASMA_dgetrs(trans,N,NRHS,A,LDA,IPIV,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in), target :: IPIV(*)
         real(kind=c_double), intent(in), target :: A(LDA,*)
         real(kind=c_double), intent(inout), target :: B(LDB,*)
         info = PLASMA_dgetrs_c(trans,N,NRHS,c_loc(A),LDA,c_loc(IPIV),c_loc(B),LDB)
      end subroutine PLASMA_dgetrs

      subroutine PLASMA_dgetrs_incpiv(trans,N,NRHS,A,LDA,L,IPIV,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: trans
         real(kind=c_double), intent(in), target :: A(LDA,*)
         real(kind=c_double), intent(inout), target :: B(LDB,*)
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L    ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dgetrs_incpiv_c(trans,N,NRHS,c_loc(A),LDA,L,IPIV,c_loc(B),LDB)
      end subroutine PLASMA_dgetrs_incpiv

#if defined(PRECISION_z) || defined(PRECISION_c)
      subroutine PLASMA_dsymm(side,uplo,M,N,alpha,A,LDA,B,LDB,beta,C,LDC,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDC
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(in) :: alpha
         real(kind=c_double), intent(in) :: beta
         real(kind=c_double), intent(in), target :: A(LDA,*)
         real(kind=c_double), intent(in), target :: B(LDB,*)
         real(kind=c_double), intent(inout), target :: C(LDC,*)
         info = PLASMA_dsymm_c(side,uplo,M,N,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
      end subroutine PLASMA_dsymm

      subroutine PLASMA_dsyrk(uplo,trans,N,K,alpha,A,LDA,beta,C,LDC,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDC
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(in), target :: A(LDA,*)
         real(kind=c_double), intent(inout), target :: C(LDC,*)
         real(kind=c_double), intent(in) :: alpha
         real(kind=c_double), intent(in) :: beta
         info = PLASMA_dsyrk_c(uplo,trans,N,K,alpha,c_loc(A),LDA,beta,c_loc(C),LDC)
      end subroutine PLASMA_dsyrk

      subroutine PLASMA_dsyr2k(uplo,trans,N,K,alpha,A,LDA,B,LDB,beta,C,LDC,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDC
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(in) :: alpha
         real(kind=c_double), intent(in), target :: A(LDA,*)
         real(kind=c_double), intent(in), target :: B(LDB,*)
         real(kind=c_double), intent(inout), target :: C(LDC,*)
         real(kind=c_double), intent(in) :: beta
         info = PLASMA_dsyr2k_c(uplo,trans,N,K,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
      end subroutine PLASMA_dsyr2k
#endif

      subroutine PLASMA_dsyev(jobz,uplo,N,A,LDA,W,T,Q,LDQ,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDQ
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(inout), target :: A(LDA,*)
         real(kind=c_double), intent(out), target :: W(*)
         real(kind=c_double), intent(out), target :: Q(LDQ,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dsyev_c(jobz,uplo,N,c_loc(A),LDA,c_loc(W),T,c_loc(Q),LDQ)
      end subroutine PLASMA_dsyev

      subroutine PLASMA_dsygv(itype,jobz,uplo,N,A,LDA,B,LDB,W,T,Q,LDQ,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDQ
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(inout), target :: A(LDA,*)
         real(kind=c_double), intent(inout), target :: B(LDB,*)
         real(kind=c_double), intent(out), target :: W(*)
         real(kind=c_double), intent(out), target :: Q(LDQ,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dsygv_c(itype,jobz,uplo,N,c_loc(A),LDA,c_loc(B),LDB,c_loc(W),T,c_loc(Q),LDQ)
      end subroutine PLASMA_dsygv

      subroutine PLASMA_dsygst(itype,uplo,N,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(in), target :: B(LDB,*)
         real(kind=c_double), intent(inout), target :: A(LDA,*)
         info = PLASMA_dsygst_c(itype,uplo,N,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine PLASMA_dsygst

      subroutine PLASMA_dsytrd(jobz,uplo,N,A,LDA,D,E,T,Q,LDQ,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDQ
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(inout), target :: A(LDA,*)
         real(kind=c_double), intent(out), target :: D(*)
         real(kind=c_double), intent(out), target :: E(*)
         real(kind=c_double), intent(out), target :: Q(LDQ,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dsytrd_c(jobz,uplo,N,c_loc(A),LDA,c_loc(D),c_loc(E),T,c_loc(Q),LDQ)
      end subroutine PLASMA_dsytrd

      function PLASMA_dlange(norm,M,N,A,LDA,work)
         use iso_c_binding
         implicit none
         real(kind=c_double) :: PLASMA_dlange
         real(kind=c_double), intent(inout), target :: work(*)
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         real(kind=c_double), intent(in), target :: A(LDA,*)
         PLASMA_dlange = PLASMA_dlange_c(norm,M,N,c_loc(A),LDA,c_loc(work))
      end function PLASMA_dlange

#if defined(PRECISION_z) || defined(PRECISION_c)
      function PLASMA_dlansy(norm,uplo,N,A,LDA,work)
         use iso_c_binding
         implicit none
         real(kind=c_double) :: PLASMA_dlansy
         real(kind=c_double), intent(inout), target :: work(*)
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         real(kind=c_double), intent(in), target :: A(LDA,*)
         PLASMA_dlansy = PLASMA_dlansy_c(norm,uplo,N,c_loc(A),LDA,c_loc(work))
      end function PLASMA_dlansy
#endif

      function PLASMA_dlansy(norm,uplo,N,A,LDA,work)
         use iso_c_binding
         implicit none
         real(kind=c_double) :: PLASMA_dlansy
         real(kind=c_double), intent(inout), target :: work(*)
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         real(kind=c_double), intent(inout), target :: A(LDA,*)
         PLASMA_dlansy = PLASMA_dlansy_c(norm,uplo,N,c_loc(A),LDA,c_loc(work))
      end function PLASMA_dlansy

      subroutine PLASMA_dlaswp(N,A,LDA,K1,K2,IPIV,INCX,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: INCX
         integer(kind=c_int), intent(in) :: K1
         integer(kind=c_int), intent(in) :: K2
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in), target :: IPIV(*)
         real(kind=c_double), intent(inout), target :: A(LDA,*)
         info = PLASMA_dlaswp_c(N,c_loc(A),LDA,K1,K2,c_loc(IPIV),INCX)
      end subroutine PLASMA_dlaswp

      subroutine PLASMA_dlauum(uplo,N,A,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(inout), target :: A(LDA,*)
         info = PLASMA_dlauum_c(uplo,N,c_loc(A),LDA)
      end subroutine PLASMA_dlauum

      subroutine PLASMA_dposv(uplo,N,NRHS,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(inout), target :: A(LDA,*)
         real(kind=c_double), intent(inout), target :: B(LDB,*)
         info = PLASMA_dposv_c(uplo,N,NRHS,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine PLASMA_dposv

      subroutine PLASMA_dpotrf(uplo,N,A,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(inout), target :: A(LDA,*)
         info = PLASMA_dpotrf_c(uplo,N,c_loc(A),LDA)
      end subroutine PLASMA_dpotrf

      subroutine PLASMA_dpotri(uplo,N,A,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(inout), target :: A(LDA,*)
         info = PLASMA_dpotri_c(uplo,N,c_loc(A),LDA)
      end subroutine PLASMA_dpotri

      subroutine PLASMA_dpotrs(uplo,N,NRHS,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(in), target :: A(LDA,*)
         real(kind=c_double), intent(inout), target :: B(LDB,*)
         info = PLASMA_dpotrs_c(uplo,N,NRHS,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine PLASMA_dpotrs

      subroutine PLASMA_dsymm(side,uplo,M,N,alpha,A,LDA,B,LDB,beta,C,LDC,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDC
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(in) :: alpha
         real(kind=c_double), intent(in) :: beta
         real(kind=c_double), intent(in), target :: A(LDA,*)
         real(kind=c_double), intent(in), target :: B(LDB,*)
         real(kind=c_double), intent(inout), target :: C(LDC,*)
         info = PLASMA_dsymm_c(side,uplo,M,N,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
      end subroutine PLASMA_dsymm

      subroutine PLASMA_dsyrk(uplo,trans,N,K,alpha,A,LDA,beta,C,LDC,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDC
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(in) :: alpha
         real(kind=c_double), intent(in) :: beta
         real(kind=c_double), intent(in), target :: A(LDA,*)
         real(kind=c_double), intent(inout), target :: C(LDC,*)
         info = PLASMA_dsyrk_c(uplo,trans,N,K,alpha,c_loc(A),LDA,beta,c_loc(C),LDC)
      end subroutine PLASMA_dsyrk

      subroutine PLASMA_dsyr2k(uplo,trans,N,K,alpha,A,LDA,B,LDB,beta,C,LDC,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDC
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(in) :: alpha
         real(kind=c_double), intent(in) :: beta
         real(kind=c_double), intent(in), target :: A(LDA,*)
         real(kind=c_double), intent(in), target :: B(LDB,*)
         real(kind=c_double), intent(inout), target :: C(LDC,*)
         info = PLASMA_dsyr2k_c(uplo,trans,N,K,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
      end subroutine PLASMA_dsyr2k

      subroutine PLASMA_dtrmm(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(in) :: alpha
         real(kind=c_double), intent(in), target :: A(LDA,*)
         real(kind=c_double), intent(inout), target :: B(LDB,*)
         info = PLASMA_dtrmm_c(side,uplo,transA,diag,N,NRHS,alpha,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine PLASMA_dtrmm

      subroutine PLASMA_dtrsm(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(in) :: alpha
         real(kind=c_double), intent(in), target :: A(LDA,*)
         real(kind=c_double), intent(inout), target :: B(LDB,*)
         info = PLASMA_dtrsm_c(side,uplo,transA,diag,N,NRHS,alpha,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine PLASMA_dtrsm

      subroutine PLASMA_dtrsmpl(N,NRHS,A,LDA,L,IPIV,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         real(kind=c_double), intent(in), target :: A(LDA,*)
         type(c_ptr), value :: L    ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         real(kind=c_double), intent(inout), target :: B(LDB,*)
         info = PLASMA_dtrsmpl_c(N,NRHS,c_loc(A),LDA,L,IPIV,c_loc(B),LDB)
      end subroutine PLASMA_dtrsmpl

      subroutine PLASMA_dtrtri(uplo,diag,N,A,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(inout), target :: A(LDA,*)
         info = PLASMA_dtrtri_c(uplo,diag,N,c_loc(A),LDA)
      end subroutine PLASMA_dtrtri

      subroutine PLASMA_dorglq(M,N,K,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         real(kind=c_double), intent(inout), target :: A(LDA,*)
         real(kind=c_double), intent(out), target :: B(LDB,*)
         info = PLASMA_dorglq_c(M,N,K,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_dorglq

      subroutine PLASMA_dorgqr(M,N,K,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         real(kind=c_double), intent(inout), target :: A(LDA,*)
         real(kind=c_double), intent(out), target :: B(LDB,*)
         info = PLASMA_dorgqr_c(M,N,K,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_dorgqr

      subroutine PLASMA_dormlq(side,trans,M,N,K,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: trans
         real(kind=c_double), intent(in), target :: A(LDA,*)
         real(kind=c_double), intent(inout), target :: B(LDB,*)
         info = PLASMA_dormlq_c(side,trans,M,N,K,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_dormlq

      subroutine PLASMA_dormqr(side,trans,M,N,K,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: trans
         real(kind=c_double), intent(in), target :: A(LDA,*)
         real(kind=c_double), intent(inout), target :: B(LDB,*)
         info = PLASMA_dormqr_c(side,trans,M,N,K,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_dormqr

      subroutine PLASMA_dgecfi(m,n,A,fin,imb,inb,fout,omb,onb,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         real(kind=c_double), intent(inout), target :: A(*)
         integer(kind=c_int), intent(in) :: fin
         integer(kind=c_int), intent(in) :: fout
         integer(kind=c_int), intent(in) :: imb
         integer(kind=c_int), intent(in) :: inb
         integer(kind=c_int), intent(in) :: omb
         integer(kind=c_int), intent(in) :: onb
         integer(kind=c_int), intent(in) :: m
         integer(kind=c_int), intent(in) :: n
         info = PLASMA_dgecfi_c(m,n,c_loc(A),fin,imb,inb,fout,omb,onb)
      end subroutine PLASMA_dgecfi

      subroutine PLASMA_dgetmi(m,n,A,fin,mb,nb,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         real(kind=c_double), intent(inout), target :: A(*)
         integer(kind=c_int), intent(in) :: fin
         integer(kind=c_int), intent(in) :: mb
         integer(kind=c_int), intent(in) :: nb
         integer(kind=c_int), intent(in) :: m
         integer(kind=c_int), intent(in) :: n
         info = PLASMA_dgetmi_c(m,n,c_loc(A),fin,mb,nb)
      end subroutine PLASMA_dgetmi

      subroutine PLASMA_dgebrd_Tile(jobu,jobvt,A,D,E,U,VT,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobu
         integer(kind=c_int), intent(in) :: jobvt
         real(kind=c_double), intent(out), target :: D(*)
         real(kind=c_double), intent(out), target :: E(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: U ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: VT ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dgebrd_Tile_c(jobu,jobvt,A,c_loc(D),c_loc(E),U,VT,T)
      end subroutine PLASMA_dgebrd_Tile

      subroutine PLASMA_dgelqf_Tile(A,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dgelqf_Tile_c(A,T)
      end subroutine PLASMA_dgelqf_Tile

      subroutine PLASMA_dgelqs_Tile(A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dgelqs_Tile_c(A,T,B)
      end subroutine PLASMA_dgelqs_Tile

      subroutine PLASMA_dgels_Tile(trans,A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dgels_Tile_c(trans,A,T,B)
      end subroutine PLASMA_dgels_Tile

      subroutine PLASMA_dgemm_Tile(transA,transB,alpha,A,B,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: transB
         real(kind=c_double), intent(in) :: alpha
         real(kind=c_double), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dgemm_Tile_c(transA,transB,alpha,A,B,beta,C)
      end subroutine PLASMA_dgemm_Tile

      subroutine PLASMA_dgeqrf_Tile(A,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dgeqrf_Tile_c(A,T)
      end subroutine PLASMA_dgeqrf_Tile

      subroutine PLASMA_dgeqrs_Tile(A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dgeqrs_Tile_c(A,T,B)
      end subroutine PLASMA_dgeqrs_Tile

      subroutine PLASMA_dgesv_Tile(A,IPIV,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: IPIV(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dgesv_Tile_c(A,c_loc(IPIV),B)
      end subroutine PLASMA_dgesv_Tile

      subroutine PLASMA_dgesv_incpiv_Tile(A,L,IPIV,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dgesv_incpiv_Tile_c(A,L,IPIV,B)
      end subroutine PLASMA_dgesv_incpiv_Tile

      subroutine PLASMA_dgesvd_Tile(jobu,jobvt,A,S,U,VT,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobu
         integer(kind=c_int), intent(in) :: jobvt
         real(kind=c_double), intent(out), target :: S(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: U ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: VT ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dgesvd_Tile_c(jobu,jobvt,A,c_loc(S),U,VT,T)
      end subroutine PLASMA_dgesvd_Tile

      subroutine PLASMA_dgetrf_Tile(A,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: IPIV(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dgetrf_Tile_c(A,c_loc(IPIV))
      end subroutine PLASMA_dgetrf_Tile

      subroutine PLASMA_dgetrf_incpiv_Tile(A,L,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dgetrf_incpiv_Tile_c(A,L,IPIV)
      end subroutine PLASMA_dgetrf_incpiv_Tile

      subroutine PLASMA_dgetrs_Tile(trans,A,IPIV,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         integer(kind=c_int), intent(in), target :: IPIV(*)
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dgetrs_Tile_c(trans,A,c_loc(IPIV),B)
      end subroutine PLASMA_dgetrs_Tile

      subroutine PLASMA_dgetrs_incpiv_Tile(A,L,IPIV,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dgetrs_incpiv_Tile_c(A,L,IPIV,B)
      end subroutine PLASMA_dgetrs_incpiv_Tile

#if defined(PRECISION_z) || defined(PRECISION_c)
      subroutine PLASMA_dsymm_Tile(side,uplo,alpha,A,B,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(in) :: alpha
         real(kind=c_double), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dsymm_Tile_c(side,uplo,alpha,A,B,beta,C)
      end subroutine PLASMA_dsymm_Tile

      subroutine PLASMA_dsyrk_Tile(uplo,trans,alpha,A,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(in) :: alpha
         real(kind=c_double), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dsyrk_Tile_c(uplo,trans,alpha,A,beta,C)
      end subroutine PLASMA_dsyrk_Tile

      subroutine PLASMA_dsyr2k_Tile(uplo,trans,alpha,A,B,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(in) :: alpha
         real(kind=c_double), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dsyr2k_Tile_c(uplo,trans,alpha,A,B,beta,C)
      end subroutine PLASMA_dsyr2k_Tile
#endif

      subroutine PLASMA_dsyev_Tile(jobz,uplo,A,W,T,Q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(out), target :: W(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: Q ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dsyev_Tile_c(jobz,uplo,A,c_loc(W),T,Q)
      end subroutine PLASMA_dsyev_Tile

      subroutine PLASMA_dsygv_Tile(itype,jobz,uplo,A,B,W,T,Q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(out), target :: W(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: Q ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dsygv_Tile_c(itype,jobz,uplo,A,B,c_loc(W),T,Q)
      end subroutine PLASMA_dsygv_Tile

      subroutine PLASMA_dsygst_Tile(itype,uplo,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dsygst_Tile_c(itype,uplo,A,B)
      end subroutine PLASMA_dsygst_Tile

      subroutine PLASMA_dsytrd_Tile(jobz,uplo,A,D,E,T,Q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(out), target :: D(*)
         real(kind=c_double), intent(out), target :: E(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: Q ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dsytrd_Tile_c(jobz,uplo,A,c_loc(D),c_loc(E),T,Q)
      end subroutine PLASMA_dsytrd_Tile

      function PLASMA_dlange_Tile(norm,A,work)
         use iso_c_binding
         implicit none
         real(kind=c_double) :: PLASMA_dlange_Tile
         integer(kind=c_int), intent(in) :: norm
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         real(kind=c_double), intent(inout), target :: work(*)
         PLASMA_dlange_Tile = PLASMA_dlange_Tile_c(norm,A,c_loc(work))
       end function PLASMA_dlange_Tile

#if defined(PRECISION_z) || defined(PRECISION_c)
      function PLASMA_dlansy_Tile(norm,uplo,A,work)
         use iso_c_binding
         implicit none
         real(kind=c_double) :: PLASMA_dlansy_Tile
         real(kind=c_double), intent(inout), target :: work(*)
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         PLASMA_dlansy_Tile = PLASMA_dlansy_Tile_c(norm,uplo,A,c_loc(work))
      end function PLASMA_dlansy_Tile
#endif

      function PLASMA_dlansy_Tile(norm,uplo,A,work)
         use iso_c_binding
         implicit none
         real(kind=c_double) :: PLASMA_dlansy_Tile
         real(kind=c_double), intent(inout), target :: work(*)
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         PLASMA_dlansy_Tile = PLASMA_dlansy_Tile_c(norm,uplo,A,c_loc(work))
      end function PLASMA_dlansy_Tile

      subroutine PLASMA_dlaswp_Tile(A,K1,K2,IPIV,INCX,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: INCX
         integer(kind=c_int), intent(in) :: K1
         integer(kind=c_int), intent(in) :: K2
         integer(kind=c_int), intent(in), target :: IPIV(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dlaswp_Tile_c(A,K1,K2,c_loc(IPIV),INCX)
      end subroutine PLASMA_dlaswp_Tile

      subroutine PLASMA_dlauum_Tile(uplo,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dlauum_Tile_c(uplo,A)
      end subroutine PLASMA_dlauum_Tile

      subroutine PLASMA_dposv_Tile(uplo,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dposv_Tile_c(uplo,A,B)
      end subroutine PLASMA_dposv_Tile

      subroutine PLASMA_dpotrf_Tile(uplo,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dpotrf_Tile_c(uplo,A)
      end subroutine PLASMA_dpotrf_Tile

      subroutine PLASMA_dpotri_Tile(uplo,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dpotri_Tile_c(uplo,A)
      end subroutine PLASMA_dpotri_Tile

      subroutine PLASMA_dpotrs_Tile(uplo,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dpotrs_Tile_c(uplo,A,B)
      end subroutine PLASMA_dpotrs_Tile

      subroutine PLASMA_dsymm_Tile(side,uplo,alpha,A,B,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(in) :: alpha
         real(kind=c_double), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dsymm_Tile_c(side,uplo,alpha,A,B,beta,C)
      end subroutine PLASMA_dsymm_Tile

      subroutine PLASMA_dsyrk_Tile(uplo,trans,alpha,A,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(in) :: alpha
         real(kind=c_double), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dsyrk_Tile_c(uplo,trans,alpha,A,beta,C)
      end subroutine PLASMA_dsyrk_Tile

      subroutine PLASMA_dsyr2k_Tile(uplo,trans,alpha,A,B,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(in) :: alpha
         real(kind=c_double), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dsyr2k_Tile_c(uplo,trans,alpha,A,B,beta,C)
      end subroutine PLASMA_dsyr2k_Tile

      subroutine PLASMA_dtrmm_Tile(side,uplo,transA,diag,alpha,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(in) :: alpha
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dtrmm_Tile_c(side,uplo,transA,diag,alpha,A,B)
      end subroutine PLASMA_dtrmm_Tile

      subroutine PLASMA_dtrsm_Tile(side,uplo,transA,diag,alpha,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(in) :: alpha
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dtrsm_Tile_c(side,uplo,transA,diag,alpha,A,B)
      end subroutine PLASMA_dtrsm_Tile

      subroutine PLASMA_dtrsmpl_Tile(A,L,IPIV,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dtrsmpl_Tile_c(A,L,IPIV,B)
      end subroutine PLASMA_dtrsmpl_Tile

      subroutine PLASMA_dtrtri_Tile(uplo,diag,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dtrtri_Tile_c(uplo,diag,A)
      end subroutine PLASMA_dtrtri_Tile

      subroutine PLASMA_dorglq_Tile(A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dorglq_Tile_c(A,T,B)
      end subroutine PLASMA_dorglq_Tile

      subroutine PLASMA_dorgqr_Tile(A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dorgqr_Tile_c(A,T,B)
      end subroutine PLASMA_dorgqr_Tile

      subroutine PLASMA_dormlq_Tile(side,trans,A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dormlq_Tile_c(side,trans,A,T,B)
      end subroutine PLASMA_dormlq_Tile

      subroutine PLASMA_dormqr_Tile(side,trans,A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dormqr_Tile_c(side,trans,A,T,B)
      end subroutine PLASMA_dormqr_Tile

      subroutine PLASMA_Alloc_Workspace_dgelqf(M,N,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         info = PLASMA_Alloc_Workspace_dgelqf_c(M,N,T)
      end subroutine PLASMA_Alloc_Workspace_dgelqf

      subroutine PLASMA_Alloc_Workspace_dgels(M,N,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         info = PLASMA_Alloc_Workspace_dgels_c(M,N,T)
      end subroutine PLASMA_Alloc_Workspace_dgels

      subroutine PLASMA_Alloc_Workspace_dgeqrf(M,N,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         info = PLASMA_Alloc_Workspace_dgeqrf_c(M,N,T)
      end subroutine PLASMA_Alloc_Workspace_dgeqrf

      subroutine PLASMA_Alloc_Workspace_dgesv_incpiv(N,L,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         type(c_ptr) :: L ! L is **, so pass by reference
         info = PLASMA_Alloc_Workspace_dgesv_incpiv_c(N,L,IPIV)
      end subroutine PLASMA_Alloc_Workspace_dgesv_incpiv

      subroutine PLASMA_Alloc_Workspace_dgetrf_incpiv(M,N,L,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         type(c_ptr) :: L ! L is **, so pass by reference
         info = PLASMA_Alloc_Workspace_dgetrf_incpiv_c(M,N,L,IPIV)
      end subroutine PLASMA_Alloc_Workspace_dgetrf_incpiv

      subroutine PLASMA_Alloc_Workspace_dgeev(N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_dgeev_c(N,descT)
      end subroutine PLASMA_Alloc_Workspace_dgeev

      subroutine PLASMA_Alloc_Workspace_dgebrd(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_dgebrd_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_dgebrd

      subroutine PLASMA_Alloc_Workspace_dgesvd(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_dgesvd_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_dgesvd

      subroutine PLASMA_Alloc_Workspace_dsyev(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_dsyev_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_dsyev

      subroutine PLASMA_Alloc_Workspace_dsygv(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_dsygv_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_dsygv

      subroutine PLASMA_Alloc_Workspace_dsytrd(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_dsytrd_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_dsytrd

      subroutine PLASMA_Alloc_Workspace_dgelqf_Tile(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_dgelqf_Tile_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_dgelqf_Tile

      subroutine PLASMA_Alloc_Workspace_dgels_Tile(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_dgels_Tile_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_dgels_Tile

      subroutine PLASMA_Alloc_Workspace_dgeqrf_Tile(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_dgeqrf_Tile_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_dgeqrf_Tile

      subroutine PLASMA_Alloc_Workspace_dgesv_incpiv_Tile(N,descL,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         type(c_ptr) :: descL ! descL is **, so pass by reference
         info = PLASMA_Alloc_Workspace_dgesv_incpiv_Tile_c(N,descL,IPIV)
      end subroutine PLASMA_Alloc_Workspace_dgesv_incpiv_Tile

      subroutine PLASMA_Alloc_Workspace_dgetrf_incpiv_Tile(N,descL,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         type(c_ptr) :: descL ! descL is **, so pass by reference
         info = PLASMA_Alloc_Workspace_dgetrf_incpiv_Tile_c(N,descL,IPIV)
      end subroutine PLASMA_Alloc_Workspace_dgetrf_incpiv_Tile

      subroutine PLASMA_dLapack_to_Tile(Af77,LDA,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         real(kind=c_double), intent(in), target :: Af77(LDA,*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dLapack_to_Tile_c(c_loc(Af77),LDA,A)
      end subroutine PLASMA_dLapack_to_Tile

      subroutine PLASMA_dTile_to_Lapack(A,Af77,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         real(kind=c_double), intent(out), target :: Af77(LDA,*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dTile_to_Lapack_c(A,c_loc(Af77),LDA)
      end subroutine PLASMA_dTile_to_Lapack

end module plasma_d
