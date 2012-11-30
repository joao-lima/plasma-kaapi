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
! @file plasma_sf90.F90
!
!  PLASMA fortran 90 interface
!  PLASMA is a software package provided by Univ. of Tennessee,
!  Univ. of California Berkeley and Univ. of Colorado Denver
!
! @version 2.4.2
! @author Numerical Algorithm Group
! @author Mathieu Faverge
! @date 2011-09-15
! @generated s Thu Sep 15 12:09:27 2011
!
#define PRECISION_s

module plasma_s

!      private PLASMA_sgebrd_c
      interface
         function PLASMA_sgebrd_c(jobu,jobvt,M,N,A,LDA,D,E,U,LDU,VT,LDVT,T) &
          & bind(c, name='PLASMA_sgebrd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgebrd_c
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
         end function PLASMA_sgebrd_c
      end interface

!      private PLASMA_sgelqf_c
      interface
         function PLASMA_sgelqf_c(M,N,A,LDA,T) &
          & bind(c, name='PLASMA_sgelqf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgelqf_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: T
         end function PLASMA_sgelqf_c
      end interface

!      private PLASMA_sgelqs_c
      interface
         function PLASMA_sgelqs_c(M,N,NRHS,A,LDA,T,B,LDB) &
          & bind(c, name='PLASMA_sgelqs')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgelqs_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_sgelqs_c
      end interface

!      private PLASMA_sgels_c
      interface
         function PLASMA_sgels_c(trans,M,N,NRHS,A,LDA,T,B,LDB) &
          & bind(c, name='PLASMA_sgels')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgels_c
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_sgels_c
      end interface

!      private PLASMA_sgemm_c
      interface
         function PLASMA_sgemm_c(transA,transB,M,N,K,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='PLASMA_sgemm')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgemm_c
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: transB
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            real(kind=c_float), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            real(kind=c_float), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
         end function PLASMA_sgemm_c
      end interface

!      private PLASMA_sgeqrf_c
      interface
         function PLASMA_sgeqrf_c(M,N,A,LDA,T) &
          & bind(c, name='PLASMA_sgeqrf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgeqrf_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: T
         end function PLASMA_sgeqrf_c
      end interface

!      private PLASMA_sgeqrs_c
      interface
         function PLASMA_sgeqrs_c(M,N,NRHS,A,LDA,T,B,LDB) &
          & bind(c, name='PLASMA_sgeqrs')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgeqrs_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_sgeqrs_c
      end interface

!      private PLASMA_sgesv_c
      interface
         function PLASMA_sgesv_c(N,NRHS,A,LDA,IPIV,B,LDB) &
          & bind(c, name='PLASMA_sgesv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgesv_c
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_sgesv_c
      end interface

!      private PLASMA_sgesv_incpiv_c
      interface
         function PLASMA_sgesv_incpiv_c(N,NRHS,A,LDA,L,IPIV,B,LDB) &
          & bind(c, name='PLASMA_sgesv_incpiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgesv_incpiv_c
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_sgesv_incpiv_c
      end interface

!      private PLASMA_sgesvd_c
      interface
         function PLASMA_sgesvd_c(jobu,jobvt,M,N,A,LDA,S,U,LDU,VT,LDVT,T) &
          & bind(c, name='PLASMA_sgesvd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgesvd_c
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
         end function PLASMA_sgesvd_c
      end interface

!      private PLASMA_sgetrf_c
      interface
         function PLASMA_sgetrf_c(M,N,A,LDA,IPIV) &
          & bind(c, name='PLASMA_sgetrf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgetrf_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: IPIV
         end function PLASMA_sgetrf_c
      end interface

!      private PLASMA_sgetrf_incpiv_c
      interface
         function PLASMA_sgetrf_incpiv_c(M,N,A,LDA,L,IPIV) &
          & bind(c, name='PLASMA_sgetrf_incpiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgetrf_incpiv_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
         end function PLASMA_sgetrf_incpiv_c
      end interface

!      private PLASMA_sgetrs_c
      interface
         function PLASMA_sgetrs_c(trans,N,NRHS,A,LDA,IPIV,B,LDB) &
          & bind(c, name='PLASMA_sgetrs')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgetrs_c
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_sgetrs_c
      end interface

!      private PLASMA_sgetrs_incpiv_c
      interface
         function PLASMA_sgetrs_incpiv_c(trans,N,NRHS,A,LDA,L,IPIV,B,LDB) &
          & bind(c, name='PLASMA_sgetrs_incpiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgetrs_incpiv_c
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_sgetrs_incpiv_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
!      private PLASMA_ssymm_c
      interface
         function PLASMA_ssymm_c(side,uplo,M,N,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='PLASMA_ssymm')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ssymm_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            real(kind=c_float), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            real(kind=c_float), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
         end function PLASMA_ssymm_c
      end interface

!      private PLASMA_ssyrk_c
      interface
         function PLASMA_ssyrk_c(uplo,trans,N,K,alpha,A,LDA,beta,C,LDC) &
          & bind(c, name='PLASMA_ssyrk')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ssyrk_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            real(kind=c_float), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            real(kind=c_float), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
         end function PLASMA_ssyrk_c
      end interface

!      private PLASMA_ssyr2k_c
      interface
         function PLASMA_ssyr2k_c(uplo,trans,N,K,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='PLASMA_ssyr2k')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ssyr2k_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            real(kind=c_float), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            real(kind=c_float), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
         end function PLASMA_ssyr2k_c
      end interface
#endif

!      private PLASMA_ssyev_c
      interface
         function PLASMA_ssyev_c(jobz,uplo,N,A,LDA,W,T,Q,LDQ) &
          & bind(c, name='PLASMA_ssyev')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ssyev_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
         end function PLASMA_ssyev_c
      end interface

!      private PLASMA_ssygv_c
      interface
         function PLASMA_ssygv_c(itype,jobz,uplo,N,A,LDA,B,LDB,W,T,Q,LDQ) &
          & bind(c, name='PLASMA_ssygv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ssygv_c
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
         end function PLASMA_ssygv_c
      end interface

!      private PLASMA_ssygst_c
      interface
         function PLASMA_ssygst_c(itype,uplo,N,A,LDA,B,LDB) &
          & bind(c, name='PLASMA_ssygst')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ssygst_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_ssygst_c
      end interface

!      private PLASMA_ssytrd_c
      interface
         function PLASMA_ssytrd_c(jobz,uplo,N,A,LDA,D,E,T,Q,LDQ) &
          & bind(c, name='PLASMA_ssytrd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ssytrd_c
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
         end function PLASMA_ssytrd_c
      end interface

!      private PLASMA_slange_c
      interface
         function PLASMA_slange_c(norm,M,N,A,LDA,work) &
          & bind(c, name='PLASMA_slange')
            use iso_c_binding
            implicit none
            real(kind=c_float) :: PLASMA_slange_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: work
         end function PLASMA_slange_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
!      private PLASMA_slansy_c
      interface
         function PLASMA_slansy_c(norm,uplo,N,A,LDA,work) &
          & bind(c, name='PLASMA_slansy')
            use iso_c_binding
            implicit none
            real(kind=c_float) :: PLASMA_slansy_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: work
         end function PLASMA_slansy_c
      end interface
#endif

!      private PLASMA_slansy_c
      interface
         function PLASMA_slansy_c(norm,uplo,N,A,LDA,work) &
          & bind(c, name='PLASMA_slansy')
            use iso_c_binding
            implicit none
            real(kind=c_float) :: PLASMA_slansy_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: work
         end function PLASMA_slansy_c
      end interface

!      private PLASMA_slaswp_c
      interface
         function PLASMA_slaswp_c(N,A,LDA,K1,K2,IPIV,INCX) &
          & bind(c, name='PLASMA_slaswp')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_slaswp_c
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            integer(kind=c_int), value :: K1
            integer(kind=c_int), value :: K2
            type(c_ptr), value :: IPIV
            integer(kind=c_int), value :: INCX
         end function PLASMA_slaswp_c
      end interface

!      private PLASMA_slauum_c
      interface
         function PLASMA_slauum_c(uplo,N,A,LDA) &
          & bind(c, name='PLASMA_slauum')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_slauum_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
         end function PLASMA_slauum_c
      end interface

! PLASMA_splgsy uses unsigned ints, not generating Fortran interface

! PLASMA_splgsy uses unsigned ints, not generating Fortran interface

! PLASMA_splrnt uses unsigned ints, not generating Fortran interface

!      private PLASMA_sposv_c
      interface
         function PLASMA_sposv_c(uplo,N,NRHS,A,LDA,B,LDB) &
          & bind(c, name='PLASMA_sposv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sposv_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_sposv_c
      end interface

!      private PLASMA_spotrf_c
      interface
         function PLASMA_spotrf_c(uplo,N,A,LDA) &
          & bind(c, name='PLASMA_spotrf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_spotrf_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
         end function PLASMA_spotrf_c
      end interface

!      private PLASMA_spotri_c
      interface
         function PLASMA_spotri_c(uplo,N,A,LDA) &
          & bind(c, name='PLASMA_spotri')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_spotri_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
         end function PLASMA_spotri_c
      end interface

!      private PLASMA_spotrs_c
      interface
         function PLASMA_spotrs_c(uplo,N,NRHS,A,LDA,B,LDB) &
          & bind(c, name='PLASMA_spotrs')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_spotrs_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_spotrs_c
      end interface

!      private PLASMA_ssymm_c
      interface
         function PLASMA_ssymm_c(side,uplo,M,N,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='PLASMA_ssymm')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ssymm_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            real(kind=c_float), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            real(kind=c_float), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
         end function PLASMA_ssymm_c
      end interface

!      private PLASMA_ssyrk_c
      interface
         function PLASMA_ssyrk_c(uplo,trans,N,K,alpha,A,LDA,beta,C,LDC) &
          & bind(c, name='PLASMA_ssyrk')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ssyrk_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            real(kind=c_float), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            real(kind=c_float), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
         end function PLASMA_ssyrk_c
      end interface

!      private PLASMA_ssyr2k_c
      interface
         function PLASMA_ssyr2k_c(uplo,trans,N,K,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='PLASMA_ssyr2k')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ssyr2k_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            real(kind=c_float), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            real(kind=c_float), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
         end function PLASMA_ssyr2k_c
      end interface

!      private PLASMA_strmm_c
      interface
         function PLASMA_strmm_c(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB) &
          & bind(c, name='PLASMA_strmm')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_strmm_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            real(kind=c_float), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_strmm_c
      end interface

!      private PLASMA_strsm_c
      interface
         function PLASMA_strsm_c(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB) &
          & bind(c, name='PLASMA_strsm')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_strsm_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            real(kind=c_float), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_strsm_c
      end interface

!      private PLASMA_strsmpl_c
      interface
         function PLASMA_strsmpl_c(N,NRHS,A,LDA,L,IPIV,B,LDB) &
          & bind(c, name='PLASMA_strsmpl')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_strsmpl_c
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_strsmpl_c
      end interface

!      private PLASMA_strtri_c
      interface
         function PLASMA_strtri_c(uplo,diag,N,A,LDA) &
          & bind(c, name='PLASMA_strtri')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_strtri_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: diag
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
         end function PLASMA_strtri_c
      end interface

!      private PLASMA_sorglq_c
      interface
         function PLASMA_sorglq_c(M,N,K,A,LDA,T,B,LDB) &
          & bind(c, name='PLASMA_sorglq')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sorglq_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_sorglq_c
      end interface

!      private PLASMA_sorgqr_c
      interface
         function PLASMA_sorgqr_c(M,N,K,A,LDA,T,B,LDB) &
          & bind(c, name='PLASMA_sorgqr')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sorgqr_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_sorgqr_c
      end interface

!      private PLASMA_sormlq_c
      interface
         function PLASMA_sormlq_c(side,trans,M,N,K,A,LDA,T,B,LDB) &
          & bind(c, name='PLASMA_sormlq')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sormlq_c
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
         end function PLASMA_sormlq_c
      end interface

!      private PLASMA_sormqr_c
      interface
         function PLASMA_sormqr_c(side,trans,M,N,K,A,LDA,T,B,LDB) &
          & bind(c, name='PLASMA_sormqr')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sormqr_c
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
         end function PLASMA_sormqr_c
      end interface

!      private PLASMA_sgecfi_c
      interface
         function PLASMA_sgecfi_c(m,n,A,fin,imb,inb,fout,omb,onb) &
          & bind(c, name='PLASMA_sgecfi')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgecfi_c
            integer(kind=c_int), value :: m
            integer(kind=c_int), value :: n
            type(c_ptr), value :: A
            integer(kind=c_int), value :: fin
            integer(kind=c_int), value :: imb
            integer(kind=c_int), value :: inb
            integer(kind=c_int), value :: fout
            integer(kind=c_int), value :: omb
            integer(kind=c_int), value :: onb
         end function PLASMA_sgecfi_c
      end interface

!      private PLASMA_sgetmi_c
      interface
         function PLASMA_sgetmi_c(m,n,A,fin,mb,nb) &
          & bind(c, name='PLASMA_sgetmi')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgetmi_c
            integer(kind=c_int), value :: m
            integer(kind=c_int), value :: n
            type(c_ptr), value :: A
            integer(kind=c_int), value :: fin
            integer(kind=c_int), value :: mb
            integer(kind=c_int), value :: nb
         end function PLASMA_sgetmi_c
      end interface

!      private PLASMA_sgebrd_Tile_c
      interface
         function PLASMA_sgebrd_Tile_c(jobu,jobvt,A,D,E,U,VT,T) &
          & bind(c, name='PLASMA_sgebrd_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgebrd_Tile_c
            integer(kind=c_int), value :: jobu
            integer(kind=c_int), value :: jobvt
            type(c_ptr), value :: A
            type(c_ptr), value :: D
            type(c_ptr), value :: E
            type(c_ptr), value :: U
            type(c_ptr), value :: VT
            type(c_ptr), value :: T
         end function PLASMA_sgebrd_Tile_c
      end interface

!      private PLASMA_sgelqf_Tile_c
      interface
         function PLASMA_sgelqf_Tile_c(A,T) &
          & bind(c, name='PLASMA_sgelqf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgelqf_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
         end function PLASMA_sgelqf_Tile_c
      end interface

!      private PLASMA_sgelqs_Tile_c
      interface
         function PLASMA_sgelqs_Tile_c(A,B,T) &
          & bind(c, name='PLASMA_sgelqs_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgelqs_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: T
         end function PLASMA_sgelqs_Tile_c
      end interface

!      private PLASMA_sgels_Tile_c
      interface
         function PLASMA_sgels_Tile_c(trans,A,T,B) &
          & bind(c, name='PLASMA_sgels_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgels_Tile_c
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: T
         end function PLASMA_sgels_Tile_c
      end interface

!      private PLASMA_sgemm_Tile_c
      interface
         function PLASMA_sgemm_Tile_c(transA,transB,alpha,A,B,beta,C) &
          & bind(c, name='PLASMA_sgemm_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgemm_Tile_c
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: transB
            real(kind=c_float), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            real(kind=c_float), value :: beta
            type(c_ptr), value :: C
         end function PLASMA_sgemm_Tile_c
      end interface

!      private PLASMA_sgeqrf_Tile_c
      interface
         function PLASMA_sgeqrf_Tile_c(A,T) &
          & bind(c, name='PLASMA_sgeqrf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgeqrf_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
         end function PLASMA_sgeqrf_Tile_c
      end interface

!      private PLASMA_sgeqrs_Tile_c
      interface
         function PLASMA_sgeqrs_Tile_c(A,T,B) &
          & bind(c, name='PLASMA_sgeqrs_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgeqrs_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: T
         end function PLASMA_sgeqrs_Tile_c
      end interface

!      private PLASMA_sgesv_Tile_c
      interface
         function PLASMA_sgesv_Tile_c(A,IPIV,B) &
          & bind(c, name='PLASMA_sgesv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgesv_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
         end function PLASMA_sgesv_Tile_c
      end interface

!      private PLASMA_sgesv_incpiv_Tile_c
      interface
         function PLASMA_sgesv_incpiv_Tile_c(A,L,IPIV,B) &
          & bind(c, name='PLASMA_sgesv_incpiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgesv_incpiv_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
         end function PLASMA_sgesv_incpiv_Tile_c
      end interface

!      private PLASMA_sgesvd_Tile_c
      interface
         function PLASMA_sgesvd_Tile_c(jobu,jobvt,A,S,U,VT,T) &
          & bind(c, name='PLASMA_sgesvd_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgesvd_Tile_c
            integer(kind=c_int), value :: jobu
            integer(kind=c_int), value :: jobvt
            type(c_ptr), value :: A
            type(c_ptr), value :: S
            type(c_ptr), value :: U
            type(c_ptr), value :: VT
            type(c_ptr), value :: T
         end function PLASMA_sgesvd_Tile_c
      end interface

!      private PLASMA_sgetrf_Tile_c
      interface
         function PLASMA_sgetrf_Tile_c(A,IPIV) &
          & bind(c, name='PLASMA_sgetrf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgetrf_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
         end function PLASMA_sgetrf_Tile_c
      end interface

!      private PLASMA_sgetrf_incpiv_Tile_c
      interface
         function PLASMA_sgetrf_incpiv_Tile_c(A,L,IPIV) &
          & bind(c, name='PLASMA_sgetrf_incpiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgetrf_incpiv_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
         end function PLASMA_sgetrf_incpiv_Tile_c
      end interface

!      private PLASMA_sgetrs_Tile_c
      interface
         function PLASMA_sgetrs_Tile_c(trans,A,IPIV,B) &
          & bind(c, name='PLASMA_sgetrs_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgetrs_Tile_c
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
         end function PLASMA_sgetrs_Tile_c
      end interface

!      private PLASMA_sgetrs_incpiv_Tile_c
      interface
         function PLASMA_sgetrs_incpiv_Tile_c(A,L,IPIV,B) &
          & bind(c, name='PLASMA_sgetrs_incpiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sgetrs_incpiv_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
         end function PLASMA_sgetrs_incpiv_Tile_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
!      private PLASMA_ssymm_Tile_c
      interface
         function PLASMA_ssymm_Tile_c(side,uplo,alpha,A,B,beta,C) &
          & bind(c, name='PLASMA_ssymm_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ssymm_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            real(kind=c_float), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            real(kind=c_float), value :: beta
            type(c_ptr), value :: C
         end function PLASMA_ssymm_Tile_c
      end interface

!      private PLASMA_ssyrk_Tile_c
      interface
         function PLASMA_ssyrk_Tile_c(uplo,trans,alpha,A,beta,C) &
          & bind(c, name='PLASMA_ssyrk_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ssyrk_Tile_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            real(kind=c_float), value :: alpha
            type(c_ptr), value :: A
            real(kind=c_float), value :: beta
            type(c_ptr), value :: C
         end function PLASMA_ssyrk_Tile_c
      end interface

!      private PLASMA_ssyr2k_Tile_c
      interface
         function PLASMA_ssyr2k_Tile_c(uplo,trans,alpha,A,B,beta,C) &
          & bind(c, name='PLASMA_ssyr2k_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ssyr2k_Tile_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            real(kind=c_float), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            real(kind=c_float), value :: beta
            type(c_ptr), value :: C
         end function PLASMA_ssyr2k_Tile_c
      end interface
#endif

!      private PLASMA_ssyev_Tile_c
      interface
         function PLASMA_ssyev_Tile_c(jobz,uplo,A,W,T,Q) &
          & bind(c, name='PLASMA_ssyev_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ssyev_Tile_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
         end function PLASMA_ssyev_Tile_c
      end interface

!      private PLASMA_ssygv_Tile_c
      interface
         function PLASMA_ssygv_Tile_c(itype,jobz,uplo,A,B,W,T,Q) &
          & bind(c, name='PLASMA_ssygv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ssygv_Tile_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
         end function PLASMA_ssygv_Tile_c
      end interface

!      private PLASMA_ssygst_Tile_c
      interface
         function PLASMA_ssygst_Tile_c(itype,uplo,A,B) &
          & bind(c, name='PLASMA_ssygst_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ssygst_Tile_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
         end function PLASMA_ssygst_Tile_c
      end interface

!      private PLASMA_ssytrd_Tile_c
      interface
         function PLASMA_ssytrd_Tile_c(jobz,uplo,A,D,E,T,Q) &
          & bind(c, name='PLASMA_ssytrd_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ssytrd_Tile_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: D
            type(c_ptr), value :: E
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
         end function PLASMA_ssytrd_Tile_c
      end interface

!      private PLASMA_slange_Tile_c
      interface
         function PLASMA_slange_Tile_c(norm,A,work) &
          & bind(c, name='PLASMA_slange_Tile')
            use iso_c_binding
            implicit none
            real(kind=c_float) :: PLASMA_slange_Tile_c
            integer(kind=c_int), value :: norm
            type(c_ptr), value :: A
            type(c_ptr), value :: work
         end function PLASMA_slange_Tile_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
!      private PLASMA_slansy_Tile_c
      interface
         function PLASMA_slansy_Tile_c(norm,uplo,A,work) &
          & bind(c, name='PLASMA_slansy_Tile')
            use iso_c_binding
            implicit none
            real(kind=c_float) :: PLASMA_slansy_Tile_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: work
         end function PLASMA_slansy_Tile_c
      end interface
#endif

!      private PLASMA_slansy_Tile_c
      interface
         function PLASMA_slansy_Tile_c(norm,uplo,A,work) &
          & bind(c, name='PLASMA_slansy_Tile')
            use iso_c_binding
            implicit none
            real(kind=c_float) :: PLASMA_slansy_Tile_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: work
         end function PLASMA_slansy_Tile_c
      end interface

!      private PLASMA_slaswp_Tile_c
      interface
         function PLASMA_slaswp_Tile_c(A,K1,K2,IPIV,INCX) &
          & bind(c, name='PLASMA_slaswp_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_slaswp_Tile_c
            type(c_ptr), value :: A
            integer(kind=c_int), value :: K1
            integer(kind=c_int), value :: K2
            type(c_ptr), value :: IPIV
            integer(kind=c_int), value :: INCX
         end function PLASMA_slaswp_Tile_c
      end interface

!      private PLASMA_slauum_Tile_c
      interface
         function PLASMA_slauum_Tile_c(uplo,A) &
          & bind(c, name='PLASMA_slauum_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_slauum_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
         end function PLASMA_slauum_Tile_c
      end interface

! PLASMA_splgsy_Tile uses unsigned ints, not generating Fortran interface

! PLASMA_splgsy_Tile uses unsigned ints, not generating Fortran interface

! PLASMA_splrnt_Tile uses unsigned ints, not generating Fortran interface

!      private PLASMA_sposv_Tile_c
      interface
         function PLASMA_sposv_Tile_c(uplo,A,B) &
          & bind(c, name='PLASMA_sposv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sposv_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
         end function PLASMA_sposv_Tile_c
      end interface

!      private PLASMA_spotrf_Tile_c
      interface
         function PLASMA_spotrf_Tile_c(uplo,A) &
          & bind(c, name='PLASMA_spotrf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_spotrf_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
         end function PLASMA_spotrf_Tile_c
      end interface

!      private PLASMA_spotri_Tile_c
      interface
         function PLASMA_spotri_Tile_c(uplo,A) &
          & bind(c, name='PLASMA_spotri_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_spotri_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
         end function PLASMA_spotri_Tile_c
      end interface

!      private PLASMA_spotrs_Tile_c
      interface
         function PLASMA_spotrs_Tile_c(uplo,A,B) &
          & bind(c, name='PLASMA_spotrs_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_spotrs_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
         end function PLASMA_spotrs_Tile_c
      end interface

!      private PLASMA_ssymm_Tile_c
      interface
         function PLASMA_ssymm_Tile_c(side,uplo,alpha,A,B,beta,C) &
          & bind(c, name='PLASMA_ssymm_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ssymm_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            real(kind=c_float), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            real(kind=c_float), value :: beta
            type(c_ptr), value :: C
         end function PLASMA_ssymm_Tile_c
      end interface

!      private PLASMA_ssyrk_Tile_c
      interface
         function PLASMA_ssyrk_Tile_c(uplo,trans,alpha,A,beta,C) &
          & bind(c, name='PLASMA_ssyrk_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ssyrk_Tile_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            real(kind=c_float), value :: alpha
            type(c_ptr), value :: A
            real(kind=c_float), value :: beta
            type(c_ptr), value :: C
         end function PLASMA_ssyrk_Tile_c
      end interface

!      private PLASMA_ssyr2k_Tile_c
      interface
         function PLASMA_ssyr2k_Tile_c(uplo,trans,alpha,A,B,beta,C) &
          & bind(c, name='PLASMA_ssyr2k_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ssyr2k_Tile_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            real(kind=c_float), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            real(kind=c_float), value :: beta
            type(c_ptr), value :: C
         end function PLASMA_ssyr2k_Tile_c
      end interface

!      private PLASMA_strmm_Tile_c
      interface
         function PLASMA_strmm_Tile_c(side,uplo,transA,diag,alpha,A,B) &
          & bind(c, name='PLASMA_strmm_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_strmm_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            real(kind=c_float), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
         end function PLASMA_strmm_Tile_c
      end interface

!      private PLASMA_strsm_Tile_c
      interface
         function PLASMA_strsm_Tile_c(side,uplo,transA,diag,alpha,A,B) &
          & bind(c, name='PLASMA_strsm_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_strsm_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            real(kind=c_float), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
         end function PLASMA_strsm_Tile_c
      end interface

!      private PLASMA_strsmpl_Tile_c
      interface
         function PLASMA_strsmpl_Tile_c(A,L,IPIV,B) &
          & bind(c, name='PLASMA_strsmpl_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_strsmpl_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
         end function PLASMA_strsmpl_Tile_c
      end interface

!      private PLASMA_strtri_Tile_c
      interface
         function PLASMA_strtri_Tile_c(uplo,diag,A) &
          & bind(c, name='PLASMA_strtri_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_strtri_Tile_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: diag
            type(c_ptr), value :: A
         end function PLASMA_strtri_Tile_c
      end interface

!      private PLASMA_sorglq_Tile_c
      interface
         function PLASMA_sorglq_Tile_c(A,T,B) &
          & bind(c, name='PLASMA_sorglq_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sorglq_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
         end function PLASMA_sorglq_Tile_c
      end interface

!      private PLASMA_sorgqr_Tile_c
      interface
         function PLASMA_sorgqr_Tile_c(A,T,B) &
          & bind(c, name='PLASMA_sorgqr_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sorgqr_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
         end function PLASMA_sorgqr_Tile_c
      end interface

!      private PLASMA_sormlq_Tile_c
      interface
         function PLASMA_sormlq_Tile_c(side,trans,A,T,B) &
          & bind(c, name='PLASMA_sormlq_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sormlq_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
         end function PLASMA_sormlq_Tile_c
      end interface

!      private PLASMA_sormqr_Tile_c
      interface
         function PLASMA_sormqr_Tile_c(side,trans,A,T,B) &
          & bind(c, name='PLASMA_sormqr_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sormqr_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
         end function PLASMA_sormqr_Tile_c
      end interface

!      private PLASMA_Alloc_Workspace_sgelqf_c
      interface
         function PLASMA_Alloc_Workspace_sgelqf_c(M,N,T) &
          & bind(c, name='PLASMA_Alloc_Workspace_sgelqf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_sgelqf_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: T ! T is **, so pass by reference
         end function PLASMA_Alloc_Workspace_sgelqf_c
      end interface

!      private PLASMA_Alloc_Workspace_sgels_c
      interface
         function PLASMA_Alloc_Workspace_sgels_c(M,N,T) &
          & bind(c, name='PLASMA_Alloc_Workspace_sgels')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_sgels_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: T ! T is **, so pass by reference
         end function PLASMA_Alloc_Workspace_sgels_c
      end interface

!      private PLASMA_Alloc_Workspace_sgeqrf_c
      interface
         function PLASMA_Alloc_Workspace_sgeqrf_c(M,N,T) &
          & bind(c, name='PLASMA_Alloc_Workspace_sgeqrf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_sgeqrf_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: T ! T is **, so pass by reference
         end function PLASMA_Alloc_Workspace_sgeqrf_c
      end interface

!      private PLASMA_Alloc_Workspace_sgesv_incpiv_c
      interface
         function PLASMA_Alloc_Workspace_sgesv_incpiv_c(N,L,IPIV) &
          & bind(c, name='PLASMA_Alloc_Workspace_sgesv_incpiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_sgesv_incpiv_c
            integer(kind=c_int), value :: N
            type(c_ptr) :: L ! L is **, so pass by reference
            type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         end function PLASMA_Alloc_Workspace_sgesv_incpiv_c
      end interface

!      private PLASMA_Alloc_Workspace_sgetrf_incpiv_c
      interface
         function PLASMA_Alloc_Workspace_sgetrf_incpiv_c(M,N,L,IPIV) &
          & bind(c, name='PLASMA_Alloc_Workspace_sgetrf_incpiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_sgetrf_incpiv_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: L ! L is **, so pass by reference
            type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         end function PLASMA_Alloc_Workspace_sgetrf_incpiv_c
      end interface

!      private PLASMA_Alloc_Workspace_sgeev_c
      interface
         function PLASMA_Alloc_Workspace_sgeev_c(N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_sgeev')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_sgeev_c
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_sgeev_c
      end interface

!      private PLASMA_Alloc_Workspace_sgebrd_c
      interface
         function PLASMA_Alloc_Workspace_sgebrd_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_sgebrd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_sgebrd_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_sgebrd_c
      end interface

!      private PLASMA_Alloc_Workspace_sgesvd_c
      interface
         function PLASMA_Alloc_Workspace_sgesvd_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_sgesvd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_sgesvd_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_sgesvd_c
      end interface

!      private PLASMA_Alloc_Workspace_ssyev_c
      interface
         function PLASMA_Alloc_Workspace_ssyev_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_ssyev')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_ssyev_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_ssyev_c
      end interface

!      private PLASMA_Alloc_Workspace_ssygv_c
      interface
         function PLASMA_Alloc_Workspace_ssygv_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_ssygv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_ssygv_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_ssygv_c
      end interface

!      private PLASMA_Alloc_Workspace_ssytrd_c
      interface
         function PLASMA_Alloc_Workspace_ssytrd_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_ssytrd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_ssytrd_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_ssytrd_c
      end interface

!      private PLASMA_Alloc_Workspace_sgelqf_Tile_c
      interface
         function PLASMA_Alloc_Workspace_sgelqf_Tile_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_sgelqf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_sgelqf_Tile_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_sgelqf_Tile_c
      end interface

!      private PLASMA_Alloc_Workspace_sgels_Tile_c
      interface
         function PLASMA_Alloc_Workspace_sgels_Tile_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_sgels_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_sgels_Tile_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_sgels_Tile_c
      end interface

!      private PLASMA_Alloc_Workspace_sgeqrf_Tile_c
      interface
         function PLASMA_Alloc_Workspace_sgeqrf_Tile_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_sgeqrf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_sgeqrf_Tile_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_sgeqrf_Tile_c
      end interface

!      private PLASMA_Alloc_Workspace_sgesv_incpiv_Tile_c
      interface
         function PLASMA_Alloc_Workspace_sgesv_incpiv_Tile_c(N,descL,IPIV) &
          & bind(c, name='PLASMA_Alloc_Workspace_sgesv_incpiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_sgesv_incpiv_Tile_c
            integer(kind=c_int), value :: N
            type(c_ptr) :: descL ! descL is **, so pass by reference
            type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         end function PLASMA_Alloc_Workspace_sgesv_incpiv_Tile_c
      end interface

!      private PLASMA_Alloc_Workspace_sgetrf_incpiv_Tile_c
      interface
         function PLASMA_Alloc_Workspace_sgetrf_incpiv_Tile_c(N,descL,IPIV) &
          & bind(c, name='PLASMA_Alloc_Workspace_sgetrf_incpiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_sgetrf_incpiv_Tile_c
            integer(kind=c_int), value :: N
            type(c_ptr) :: descL ! descL is **, so pass by reference
            type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         end function PLASMA_Alloc_Workspace_sgetrf_incpiv_Tile_c
      end interface

!      private PLASMA_sLapack_to_Tile_c
      interface
         function PLASMA_sLapack_to_Tile_c(Af77,LDA,A) &
          & bind(c, name='PLASMA_sLapack_to_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sLapack_to_Tile_c
            type(c_ptr), value :: Af77
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: A
         end function PLASMA_sLapack_to_Tile_c
      end interface

!      private PLASMA_sTile_to_Lapack_c
      interface
         function PLASMA_sTile_to_Lapack_c(A,Af77,LDA) &
          & bind(c, name='PLASMA_sTile_to_Lapack')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_sTile_to_Lapack_c
            type(c_ptr), value :: A
            type(c_ptr), value :: Af77
            integer(kind=c_int), value :: LDA
         end function PLASMA_sTile_to_Lapack_c
      end interface

  contains

      subroutine PLASMA_sgebrd(jobu,jobvt,M,N,A,LDA,D,E,U,LDU,VT,LDVT,T,info)
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
         real(kind=c_float), intent(inout), target :: A(LDA,*)
         real(kind=c_float), intent(out), target :: D(*)
         real(kind=c_float), intent(out), target :: E(*)
         real(kind=c_float), intent(out), target :: U(LDU,*)
         real(kind=c_float), intent(out), target :: VT(LDVT,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sgebrd_c(jobu,jobvt,M,N,c_loc(A),LDA,c_loc(D),c_loc(E),c_loc(U),LDU,c_loc(VT),LDVT,T)
      end subroutine PLASMA_sgebrd

      subroutine PLASMA_sgelqf(M,N,A,LDA,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         real(kind=c_float), intent(inout), target :: A(LDA,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sgelqf_c(M,N,c_loc(A),LDA,T)
      end subroutine PLASMA_sgelqf

      subroutine PLASMA_sgelqs(M,N,NRHS,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         real(kind=c_float), intent(in), target :: A(LDA,*)
         real(kind=c_float), intent(inout), target :: B(LDB,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sgelqs_c(M,N,NRHS,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_sgelqs

      subroutine PLASMA_sgels(trans,M,N,NRHS,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: trans
         real(kind=c_float), intent(inout), target :: A(LDA,*)
         real(kind=c_float), intent(inout), target :: B(LDB,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sgels_c(trans,M,N,NRHS,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_sgels

      subroutine PLASMA_sgemm(transA,transB,M,N,K,alpha,A,LDA,B,LDB,beta,C,LDC,info)
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
         real(kind=c_float), intent(in) :: alpha
         real(kind=c_float), intent(in) :: beta
         real(kind=c_float), intent(in), target :: A(LDA,*)
         real(kind=c_float), intent(in), target :: B(LDB,*)
         real(kind=c_float), intent(inout), target :: C(LDC,*)
         info = PLASMA_sgemm_c(transA,transB,M,N,K,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
      end subroutine PLASMA_sgemm

      subroutine PLASMA_sgeqrf(M,N,A,LDA,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         real(kind=c_float), intent(inout), target :: A(LDA,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sgeqrf_c(M,N,c_loc(A),LDA,T)
      end subroutine PLASMA_sgeqrf

      subroutine PLASMA_sgeqrs(M,N,NRHS,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         real(kind=c_float), intent(inout), target :: A(LDA,*)
         real(kind=c_float), intent(inout), target :: B(LDB,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sgeqrs_c(M,N,NRHS,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_sgeqrs

      subroutine PLASMA_sgesv(N,NRHS,A,LDA,IPIV,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(out), target :: IPIV(*)
         real(kind=c_float), intent(inout), target :: A(LDA,*)
         real(kind=c_float), intent(inout), target :: B(LDB,*)
         info = PLASMA_sgesv_c(N,NRHS,c_loc(A),LDA,c_loc(IPIV),c_loc(B),LDB)
      end subroutine PLASMA_sgesv

      subroutine PLASMA_sgesv_incpiv(N,NRHS,A,LDA,L,IPIV,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         real(kind=c_float), intent(inout), target :: A(LDA,*)
         real(kind=c_float), intent(inout), target :: B(LDB,*)
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sgesv_incpiv_c(N,NRHS,c_loc(A),LDA,L,IPIV,c_loc(B),LDB)
      end subroutine PLASMA_sgesv_incpiv

      subroutine PLASMA_sgesvd(jobu,jobvt,M,N,A,LDA,S,U,LDU,VT,LDVT,T,info)
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
         real(kind=c_float), intent(inout), target :: A(LDA,*)
         real(kind=c_float), intent(out), target :: S(*)
         real(kind=c_float), intent(out), target :: U(LDU,*)
         real(kind=c_float), intent(out), target :: VT(LDVT,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sgesvd_c(jobu,jobvt,M,N,c_loc(A),LDA,c_loc(S),c_loc(U),LDU,c_loc(VT),LDVT,T)
      end subroutine PLASMA_sgesvd

      subroutine PLASMA_sgetrf(M,N,A,LDA,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(out), target :: IPIV(*)
         real(kind=c_float), intent(inout), target :: A(LDA,*)
         info = PLASMA_sgetrf_c(M,N,c_loc(A),LDA,c_loc(IPIV))
      end subroutine PLASMA_sgetrf

      subroutine PLASMA_sgetrf_incpiv(M,N,A,LDA,L,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         real(kind=c_float), intent(inout), target :: A(LDA,*)
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L    ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sgetrf_incpiv_c(M,N,c_loc(A),LDA,L,IPIV)
      end subroutine PLASMA_sgetrf_incpiv

      subroutine PLASMA_sgetrs(trans,N,NRHS,A,LDA,IPIV,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in), target :: IPIV(*)
         real(kind=c_float), intent(in), target :: A(LDA,*)
         real(kind=c_float), intent(inout), target :: B(LDB,*)
         info = PLASMA_sgetrs_c(trans,N,NRHS,c_loc(A),LDA,c_loc(IPIV),c_loc(B),LDB)
      end subroutine PLASMA_sgetrs

      subroutine PLASMA_sgetrs_incpiv(trans,N,NRHS,A,LDA,L,IPIV,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: trans
         real(kind=c_float), intent(in), target :: A(LDA,*)
         real(kind=c_float), intent(inout), target :: B(LDB,*)
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L    ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sgetrs_incpiv_c(trans,N,NRHS,c_loc(A),LDA,L,IPIV,c_loc(B),LDB)
      end subroutine PLASMA_sgetrs_incpiv

#if defined(PRECISION_z) || defined(PRECISION_c)
      subroutine PLASMA_ssymm(side,uplo,M,N,alpha,A,LDA,B,LDB,beta,C,LDC,info)
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
         real(kind=c_float), intent(in) :: alpha
         real(kind=c_float), intent(in) :: beta
         real(kind=c_float), intent(in), target :: A(LDA,*)
         real(kind=c_float), intent(in), target :: B(LDB,*)
         real(kind=c_float), intent(inout), target :: C(LDC,*)
         info = PLASMA_ssymm_c(side,uplo,M,N,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
      end subroutine PLASMA_ssymm

      subroutine PLASMA_ssyrk(uplo,trans,N,K,alpha,A,LDA,beta,C,LDC,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDC
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_float), intent(in), target :: A(LDA,*)
         real(kind=c_float), intent(inout), target :: C(LDC,*)
         real(kind=c_float), intent(in) :: alpha
         real(kind=c_float), intent(in) :: beta
         info = PLASMA_ssyrk_c(uplo,trans,N,K,alpha,c_loc(A),LDA,beta,c_loc(C),LDC)
      end subroutine PLASMA_ssyrk

      subroutine PLASMA_ssyr2k(uplo,trans,N,K,alpha,A,LDA,B,LDB,beta,C,LDC,info)
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
         real(kind=c_float), intent(in) :: alpha
         real(kind=c_float), intent(in), target :: A(LDA,*)
         real(kind=c_float), intent(in), target :: B(LDB,*)
         real(kind=c_float), intent(inout), target :: C(LDC,*)
         real(kind=c_float), intent(in) :: beta
         info = PLASMA_ssyr2k_c(uplo,trans,N,K,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
      end subroutine PLASMA_ssyr2k
#endif

      subroutine PLASMA_ssyev(jobz,uplo,N,A,LDA,W,T,Q,LDQ,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDQ
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_float), intent(inout), target :: A(LDA,*)
         real(kind=c_float), intent(out), target :: W(*)
         real(kind=c_float), intent(out), target :: Q(LDQ,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ssyev_c(jobz,uplo,N,c_loc(A),LDA,c_loc(W),T,c_loc(Q),LDQ)
      end subroutine PLASMA_ssyev

      subroutine PLASMA_ssygv(itype,jobz,uplo,N,A,LDA,B,LDB,W,T,Q,LDQ,info)
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
         real(kind=c_float), intent(inout), target :: A(LDA,*)
         real(kind=c_float), intent(inout), target :: B(LDB,*)
         real(kind=c_float), intent(out), target :: W(*)
         real(kind=c_float), intent(out), target :: Q(LDQ,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ssygv_c(itype,jobz,uplo,N,c_loc(A),LDA,c_loc(B),LDB,c_loc(W),T,c_loc(Q),LDQ)
      end subroutine PLASMA_ssygv

      subroutine PLASMA_ssygst(itype,uplo,N,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_float), intent(in), target :: B(LDB,*)
         real(kind=c_float), intent(inout), target :: A(LDA,*)
         info = PLASMA_ssygst_c(itype,uplo,N,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine PLASMA_ssygst

      subroutine PLASMA_ssytrd(jobz,uplo,N,A,LDA,D,E,T,Q,LDQ,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDQ
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_float), intent(inout), target :: A(LDA,*)
         real(kind=c_float), intent(out), target :: D(*)
         real(kind=c_float), intent(out), target :: E(*)
         real(kind=c_float), intent(out), target :: Q(LDQ,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ssytrd_c(jobz,uplo,N,c_loc(A),LDA,c_loc(D),c_loc(E),T,c_loc(Q),LDQ)
      end subroutine PLASMA_ssytrd

      function PLASMA_slange(norm,M,N,A,LDA,work)
         use iso_c_binding
         implicit none
         real(kind=c_float) :: PLASMA_slange
         real(kind=c_float), intent(inout), target :: work(*)
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         real(kind=c_float), intent(in), target :: A(LDA,*)
         PLASMA_slange = PLASMA_slange_c(norm,M,N,c_loc(A),LDA,c_loc(work))
      end function PLASMA_slange

#if defined(PRECISION_z) || defined(PRECISION_c)
      function PLASMA_slansy(norm,uplo,N,A,LDA,work)
         use iso_c_binding
         implicit none
         real(kind=c_float) :: PLASMA_slansy
         real(kind=c_float), intent(inout), target :: work(*)
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         real(kind=c_float), intent(in), target :: A(LDA,*)
         PLASMA_slansy = PLASMA_slansy_c(norm,uplo,N,c_loc(A),LDA,c_loc(work))
      end function PLASMA_slansy
#endif

      function PLASMA_slansy(norm,uplo,N,A,LDA,work)
         use iso_c_binding
         implicit none
         real(kind=c_float) :: PLASMA_slansy
         real(kind=c_float), intent(inout), target :: work(*)
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         real(kind=c_float), intent(inout), target :: A(LDA,*)
         PLASMA_slansy = PLASMA_slansy_c(norm,uplo,N,c_loc(A),LDA,c_loc(work))
      end function PLASMA_slansy

      subroutine PLASMA_slaswp(N,A,LDA,K1,K2,IPIV,INCX,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: INCX
         integer(kind=c_int), intent(in) :: K1
         integer(kind=c_int), intent(in) :: K2
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in), target :: IPIV(*)
         real(kind=c_float), intent(inout), target :: A(LDA,*)
         info = PLASMA_slaswp_c(N,c_loc(A),LDA,K1,K2,c_loc(IPIV),INCX)
      end subroutine PLASMA_slaswp

      subroutine PLASMA_slauum(uplo,N,A,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_float), intent(inout), target :: A(LDA,*)
         info = PLASMA_slauum_c(uplo,N,c_loc(A),LDA)
      end subroutine PLASMA_slauum

      subroutine PLASMA_sposv(uplo,N,NRHS,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_float), intent(inout), target :: A(LDA,*)
         real(kind=c_float), intent(inout), target :: B(LDB,*)
         info = PLASMA_sposv_c(uplo,N,NRHS,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine PLASMA_sposv

      subroutine PLASMA_spotrf(uplo,N,A,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_float), intent(inout), target :: A(LDA,*)
         info = PLASMA_spotrf_c(uplo,N,c_loc(A),LDA)
      end subroutine PLASMA_spotrf

      subroutine PLASMA_spotri(uplo,N,A,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_float), intent(inout), target :: A(LDA,*)
         info = PLASMA_spotri_c(uplo,N,c_loc(A),LDA)
      end subroutine PLASMA_spotri

      subroutine PLASMA_spotrs(uplo,N,NRHS,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_float), intent(in), target :: A(LDA,*)
         real(kind=c_float), intent(inout), target :: B(LDB,*)
         info = PLASMA_spotrs_c(uplo,N,NRHS,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine PLASMA_spotrs

      subroutine PLASMA_ssymm(side,uplo,M,N,alpha,A,LDA,B,LDB,beta,C,LDC,info)
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
         real(kind=c_float), intent(in) :: alpha
         real(kind=c_float), intent(in) :: beta
         real(kind=c_float), intent(in), target :: A(LDA,*)
         real(kind=c_float), intent(in), target :: B(LDB,*)
         real(kind=c_float), intent(inout), target :: C(LDC,*)
         info = PLASMA_ssymm_c(side,uplo,M,N,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
      end subroutine PLASMA_ssymm

      subroutine PLASMA_ssyrk(uplo,trans,N,K,alpha,A,LDA,beta,C,LDC,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDC
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_float), intent(in) :: alpha
         real(kind=c_float), intent(in) :: beta
         real(kind=c_float), intent(in), target :: A(LDA,*)
         real(kind=c_float), intent(inout), target :: C(LDC,*)
         info = PLASMA_ssyrk_c(uplo,trans,N,K,alpha,c_loc(A),LDA,beta,c_loc(C),LDC)
      end subroutine PLASMA_ssyrk

      subroutine PLASMA_ssyr2k(uplo,trans,N,K,alpha,A,LDA,B,LDB,beta,C,LDC,info)
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
         real(kind=c_float), intent(in) :: alpha
         real(kind=c_float), intent(in) :: beta
         real(kind=c_float), intent(in), target :: A(LDA,*)
         real(kind=c_float), intent(in), target :: B(LDB,*)
         real(kind=c_float), intent(inout), target :: C(LDC,*)
         info = PLASMA_ssyr2k_c(uplo,trans,N,K,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
      end subroutine PLASMA_ssyr2k

      subroutine PLASMA_strmm(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB,info)
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
         real(kind=c_float), intent(in) :: alpha
         real(kind=c_float), intent(in), target :: A(LDA,*)
         real(kind=c_float), intent(inout), target :: B(LDB,*)
         info = PLASMA_strmm_c(side,uplo,transA,diag,N,NRHS,alpha,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine PLASMA_strmm

      subroutine PLASMA_strsm(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB,info)
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
         real(kind=c_float), intent(in) :: alpha
         real(kind=c_float), intent(in), target :: A(LDA,*)
         real(kind=c_float), intent(inout), target :: B(LDB,*)
         info = PLASMA_strsm_c(side,uplo,transA,diag,N,NRHS,alpha,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine PLASMA_strsm

      subroutine PLASMA_strsmpl(N,NRHS,A,LDA,L,IPIV,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         real(kind=c_float), intent(in), target :: A(LDA,*)
         type(c_ptr), value :: L    ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         real(kind=c_float), intent(inout), target :: B(LDB,*)
         info = PLASMA_strsmpl_c(N,NRHS,c_loc(A),LDA,L,IPIV,c_loc(B),LDB)
      end subroutine PLASMA_strsmpl

      subroutine PLASMA_strtri(uplo,diag,N,A,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_float), intent(inout), target :: A(LDA,*)
         info = PLASMA_strtri_c(uplo,diag,N,c_loc(A),LDA)
      end subroutine PLASMA_strtri

      subroutine PLASMA_sorglq(M,N,K,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         real(kind=c_float), intent(inout), target :: A(LDA,*)
         real(kind=c_float), intent(out), target :: B(LDB,*)
         info = PLASMA_sorglq_c(M,N,K,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_sorglq

      subroutine PLASMA_sorgqr(M,N,K,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         real(kind=c_float), intent(inout), target :: A(LDA,*)
         real(kind=c_float), intent(out), target :: B(LDB,*)
         info = PLASMA_sorgqr_c(M,N,K,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_sorgqr

      subroutine PLASMA_sormlq(side,trans,M,N,K,A,LDA,T,B,LDB,info)
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
         real(kind=c_float), intent(in), target :: A(LDA,*)
         real(kind=c_float), intent(inout), target :: B(LDB,*)
         info = PLASMA_sormlq_c(side,trans,M,N,K,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_sormlq

      subroutine PLASMA_sormqr(side,trans,M,N,K,A,LDA,T,B,LDB,info)
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
         real(kind=c_float), intent(in), target :: A(LDA,*)
         real(kind=c_float), intent(inout), target :: B(LDB,*)
         info = PLASMA_sormqr_c(side,trans,M,N,K,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_sormqr

      subroutine PLASMA_sgecfi(m,n,A,fin,imb,inb,fout,omb,onb,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         real(kind=c_float), intent(inout), target :: A(*)
         integer(kind=c_int), intent(in) :: fin
         integer(kind=c_int), intent(in) :: fout
         integer(kind=c_int), intent(in) :: imb
         integer(kind=c_int), intent(in) :: inb
         integer(kind=c_int), intent(in) :: omb
         integer(kind=c_int), intent(in) :: onb
         integer(kind=c_int), intent(in) :: m
         integer(kind=c_int), intent(in) :: n
         info = PLASMA_sgecfi_c(m,n,c_loc(A),fin,imb,inb,fout,omb,onb)
      end subroutine PLASMA_sgecfi

      subroutine PLASMA_sgetmi(m,n,A,fin,mb,nb,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         real(kind=c_float), intent(inout), target :: A(*)
         integer(kind=c_int), intent(in) :: fin
         integer(kind=c_int), intent(in) :: mb
         integer(kind=c_int), intent(in) :: nb
         integer(kind=c_int), intent(in) :: m
         integer(kind=c_int), intent(in) :: n
         info = PLASMA_sgetmi_c(m,n,c_loc(A),fin,mb,nb)
      end subroutine PLASMA_sgetmi

      subroutine PLASMA_sgebrd_Tile(jobu,jobvt,A,D,E,U,VT,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobu
         integer(kind=c_int), intent(in) :: jobvt
         real(kind=c_float), intent(out), target :: D(*)
         real(kind=c_float), intent(out), target :: E(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: U ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: VT ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sgebrd_Tile_c(jobu,jobvt,A,c_loc(D),c_loc(E),U,VT,T)
      end subroutine PLASMA_sgebrd_Tile

      subroutine PLASMA_sgelqf_Tile(A,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sgelqf_Tile_c(A,T)
      end subroutine PLASMA_sgelqf_Tile

      subroutine PLASMA_sgelqs_Tile(A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sgelqs_Tile_c(A,T,B)
      end subroutine PLASMA_sgelqs_Tile

      subroutine PLASMA_sgels_Tile(trans,A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sgels_Tile_c(trans,A,T,B)
      end subroutine PLASMA_sgels_Tile

      subroutine PLASMA_sgemm_Tile(transA,transB,alpha,A,B,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: transB
         real(kind=c_float), intent(in) :: alpha
         real(kind=c_float), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sgemm_Tile_c(transA,transB,alpha,A,B,beta,C)
      end subroutine PLASMA_sgemm_Tile

      subroutine PLASMA_sgeqrf_Tile(A,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sgeqrf_Tile_c(A,T)
      end subroutine PLASMA_sgeqrf_Tile

      subroutine PLASMA_sgeqrs_Tile(A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sgeqrs_Tile_c(A,T,B)
      end subroutine PLASMA_sgeqrs_Tile

      subroutine PLASMA_sgesv_Tile(A,IPIV,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: IPIV(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sgesv_Tile_c(A,c_loc(IPIV),B)
      end subroutine PLASMA_sgesv_Tile

      subroutine PLASMA_sgesv_incpiv_Tile(A,L,IPIV,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sgesv_incpiv_Tile_c(A,L,IPIV,B)
      end subroutine PLASMA_sgesv_incpiv_Tile

      subroutine PLASMA_sgesvd_Tile(jobu,jobvt,A,S,U,VT,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobu
         integer(kind=c_int), intent(in) :: jobvt
         real(kind=c_float), intent(out), target :: S(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: U ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: VT ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sgesvd_Tile_c(jobu,jobvt,A,c_loc(S),U,VT,T)
      end subroutine PLASMA_sgesvd_Tile

      subroutine PLASMA_sgetrf_Tile(A,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: IPIV(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sgetrf_Tile_c(A,c_loc(IPIV))
      end subroutine PLASMA_sgetrf_Tile

      subroutine PLASMA_sgetrf_incpiv_Tile(A,L,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sgetrf_incpiv_Tile_c(A,L,IPIV)
      end subroutine PLASMA_sgetrf_incpiv_Tile

      subroutine PLASMA_sgetrs_Tile(trans,A,IPIV,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         integer(kind=c_int), intent(in), target :: IPIV(*)
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sgetrs_Tile_c(trans,A,c_loc(IPIV),B)
      end subroutine PLASMA_sgetrs_Tile

      subroutine PLASMA_sgetrs_incpiv_Tile(A,L,IPIV,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sgetrs_incpiv_Tile_c(A,L,IPIV,B)
      end subroutine PLASMA_sgetrs_incpiv_Tile

#if defined(PRECISION_z) || defined(PRECISION_c)
      subroutine PLASMA_ssymm_Tile(side,uplo,alpha,A,B,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_float), intent(in) :: alpha
         real(kind=c_float), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ssymm_Tile_c(side,uplo,alpha,A,B,beta,C)
      end subroutine PLASMA_ssymm_Tile

      subroutine PLASMA_ssyrk_Tile(uplo,trans,alpha,A,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_float), intent(in) :: alpha
         real(kind=c_float), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ssyrk_Tile_c(uplo,trans,alpha,A,beta,C)
      end subroutine PLASMA_ssyrk_Tile

      subroutine PLASMA_ssyr2k_Tile(uplo,trans,alpha,A,B,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_float), intent(in) :: alpha
         real(kind=c_float), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ssyr2k_Tile_c(uplo,trans,alpha,A,B,beta,C)
      end subroutine PLASMA_ssyr2k_Tile
#endif

      subroutine PLASMA_ssyev_Tile(jobz,uplo,A,W,T,Q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_float), intent(out), target :: W(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: Q ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ssyev_Tile_c(jobz,uplo,A,c_loc(W),T,Q)
      end subroutine PLASMA_ssyev_Tile

      subroutine PLASMA_ssygv_Tile(itype,jobz,uplo,A,B,W,T,Q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_float), intent(out), target :: W(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: Q ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ssygv_Tile_c(itype,jobz,uplo,A,B,c_loc(W),T,Q)
      end subroutine PLASMA_ssygv_Tile

      subroutine PLASMA_ssygst_Tile(itype,uplo,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ssygst_Tile_c(itype,uplo,A,B)
      end subroutine PLASMA_ssygst_Tile

      subroutine PLASMA_ssytrd_Tile(jobz,uplo,A,D,E,T,Q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_float), intent(out), target :: D(*)
         real(kind=c_float), intent(out), target :: E(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: Q ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ssytrd_Tile_c(jobz,uplo,A,c_loc(D),c_loc(E),T,Q)
      end subroutine PLASMA_ssytrd_Tile

      function PLASMA_slange_Tile(norm,A,work)
         use iso_c_binding
         implicit none
         real(kind=c_float) :: PLASMA_slange_Tile
         integer(kind=c_int), intent(in) :: norm
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         real(kind=c_float), intent(inout), target :: work(*)
         PLASMA_slange_Tile = PLASMA_slange_Tile_c(norm,A,c_loc(work))
       end function PLASMA_slange_Tile

#if defined(PRECISION_z) || defined(PRECISION_c)
      function PLASMA_slansy_Tile(norm,uplo,A,work)
         use iso_c_binding
         implicit none
         real(kind=c_float) :: PLASMA_slansy_Tile
         real(kind=c_float), intent(inout), target :: work(*)
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         PLASMA_slansy_Tile = PLASMA_slansy_Tile_c(norm,uplo,A,c_loc(work))
      end function PLASMA_slansy_Tile
#endif

      function PLASMA_slansy_Tile(norm,uplo,A,work)
         use iso_c_binding
         implicit none
         real(kind=c_float) :: PLASMA_slansy_Tile
         real(kind=c_float), intent(inout), target :: work(*)
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         PLASMA_slansy_Tile = PLASMA_slansy_Tile_c(norm,uplo,A,c_loc(work))
      end function PLASMA_slansy_Tile

      subroutine PLASMA_slaswp_Tile(A,K1,K2,IPIV,INCX,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: INCX
         integer(kind=c_int), intent(in) :: K1
         integer(kind=c_int), intent(in) :: K2
         integer(kind=c_int), intent(in), target :: IPIV(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_slaswp_Tile_c(A,K1,K2,c_loc(IPIV),INCX)
      end subroutine PLASMA_slaswp_Tile

      subroutine PLASMA_slauum_Tile(uplo,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_slauum_Tile_c(uplo,A)
      end subroutine PLASMA_slauum_Tile

      subroutine PLASMA_sposv_Tile(uplo,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sposv_Tile_c(uplo,A,B)
      end subroutine PLASMA_sposv_Tile

      subroutine PLASMA_spotrf_Tile(uplo,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_spotrf_Tile_c(uplo,A)
      end subroutine PLASMA_spotrf_Tile

      subroutine PLASMA_spotri_Tile(uplo,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_spotri_Tile_c(uplo,A)
      end subroutine PLASMA_spotri_Tile

      subroutine PLASMA_spotrs_Tile(uplo,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_spotrs_Tile_c(uplo,A,B)
      end subroutine PLASMA_spotrs_Tile

      subroutine PLASMA_ssymm_Tile(side,uplo,alpha,A,B,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_float), intent(in) :: alpha
         real(kind=c_float), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ssymm_Tile_c(side,uplo,alpha,A,B,beta,C)
      end subroutine PLASMA_ssymm_Tile

      subroutine PLASMA_ssyrk_Tile(uplo,trans,alpha,A,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_float), intent(in) :: alpha
         real(kind=c_float), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ssyrk_Tile_c(uplo,trans,alpha,A,beta,C)
      end subroutine PLASMA_ssyrk_Tile

      subroutine PLASMA_ssyr2k_Tile(uplo,trans,alpha,A,B,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_float), intent(in) :: alpha
         real(kind=c_float), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ssyr2k_Tile_c(uplo,trans,alpha,A,B,beta,C)
      end subroutine PLASMA_ssyr2k_Tile

      subroutine PLASMA_strmm_Tile(side,uplo,transA,diag,alpha,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_float), intent(in) :: alpha
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_strmm_Tile_c(side,uplo,transA,diag,alpha,A,B)
      end subroutine PLASMA_strmm_Tile

      subroutine PLASMA_strsm_Tile(side,uplo,transA,diag,alpha,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_float), intent(in) :: alpha
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_strsm_Tile_c(side,uplo,transA,diag,alpha,A,B)
      end subroutine PLASMA_strsm_Tile

      subroutine PLASMA_strsmpl_Tile(A,L,IPIV,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_strsmpl_Tile_c(A,L,IPIV,B)
      end subroutine PLASMA_strsmpl_Tile

      subroutine PLASMA_strtri_Tile(uplo,diag,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_strtri_Tile_c(uplo,diag,A)
      end subroutine PLASMA_strtri_Tile

      subroutine PLASMA_sorglq_Tile(A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sorglq_Tile_c(A,T,B)
      end subroutine PLASMA_sorglq_Tile

      subroutine PLASMA_sorgqr_Tile(A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sorgqr_Tile_c(A,T,B)
      end subroutine PLASMA_sorgqr_Tile

      subroutine PLASMA_sormlq_Tile(side,trans,A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sormlq_Tile_c(side,trans,A,T,B)
      end subroutine PLASMA_sormlq_Tile

      subroutine PLASMA_sormqr_Tile(side,trans,A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sormqr_Tile_c(side,trans,A,T,B)
      end subroutine PLASMA_sormqr_Tile

      subroutine PLASMA_Alloc_Workspace_sgelqf(M,N,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         info = PLASMA_Alloc_Workspace_sgelqf_c(M,N,T)
      end subroutine PLASMA_Alloc_Workspace_sgelqf

      subroutine PLASMA_Alloc_Workspace_sgels(M,N,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         info = PLASMA_Alloc_Workspace_sgels_c(M,N,T)
      end subroutine PLASMA_Alloc_Workspace_sgels

      subroutine PLASMA_Alloc_Workspace_sgeqrf(M,N,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         info = PLASMA_Alloc_Workspace_sgeqrf_c(M,N,T)
      end subroutine PLASMA_Alloc_Workspace_sgeqrf

      subroutine PLASMA_Alloc_Workspace_sgesv_incpiv(N,L,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         type(c_ptr) :: L ! L is **, so pass by reference
         info = PLASMA_Alloc_Workspace_sgesv_incpiv_c(N,L,IPIV)
      end subroutine PLASMA_Alloc_Workspace_sgesv_incpiv

      subroutine PLASMA_Alloc_Workspace_sgetrf_incpiv(M,N,L,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         type(c_ptr) :: L ! L is **, so pass by reference
         info = PLASMA_Alloc_Workspace_sgetrf_incpiv_c(M,N,L,IPIV)
      end subroutine PLASMA_Alloc_Workspace_sgetrf_incpiv

      subroutine PLASMA_Alloc_Workspace_sgeev(N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_sgeev_c(N,descT)
      end subroutine PLASMA_Alloc_Workspace_sgeev

      subroutine PLASMA_Alloc_Workspace_sgebrd(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_sgebrd_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_sgebrd

      subroutine PLASMA_Alloc_Workspace_sgesvd(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_sgesvd_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_sgesvd

      subroutine PLASMA_Alloc_Workspace_ssyev(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_ssyev_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_ssyev

      subroutine PLASMA_Alloc_Workspace_ssygv(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_ssygv_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_ssygv

      subroutine PLASMA_Alloc_Workspace_ssytrd(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_ssytrd_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_ssytrd

      subroutine PLASMA_Alloc_Workspace_sgelqf_Tile(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_sgelqf_Tile_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_sgelqf_Tile

      subroutine PLASMA_Alloc_Workspace_sgels_Tile(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_sgels_Tile_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_sgels_Tile

      subroutine PLASMA_Alloc_Workspace_sgeqrf_Tile(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_sgeqrf_Tile_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_sgeqrf_Tile

      subroutine PLASMA_Alloc_Workspace_sgesv_incpiv_Tile(N,descL,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         type(c_ptr) :: descL ! descL is **, so pass by reference
         info = PLASMA_Alloc_Workspace_sgesv_incpiv_Tile_c(N,descL,IPIV)
      end subroutine PLASMA_Alloc_Workspace_sgesv_incpiv_Tile

      subroutine PLASMA_Alloc_Workspace_sgetrf_incpiv_Tile(N,descL,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         type(c_ptr) :: descL ! descL is **, so pass by reference
         info = PLASMA_Alloc_Workspace_sgetrf_incpiv_Tile_c(N,descL,IPIV)
      end subroutine PLASMA_Alloc_Workspace_sgetrf_incpiv_Tile

      subroutine PLASMA_sLapack_to_Tile(Af77,LDA,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         real(kind=c_float), intent(in), target :: Af77(LDA,*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sLapack_to_Tile_c(c_loc(Af77),LDA,A)
      end subroutine PLASMA_sLapack_to_Tile

      subroutine PLASMA_sTile_to_Lapack(A,Af77,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         real(kind=c_float), intent(out), target :: Af77(LDA,*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_sTile_to_Lapack_c(A,c_loc(Af77),LDA)
      end subroutine PLASMA_sTile_to_Lapack

end module plasma_s
