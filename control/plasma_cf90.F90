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
! @file plasma_cf90.F90
!
!  PLASMA fortran 90 interface
!  PLASMA is a software package provided by Univ. of Tennessee,
!  Univ. of California Berkeley and Univ. of Colorado Denver
!
! @version 2.4.2
! @author Numerical Algorithm Group
! @author Mathieu Faverge
! @date 2011-09-15
! @generated c Thu Sep 15 12:09:27 2011
!
#define PRECISION_c

module plasma_c

!      private PLASMA_cgebrd_c
      interface
         function PLASMA_cgebrd_c(jobu,jobvt,M,N,A,LDA,D,E,U,LDU,VT,LDVT,T) &
          & bind(c, name='PLASMA_cgebrd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgebrd_c
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
         end function PLASMA_cgebrd_c
      end interface

!      private PLASMA_cgelqf_c
      interface
         function PLASMA_cgelqf_c(M,N,A,LDA,T) &
          & bind(c, name='PLASMA_cgelqf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgelqf_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: T
         end function PLASMA_cgelqf_c
      end interface

!      private PLASMA_cgelqs_c
      interface
         function PLASMA_cgelqs_c(M,N,NRHS,A,LDA,T,B,LDB) &
          & bind(c, name='PLASMA_cgelqs')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgelqs_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_cgelqs_c
      end interface

!      private PLASMA_cgels_c
      interface
         function PLASMA_cgels_c(trans,M,N,NRHS,A,LDA,T,B,LDB) &
          & bind(c, name='PLASMA_cgels')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgels_c
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_cgels_c
      end interface

!      private PLASMA_cgemm_c
      interface
         function PLASMA_cgemm_c(transA,transB,M,N,K,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='PLASMA_cgemm')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgemm_c
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: transB
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            complex(kind=c_float_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            complex(kind=c_float_complex), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
         end function PLASMA_cgemm_c
      end interface

!      private PLASMA_cgeqrf_c
      interface
         function PLASMA_cgeqrf_c(M,N,A,LDA,T) &
          & bind(c, name='PLASMA_cgeqrf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgeqrf_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: T
         end function PLASMA_cgeqrf_c
      end interface

!      private PLASMA_cgeqrs_c
      interface
         function PLASMA_cgeqrs_c(M,N,NRHS,A,LDA,T,B,LDB) &
          & bind(c, name='PLASMA_cgeqrs')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgeqrs_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_cgeqrs_c
      end interface

!      private PLASMA_cgesv_c
      interface
         function PLASMA_cgesv_c(N,NRHS,A,LDA,IPIV,B,LDB) &
          & bind(c, name='PLASMA_cgesv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgesv_c
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_cgesv_c
      end interface

!      private PLASMA_cgesv_incpiv_c
      interface
         function PLASMA_cgesv_incpiv_c(N,NRHS,A,LDA,L,IPIV,B,LDB) &
          & bind(c, name='PLASMA_cgesv_incpiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgesv_incpiv_c
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_cgesv_incpiv_c
      end interface

!      private PLASMA_cgesvd_c
      interface
         function PLASMA_cgesvd_c(jobu,jobvt,M,N,A,LDA,S,U,LDU,VT,LDVT,T) &
          & bind(c, name='PLASMA_cgesvd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgesvd_c
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
         end function PLASMA_cgesvd_c
      end interface

!      private PLASMA_cgetrf_c
      interface
         function PLASMA_cgetrf_c(M,N,A,LDA,IPIV) &
          & bind(c, name='PLASMA_cgetrf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgetrf_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: IPIV
         end function PLASMA_cgetrf_c
      end interface

!      private PLASMA_cgetrf_incpiv_c
      interface
         function PLASMA_cgetrf_incpiv_c(M,N,A,LDA,L,IPIV) &
          & bind(c, name='PLASMA_cgetrf_incpiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgetrf_incpiv_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
         end function PLASMA_cgetrf_incpiv_c
      end interface

!      private PLASMA_cgetrs_c
      interface
         function PLASMA_cgetrs_c(trans,N,NRHS,A,LDA,IPIV,B,LDB) &
          & bind(c, name='PLASMA_cgetrs')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgetrs_c
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_cgetrs_c
      end interface

!      private PLASMA_cgetrs_incpiv_c
      interface
         function PLASMA_cgetrs_incpiv_c(trans,N,NRHS,A,LDA,L,IPIV,B,LDB) &
          & bind(c, name='PLASMA_cgetrs_incpiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgetrs_incpiv_c
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_cgetrs_incpiv_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
!      private PLASMA_chemm_c
      interface
         function PLASMA_chemm_c(side,uplo,M,N,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='PLASMA_chemm')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_chemm_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            complex(kind=c_float_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            complex(kind=c_float_complex), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
         end function PLASMA_chemm_c
      end interface

!      private PLASMA_cherk_c
      interface
         function PLASMA_cherk_c(uplo,trans,N,K,alpha,A,LDA,beta,C,LDC) &
          & bind(c, name='PLASMA_cherk')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cherk_c
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
         end function PLASMA_cherk_c
      end interface

!      private PLASMA_cher2k_c
      interface
         function PLASMA_cher2k_c(uplo,trans,N,K,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='PLASMA_cher2k')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cher2k_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            complex(kind=c_float_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            real(kind=c_float), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
         end function PLASMA_cher2k_c
      end interface
#endif

!      private PLASMA_cheev_c
      interface
         function PLASMA_cheev_c(jobz,uplo,N,A,LDA,W,T,Q,LDQ) &
          & bind(c, name='PLASMA_cheev')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cheev_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
         end function PLASMA_cheev_c
      end interface

!      private PLASMA_chegv_c
      interface
         function PLASMA_chegv_c(itype,jobz,uplo,N,A,LDA,B,LDB,W,T,Q,LDQ) &
          & bind(c, name='PLASMA_chegv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_chegv_c
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
         end function PLASMA_chegv_c
      end interface

!      private PLASMA_chegst_c
      interface
         function PLASMA_chegst_c(itype,uplo,N,A,LDA,B,LDB) &
          & bind(c, name='PLASMA_chegst')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_chegst_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_chegst_c
      end interface

!      private PLASMA_chetrd_c
      interface
         function PLASMA_chetrd_c(jobz,uplo,N,A,LDA,D,E,T,Q,LDQ) &
          & bind(c, name='PLASMA_chetrd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_chetrd_c
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
         end function PLASMA_chetrd_c
      end interface

!      private PLASMA_clange_c
      interface
         function PLASMA_clange_c(norm,M,N,A,LDA,work) &
          & bind(c, name='PLASMA_clange')
            use iso_c_binding
            implicit none
            real(kind=c_float) :: PLASMA_clange_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: work
         end function PLASMA_clange_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
!      private PLASMA_clanhe_c
      interface
         function PLASMA_clanhe_c(norm,uplo,N,A,LDA,work) &
          & bind(c, name='PLASMA_clanhe')
            use iso_c_binding
            implicit none
            real(kind=c_float) :: PLASMA_clanhe_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: work
         end function PLASMA_clanhe_c
      end interface
#endif

!      private PLASMA_clansy_c
      interface
         function PLASMA_clansy_c(norm,uplo,N,A,LDA,work) &
          & bind(c, name='PLASMA_clansy')
            use iso_c_binding
            implicit none
            real(kind=c_float) :: PLASMA_clansy_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: work
         end function PLASMA_clansy_c
      end interface

!      private PLASMA_claswp_c
      interface
         function PLASMA_claswp_c(N,A,LDA,K1,K2,IPIV,INCX) &
          & bind(c, name='PLASMA_claswp')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_claswp_c
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            integer(kind=c_int), value :: K1
            integer(kind=c_int), value :: K2
            type(c_ptr), value :: IPIV
            integer(kind=c_int), value :: INCX
         end function PLASMA_claswp_c
      end interface

!      private PLASMA_clauum_c
      interface
         function PLASMA_clauum_c(uplo,N,A,LDA) &
          & bind(c, name='PLASMA_clauum')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_clauum_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
         end function PLASMA_clauum_c
      end interface

! PLASMA_cplghe uses unsigned ints, not generating Fortran interface

! PLASMA_cplgsy uses unsigned ints, not generating Fortran interface

! PLASMA_cplrnt uses unsigned ints, not generating Fortran interface

!      private PLASMA_cposv_c
      interface
         function PLASMA_cposv_c(uplo,N,NRHS,A,LDA,B,LDB) &
          & bind(c, name='PLASMA_cposv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cposv_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_cposv_c
      end interface

!      private PLASMA_cpotrf_c
      interface
         function PLASMA_cpotrf_c(uplo,N,A,LDA) &
          & bind(c, name='PLASMA_cpotrf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cpotrf_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
         end function PLASMA_cpotrf_c
      end interface

!      private PLASMA_cpotri_c
      interface
         function PLASMA_cpotri_c(uplo,N,A,LDA) &
          & bind(c, name='PLASMA_cpotri')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cpotri_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
         end function PLASMA_cpotri_c
      end interface

!      private PLASMA_cpotrs_c
      interface
         function PLASMA_cpotrs_c(uplo,N,NRHS,A,LDA,B,LDB) &
          & bind(c, name='PLASMA_cpotrs')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cpotrs_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_cpotrs_c
      end interface

!      private PLASMA_csymm_c
      interface
         function PLASMA_csymm_c(side,uplo,M,N,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='PLASMA_csymm')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_csymm_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            complex(kind=c_float_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            complex(kind=c_float_complex), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
         end function PLASMA_csymm_c
      end interface

!      private PLASMA_csyrk_c
      interface
         function PLASMA_csyrk_c(uplo,trans,N,K,alpha,A,LDA,beta,C,LDC) &
          & bind(c, name='PLASMA_csyrk')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_csyrk_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            complex(kind=c_float_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            complex(kind=c_float_complex), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
         end function PLASMA_csyrk_c
      end interface

!      private PLASMA_csyr2k_c
      interface
         function PLASMA_csyr2k_c(uplo,trans,N,K,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='PLASMA_csyr2k')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_csyr2k_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            complex(kind=c_float_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            complex(kind=c_float_complex), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
         end function PLASMA_csyr2k_c
      end interface

!      private PLASMA_ctrmm_c
      interface
         function PLASMA_ctrmm_c(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB) &
          & bind(c, name='PLASMA_ctrmm')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ctrmm_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            complex(kind=c_float_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_ctrmm_c
      end interface

!      private PLASMA_ctrsm_c
      interface
         function PLASMA_ctrsm_c(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB) &
          & bind(c, name='PLASMA_ctrsm')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ctrsm_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            complex(kind=c_float_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_ctrsm_c
      end interface

!      private PLASMA_ctrsmpl_c
      interface
         function PLASMA_ctrsmpl_c(N,NRHS,A,LDA,L,IPIV,B,LDB) &
          & bind(c, name='PLASMA_ctrsmpl')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ctrsmpl_c
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_ctrsmpl_c
      end interface

!      private PLASMA_ctrtri_c
      interface
         function PLASMA_ctrtri_c(uplo,diag,N,A,LDA) &
          & bind(c, name='PLASMA_ctrtri')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ctrtri_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: diag
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
         end function PLASMA_ctrtri_c
      end interface

!      private PLASMA_cunglq_c
      interface
         function PLASMA_cunglq_c(M,N,K,A,LDA,T,B,LDB) &
          & bind(c, name='PLASMA_cunglq')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cunglq_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_cunglq_c
      end interface

!      private PLASMA_cungqr_c
      interface
         function PLASMA_cungqr_c(M,N,K,A,LDA,T,B,LDB) &
          & bind(c, name='PLASMA_cungqr')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cungqr_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
         end function PLASMA_cungqr_c
      end interface

!      private PLASMA_cunmlq_c
      interface
         function PLASMA_cunmlq_c(side,trans,M,N,K,A,LDA,T,B,LDB) &
          & bind(c, name='PLASMA_cunmlq')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cunmlq_c
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
         end function PLASMA_cunmlq_c
      end interface

!      private PLASMA_cunmqr_c
      interface
         function PLASMA_cunmqr_c(side,trans,M,N,K,A,LDA,T,B,LDB) &
          & bind(c, name='PLASMA_cunmqr')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cunmqr_c
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
         end function PLASMA_cunmqr_c
      end interface

!      private PLASMA_cgecfi_c
      interface
         function PLASMA_cgecfi_c(m,n,A,fin,imb,inb,fout,omb,onb) &
          & bind(c, name='PLASMA_cgecfi')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgecfi_c
            integer(kind=c_int), value :: m
            integer(kind=c_int), value :: n
            type(c_ptr), value :: A
            integer(kind=c_int), value :: fin
            integer(kind=c_int), value :: imb
            integer(kind=c_int), value :: inb
            integer(kind=c_int), value :: fout
            integer(kind=c_int), value :: omb
            integer(kind=c_int), value :: onb
         end function PLASMA_cgecfi_c
      end interface

!      private PLASMA_cgetmi_c
      interface
         function PLASMA_cgetmi_c(m,n,A,fin,mb,nb) &
          & bind(c, name='PLASMA_cgetmi')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgetmi_c
            integer(kind=c_int), value :: m
            integer(kind=c_int), value :: n
            type(c_ptr), value :: A
            integer(kind=c_int), value :: fin
            integer(kind=c_int), value :: mb
            integer(kind=c_int), value :: nb
         end function PLASMA_cgetmi_c
      end interface

!      private PLASMA_cgebrd_Tile_c
      interface
         function PLASMA_cgebrd_Tile_c(jobu,jobvt,A,D,E,U,VT,T) &
          & bind(c, name='PLASMA_cgebrd_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgebrd_Tile_c
            integer(kind=c_int), value :: jobu
            integer(kind=c_int), value :: jobvt
            type(c_ptr), value :: A
            type(c_ptr), value :: D
            type(c_ptr), value :: E
            type(c_ptr), value :: U
            type(c_ptr), value :: VT
            type(c_ptr), value :: T
         end function PLASMA_cgebrd_Tile_c
      end interface

!      private PLASMA_cgelqf_Tile_c
      interface
         function PLASMA_cgelqf_Tile_c(A,T) &
          & bind(c, name='PLASMA_cgelqf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgelqf_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
         end function PLASMA_cgelqf_Tile_c
      end interface

!      private PLASMA_cgelqs_Tile_c
      interface
         function PLASMA_cgelqs_Tile_c(A,B,T) &
          & bind(c, name='PLASMA_cgelqs_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgelqs_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: T
         end function PLASMA_cgelqs_Tile_c
      end interface

!      private PLASMA_cgels_Tile_c
      interface
         function PLASMA_cgels_Tile_c(trans,A,T,B) &
          & bind(c, name='PLASMA_cgels_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgels_Tile_c
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: T
         end function PLASMA_cgels_Tile_c
      end interface

!      private PLASMA_cgemm_Tile_c
      interface
         function PLASMA_cgemm_Tile_c(transA,transB,alpha,A,B,beta,C) &
          & bind(c, name='PLASMA_cgemm_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgemm_Tile_c
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: transB
            complex(kind=c_float_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            complex(kind=c_float_complex), value :: beta
            type(c_ptr), value :: C
         end function PLASMA_cgemm_Tile_c
      end interface

!      private PLASMA_cgeqrf_Tile_c
      interface
         function PLASMA_cgeqrf_Tile_c(A,T) &
          & bind(c, name='PLASMA_cgeqrf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgeqrf_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
         end function PLASMA_cgeqrf_Tile_c
      end interface

!      private PLASMA_cgeqrs_Tile_c
      interface
         function PLASMA_cgeqrs_Tile_c(A,T,B) &
          & bind(c, name='PLASMA_cgeqrs_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgeqrs_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: T
         end function PLASMA_cgeqrs_Tile_c
      end interface

!      private PLASMA_cgesv_Tile_c
      interface
         function PLASMA_cgesv_Tile_c(A,IPIV,B) &
          & bind(c, name='PLASMA_cgesv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgesv_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
         end function PLASMA_cgesv_Tile_c
      end interface

!      private PLASMA_cgesv_incpiv_Tile_c
      interface
         function PLASMA_cgesv_incpiv_Tile_c(A,L,IPIV,B) &
          & bind(c, name='PLASMA_cgesv_incpiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgesv_incpiv_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
         end function PLASMA_cgesv_incpiv_Tile_c
      end interface

!      private PLASMA_cgesvd_Tile_c
      interface
         function PLASMA_cgesvd_Tile_c(jobu,jobvt,A,S,U,VT,T) &
          & bind(c, name='PLASMA_cgesvd_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgesvd_Tile_c
            integer(kind=c_int), value :: jobu
            integer(kind=c_int), value :: jobvt
            type(c_ptr), value :: A
            type(c_ptr), value :: S
            type(c_ptr), value :: U
            type(c_ptr), value :: VT
            type(c_ptr), value :: T
         end function PLASMA_cgesvd_Tile_c
      end interface

!      private PLASMA_cgetrf_Tile_c
      interface
         function PLASMA_cgetrf_Tile_c(A,IPIV) &
          & bind(c, name='PLASMA_cgetrf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgetrf_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
         end function PLASMA_cgetrf_Tile_c
      end interface

!      private PLASMA_cgetrf_incpiv_Tile_c
      interface
         function PLASMA_cgetrf_incpiv_Tile_c(A,L,IPIV) &
          & bind(c, name='PLASMA_cgetrf_incpiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgetrf_incpiv_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
         end function PLASMA_cgetrf_incpiv_Tile_c
      end interface

!      private PLASMA_cgetrs_Tile_c
      interface
         function PLASMA_cgetrs_Tile_c(trans,A,IPIV,B) &
          & bind(c, name='PLASMA_cgetrs_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgetrs_Tile_c
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
         end function PLASMA_cgetrs_Tile_c
      end interface

!      private PLASMA_cgetrs_incpiv_Tile_c
      interface
         function PLASMA_cgetrs_incpiv_Tile_c(A,L,IPIV,B) &
          & bind(c, name='PLASMA_cgetrs_incpiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cgetrs_incpiv_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
         end function PLASMA_cgetrs_incpiv_Tile_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
!      private PLASMA_chemm_Tile_c
      interface
         function PLASMA_chemm_Tile_c(side,uplo,alpha,A,B,beta,C) &
          & bind(c, name='PLASMA_chemm_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_chemm_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            complex(kind=c_float_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            complex(kind=c_float_complex), value :: beta
            type(c_ptr), value :: C
         end function PLASMA_chemm_Tile_c
      end interface

!      private PLASMA_cherk_Tile_c
      interface
         function PLASMA_cherk_Tile_c(uplo,trans,alpha,A,beta,C) &
          & bind(c, name='PLASMA_cherk_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cherk_Tile_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            real(kind=c_float), value :: alpha
            type(c_ptr), value :: A
            real(kind=c_float), value :: beta
            type(c_ptr), value :: C
         end function PLASMA_cherk_Tile_c
      end interface

!      private PLASMA_cher2k_Tile_c
      interface
         function PLASMA_cher2k_Tile_c(uplo,trans,alpha,A,B,beta,C) &
          & bind(c, name='PLASMA_cher2k_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cher2k_Tile_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            complex(kind=c_float_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            real(kind=c_float), value :: beta
            type(c_ptr), value :: C
         end function PLASMA_cher2k_Tile_c
      end interface
#endif

!      private PLASMA_cheev_Tile_c
      interface
         function PLASMA_cheev_Tile_c(jobz,uplo,A,W,T,Q) &
          & bind(c, name='PLASMA_cheev_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cheev_Tile_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
         end function PLASMA_cheev_Tile_c
      end interface

!      private PLASMA_chegv_Tile_c
      interface
         function PLASMA_chegv_Tile_c(itype,jobz,uplo,A,B,W,T,Q) &
          & bind(c, name='PLASMA_chegv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_chegv_Tile_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
         end function PLASMA_chegv_Tile_c
      end interface

!      private PLASMA_chegst_Tile_c
      interface
         function PLASMA_chegst_Tile_c(itype,uplo,A,B) &
          & bind(c, name='PLASMA_chegst_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_chegst_Tile_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
         end function PLASMA_chegst_Tile_c
      end interface

!      private PLASMA_chetrd_Tile_c
      interface
         function PLASMA_chetrd_Tile_c(jobz,uplo,A,D,E,T,Q) &
          & bind(c, name='PLASMA_chetrd_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_chetrd_Tile_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: D
            type(c_ptr), value :: E
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
         end function PLASMA_chetrd_Tile_c
      end interface

!      private PLASMA_clange_Tile_c
      interface
         function PLASMA_clange_Tile_c(norm,A,work) &
          & bind(c, name='PLASMA_clange_Tile')
            use iso_c_binding
            implicit none
            real(kind=c_float) :: PLASMA_clange_Tile_c
            integer(kind=c_int), value :: norm
            type(c_ptr), value :: A
            type(c_ptr), value :: work
         end function PLASMA_clange_Tile_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
!      private PLASMA_clanhe_Tile_c
      interface
         function PLASMA_clanhe_Tile_c(norm,uplo,A,work) &
          & bind(c, name='PLASMA_clanhe_Tile')
            use iso_c_binding
            implicit none
            real(kind=c_float) :: PLASMA_clanhe_Tile_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: work
         end function PLASMA_clanhe_Tile_c
      end interface
#endif

!      private PLASMA_clansy_Tile_c
      interface
         function PLASMA_clansy_Tile_c(norm,uplo,A,work) &
          & bind(c, name='PLASMA_clansy_Tile')
            use iso_c_binding
            implicit none
            real(kind=c_float) :: PLASMA_clansy_Tile_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: work
         end function PLASMA_clansy_Tile_c
      end interface

!      private PLASMA_claswp_Tile_c
      interface
         function PLASMA_claswp_Tile_c(A,K1,K2,IPIV,INCX) &
          & bind(c, name='PLASMA_claswp_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_claswp_Tile_c
            type(c_ptr), value :: A
            integer(kind=c_int), value :: K1
            integer(kind=c_int), value :: K2
            type(c_ptr), value :: IPIV
            integer(kind=c_int), value :: INCX
         end function PLASMA_claswp_Tile_c
      end interface

!      private PLASMA_clauum_Tile_c
      interface
         function PLASMA_clauum_Tile_c(uplo,A) &
          & bind(c, name='PLASMA_clauum_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_clauum_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
         end function PLASMA_clauum_Tile_c
      end interface

! PLASMA_cplghe_Tile uses unsigned ints, not generating Fortran interface

! PLASMA_cplgsy_Tile uses unsigned ints, not generating Fortran interface

! PLASMA_cplrnt_Tile uses unsigned ints, not generating Fortran interface

!      private PLASMA_cposv_Tile_c
      interface
         function PLASMA_cposv_Tile_c(uplo,A,B) &
          & bind(c, name='PLASMA_cposv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cposv_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
         end function PLASMA_cposv_Tile_c
      end interface

!      private PLASMA_cpotrf_Tile_c
      interface
         function PLASMA_cpotrf_Tile_c(uplo,A) &
          & bind(c, name='PLASMA_cpotrf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cpotrf_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
         end function PLASMA_cpotrf_Tile_c
      end interface

!      private PLASMA_cpotri_Tile_c
      interface
         function PLASMA_cpotri_Tile_c(uplo,A) &
          & bind(c, name='PLASMA_cpotri_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cpotri_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
         end function PLASMA_cpotri_Tile_c
      end interface

!      private PLASMA_cpotrs_Tile_c
      interface
         function PLASMA_cpotrs_Tile_c(uplo,A,B) &
          & bind(c, name='PLASMA_cpotrs_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cpotrs_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
         end function PLASMA_cpotrs_Tile_c
      end interface

!      private PLASMA_csymm_Tile_c
      interface
         function PLASMA_csymm_Tile_c(side,uplo,alpha,A,B,beta,C) &
          & bind(c, name='PLASMA_csymm_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_csymm_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            complex(kind=c_float_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            complex(kind=c_float_complex), value :: beta
            type(c_ptr), value :: C
         end function PLASMA_csymm_Tile_c
      end interface

!      private PLASMA_csyrk_Tile_c
      interface
         function PLASMA_csyrk_Tile_c(uplo,trans,alpha,A,beta,C) &
          & bind(c, name='PLASMA_csyrk_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_csyrk_Tile_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            complex(kind=c_float_complex), value :: alpha
            type(c_ptr), value :: A
            complex(kind=c_float_complex), value :: beta
            type(c_ptr), value :: C
         end function PLASMA_csyrk_Tile_c
      end interface

!      private PLASMA_csyr2k_Tile_c
      interface
         function PLASMA_csyr2k_Tile_c(uplo,trans,alpha,A,B,beta,C) &
          & bind(c, name='PLASMA_csyr2k_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_csyr2k_Tile_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            complex(kind=c_float_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            complex(kind=c_float_complex), value :: beta
            type(c_ptr), value :: C
         end function PLASMA_csyr2k_Tile_c
      end interface

!      private PLASMA_ctrmm_Tile_c
      interface
         function PLASMA_ctrmm_Tile_c(side,uplo,transA,diag,alpha,A,B) &
          & bind(c, name='PLASMA_ctrmm_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ctrmm_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            complex(kind=c_float_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
         end function PLASMA_ctrmm_Tile_c
      end interface

!      private PLASMA_ctrsm_Tile_c
      interface
         function PLASMA_ctrsm_Tile_c(side,uplo,transA,diag,alpha,A,B) &
          & bind(c, name='PLASMA_ctrsm_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ctrsm_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            complex(kind=c_float_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
         end function PLASMA_ctrsm_Tile_c
      end interface

!      private PLASMA_ctrsmpl_Tile_c
      interface
         function PLASMA_ctrsmpl_Tile_c(A,L,IPIV,B) &
          & bind(c, name='PLASMA_ctrsmpl_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ctrsmpl_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
         end function PLASMA_ctrsmpl_Tile_c
      end interface

!      private PLASMA_ctrtri_Tile_c
      interface
         function PLASMA_ctrtri_Tile_c(uplo,diag,A) &
          & bind(c, name='PLASMA_ctrtri_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_ctrtri_Tile_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: diag
            type(c_ptr), value :: A
         end function PLASMA_ctrtri_Tile_c
      end interface

!      private PLASMA_cunglq_Tile_c
      interface
         function PLASMA_cunglq_Tile_c(A,T,B) &
          & bind(c, name='PLASMA_cunglq_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cunglq_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
         end function PLASMA_cunglq_Tile_c
      end interface

!      private PLASMA_cungqr_Tile_c
      interface
         function PLASMA_cungqr_Tile_c(A,T,B) &
          & bind(c, name='PLASMA_cungqr_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cungqr_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
         end function PLASMA_cungqr_Tile_c
      end interface

!      private PLASMA_cunmlq_Tile_c
      interface
         function PLASMA_cunmlq_Tile_c(side,trans,A,T,B) &
          & bind(c, name='PLASMA_cunmlq_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cunmlq_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
         end function PLASMA_cunmlq_Tile_c
      end interface

!      private PLASMA_cunmqr_Tile_c
      interface
         function PLASMA_cunmqr_Tile_c(side,trans,A,T,B) &
          & bind(c, name='PLASMA_cunmqr_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cunmqr_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
         end function PLASMA_cunmqr_Tile_c
      end interface

!      private PLASMA_Alloc_Workspace_cgelqf_c
      interface
         function PLASMA_Alloc_Workspace_cgelqf_c(M,N,T) &
          & bind(c, name='PLASMA_Alloc_Workspace_cgelqf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_cgelqf_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: T ! T is **, so pass by reference
         end function PLASMA_Alloc_Workspace_cgelqf_c
      end interface

!      private PLASMA_Alloc_Workspace_cgels_c
      interface
         function PLASMA_Alloc_Workspace_cgels_c(M,N,T) &
          & bind(c, name='PLASMA_Alloc_Workspace_cgels')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_cgels_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: T ! T is **, so pass by reference
         end function PLASMA_Alloc_Workspace_cgels_c
      end interface

!      private PLASMA_Alloc_Workspace_cgeqrf_c
      interface
         function PLASMA_Alloc_Workspace_cgeqrf_c(M,N,T) &
          & bind(c, name='PLASMA_Alloc_Workspace_cgeqrf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_cgeqrf_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: T ! T is **, so pass by reference
         end function PLASMA_Alloc_Workspace_cgeqrf_c
      end interface

!      private PLASMA_Alloc_Workspace_cgesv_incpiv_c
      interface
         function PLASMA_Alloc_Workspace_cgesv_incpiv_c(N,L,IPIV) &
          & bind(c, name='PLASMA_Alloc_Workspace_cgesv_incpiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_cgesv_incpiv_c
            integer(kind=c_int), value :: N
            type(c_ptr) :: L ! L is **, so pass by reference
            type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         end function PLASMA_Alloc_Workspace_cgesv_incpiv_c
      end interface

!      private PLASMA_Alloc_Workspace_cgetrf_incpiv_c
      interface
         function PLASMA_Alloc_Workspace_cgetrf_incpiv_c(M,N,L,IPIV) &
          & bind(c, name='PLASMA_Alloc_Workspace_cgetrf_incpiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_cgetrf_incpiv_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: L ! L is **, so pass by reference
            type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         end function PLASMA_Alloc_Workspace_cgetrf_incpiv_c
      end interface

!      private PLASMA_Alloc_Workspace_cgeev_c
      interface
         function PLASMA_Alloc_Workspace_cgeev_c(N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_cgeev')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_cgeev_c
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_cgeev_c
      end interface

!      private PLASMA_Alloc_Workspace_cgebrd_c
      interface
         function PLASMA_Alloc_Workspace_cgebrd_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_cgebrd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_cgebrd_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_cgebrd_c
      end interface

!      private PLASMA_Alloc_Workspace_cgesvd_c
      interface
         function PLASMA_Alloc_Workspace_cgesvd_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_cgesvd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_cgesvd_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_cgesvd_c
      end interface

!      private PLASMA_Alloc_Workspace_cheev_c
      interface
         function PLASMA_Alloc_Workspace_cheev_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_cheev')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_cheev_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_cheev_c
      end interface

!      private PLASMA_Alloc_Workspace_chegv_c
      interface
         function PLASMA_Alloc_Workspace_chegv_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_chegv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_chegv_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_chegv_c
      end interface

!      private PLASMA_Alloc_Workspace_chetrd_c
      interface
         function PLASMA_Alloc_Workspace_chetrd_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_chetrd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_chetrd_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_chetrd_c
      end interface

!      private PLASMA_Alloc_Workspace_cgelqf_Tile_c
      interface
         function PLASMA_Alloc_Workspace_cgelqf_Tile_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_cgelqf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_cgelqf_Tile_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_cgelqf_Tile_c
      end interface

!      private PLASMA_Alloc_Workspace_cgels_Tile_c
      interface
         function PLASMA_Alloc_Workspace_cgels_Tile_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_cgels_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_cgels_Tile_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_cgels_Tile_c
      end interface

!      private PLASMA_Alloc_Workspace_cgeqrf_Tile_c
      interface
         function PLASMA_Alloc_Workspace_cgeqrf_Tile_c(M,N,descT) &
          & bind(c, name='PLASMA_Alloc_Workspace_cgeqrf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_cgeqrf_Tile_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
         end function PLASMA_Alloc_Workspace_cgeqrf_Tile_c
      end interface

!      private PLASMA_Alloc_Workspace_cgesv_incpiv_Tile_c
      interface
         function PLASMA_Alloc_Workspace_cgesv_incpiv_Tile_c(N,descL,IPIV) &
          & bind(c, name='PLASMA_Alloc_Workspace_cgesv_incpiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_cgesv_incpiv_Tile_c
            integer(kind=c_int), value :: N
            type(c_ptr) :: descL ! descL is **, so pass by reference
            type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         end function PLASMA_Alloc_Workspace_cgesv_incpiv_Tile_c
      end interface

!      private PLASMA_Alloc_Workspace_cgetrf_incpiv_Tile_c
      interface
         function PLASMA_Alloc_Workspace_cgetrf_incpiv_Tile_c(N,descL,IPIV) &
          & bind(c, name='PLASMA_Alloc_Workspace_cgetrf_incpiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_Alloc_Workspace_cgetrf_incpiv_Tile_c
            integer(kind=c_int), value :: N
            type(c_ptr) :: descL ! descL is **, so pass by reference
            type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         end function PLASMA_Alloc_Workspace_cgetrf_incpiv_Tile_c
      end interface

!      private PLASMA_cLapack_to_Tile_c
      interface
         function PLASMA_cLapack_to_Tile_c(Af77,LDA,A) &
          & bind(c, name='PLASMA_cLapack_to_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cLapack_to_Tile_c
            type(c_ptr), value :: Af77
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: A
         end function PLASMA_cLapack_to_Tile_c
      end interface

!      private PLASMA_cTile_to_Lapack_c
      interface
         function PLASMA_cTile_to_Lapack_c(A,Af77,LDA) &
          & bind(c, name='PLASMA_cTile_to_Lapack')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_cTile_to_Lapack_c
            type(c_ptr), value :: A
            type(c_ptr), value :: Af77
            integer(kind=c_int), value :: LDA
         end function PLASMA_cTile_to_Lapack_c
      end interface

  contains

      subroutine PLASMA_cgebrd(jobu,jobvt,M,N,A,LDA,D,E,U,LDU,VT,LDVT,T,info)
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
         complex(kind=c_float_complex), intent(inout), target :: A(LDA,*)
         real(kind=c_float), intent(out), target :: D(*)
         real(kind=c_float), intent(out), target :: E(*)
         complex(kind=c_float_complex), intent(out), target :: U(LDU,*)
         complex(kind=c_float_complex), intent(out), target :: VT(LDVT,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cgebrd_c(jobu,jobvt,M,N,c_loc(A),LDA,c_loc(D),c_loc(E),c_loc(U),LDU,c_loc(VT),LDVT,T)
      end subroutine PLASMA_cgebrd

      subroutine PLASMA_cgelqf(M,N,A,LDA,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_float_complex), intent(inout), target :: A(LDA,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cgelqf_c(M,N,c_loc(A),LDA,T)
      end subroutine PLASMA_cgelqf

      subroutine PLASMA_cgelqs(M,N,NRHS,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         complex(kind=c_float_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_float_complex), intent(inout), target :: B(LDB,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cgelqs_c(M,N,NRHS,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_cgelqs

      subroutine PLASMA_cgels(trans,M,N,NRHS,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: trans
         complex(kind=c_float_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_float_complex), intent(inout), target :: B(LDB,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cgels_c(trans,M,N,NRHS,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_cgels

      subroutine PLASMA_cgemm(transA,transB,M,N,K,alpha,A,LDA,B,LDB,beta,C,LDC,info)
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
         complex(kind=c_float_complex), intent(in) :: alpha
         complex(kind=c_float_complex), intent(in) :: beta
         complex(kind=c_float_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_float_complex), intent(in), target :: B(LDB,*)
         complex(kind=c_float_complex), intent(inout), target :: C(LDC,*)
         info = PLASMA_cgemm_c(transA,transB,M,N,K,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
      end subroutine PLASMA_cgemm

      subroutine PLASMA_cgeqrf(M,N,A,LDA,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_float_complex), intent(inout), target :: A(LDA,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cgeqrf_c(M,N,c_loc(A),LDA,T)
      end subroutine PLASMA_cgeqrf

      subroutine PLASMA_cgeqrs(M,N,NRHS,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         complex(kind=c_float_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_float_complex), intent(inout), target :: B(LDB,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cgeqrs_c(M,N,NRHS,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_cgeqrs

      subroutine PLASMA_cgesv(N,NRHS,A,LDA,IPIV,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(out), target :: IPIV(*)
         complex(kind=c_float_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_float_complex), intent(inout), target :: B(LDB,*)
         info = PLASMA_cgesv_c(N,NRHS,c_loc(A),LDA,c_loc(IPIV),c_loc(B),LDB)
      end subroutine PLASMA_cgesv

      subroutine PLASMA_cgesv_incpiv(N,NRHS,A,LDA,L,IPIV,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         complex(kind=c_float_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_float_complex), intent(inout), target :: B(LDB,*)
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cgesv_incpiv_c(N,NRHS,c_loc(A),LDA,L,IPIV,c_loc(B),LDB)
      end subroutine PLASMA_cgesv_incpiv

      subroutine PLASMA_cgesvd(jobu,jobvt,M,N,A,LDA,S,U,LDU,VT,LDVT,T,info)
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
         complex(kind=c_float_complex), intent(inout), target :: A(LDA,*)
         real(kind=c_float), intent(out), target :: S(*)
         complex(kind=c_float_complex), intent(out), target :: U(LDU,*)
         complex(kind=c_float_complex), intent(out), target :: VT(LDVT,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cgesvd_c(jobu,jobvt,M,N,c_loc(A),LDA,c_loc(S),c_loc(U),LDU,c_loc(VT),LDVT,T)
      end subroutine PLASMA_cgesvd

      subroutine PLASMA_cgetrf(M,N,A,LDA,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(out), target :: IPIV(*)
         complex(kind=c_float_complex), intent(inout), target :: A(LDA,*)
         info = PLASMA_cgetrf_c(M,N,c_loc(A),LDA,c_loc(IPIV))
      end subroutine PLASMA_cgetrf

      subroutine PLASMA_cgetrf_incpiv(M,N,A,LDA,L,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_float_complex), intent(inout), target :: A(LDA,*)
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L    ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cgetrf_incpiv_c(M,N,c_loc(A),LDA,L,IPIV)
      end subroutine PLASMA_cgetrf_incpiv

      subroutine PLASMA_cgetrs(trans,N,NRHS,A,LDA,IPIV,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in), target :: IPIV(*)
         complex(kind=c_float_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_float_complex), intent(inout), target :: B(LDB,*)
         info = PLASMA_cgetrs_c(trans,N,NRHS,c_loc(A),LDA,c_loc(IPIV),c_loc(B),LDB)
      end subroutine PLASMA_cgetrs

      subroutine PLASMA_cgetrs_incpiv(trans,N,NRHS,A,LDA,L,IPIV,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: trans
         complex(kind=c_float_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_float_complex), intent(inout), target :: B(LDB,*)
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L    ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cgetrs_incpiv_c(trans,N,NRHS,c_loc(A),LDA,L,IPIV,c_loc(B),LDB)
      end subroutine PLASMA_cgetrs_incpiv

#if defined(PRECISION_z) || defined(PRECISION_c)
      subroutine PLASMA_chemm(side,uplo,M,N,alpha,A,LDA,B,LDB,beta,C,LDC,info)
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
         complex(kind=c_float_complex), intent(in) :: alpha
         complex(kind=c_float_complex), intent(in) :: beta
         complex(kind=c_float_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_float_complex), intent(in), target :: B(LDB,*)
         complex(kind=c_float_complex), intent(inout), target :: C(LDC,*)
         info = PLASMA_chemm_c(side,uplo,M,N,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
      end subroutine PLASMA_chemm

      subroutine PLASMA_cherk(uplo,trans,N,K,alpha,A,LDA,beta,C,LDC,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDC
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_float_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_float_complex), intent(inout), target :: C(LDC,*)
         real(kind=c_float), intent(in) :: alpha
         real(kind=c_float), intent(in) :: beta
         info = PLASMA_cherk_c(uplo,trans,N,K,alpha,c_loc(A),LDA,beta,c_loc(C),LDC)
      end subroutine PLASMA_cherk

      subroutine PLASMA_cher2k(uplo,trans,N,K,alpha,A,LDA,B,LDB,beta,C,LDC,info)
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
         complex(kind=c_float_complex), intent(in) :: alpha
         complex(kind=c_float_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_float_complex), intent(in), target :: B(LDB,*)
         complex(kind=c_float_complex), intent(inout), target :: C(LDC,*)
         real(kind=c_float), intent(in) :: beta
         info = PLASMA_cher2k_c(uplo,trans,N,K,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
      end subroutine PLASMA_cher2k
#endif

      subroutine PLASMA_cheev(jobz,uplo,N,A,LDA,W,T,Q,LDQ,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDQ
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_float_complex), intent(inout), target :: A(LDA,*)
         real(kind=c_float), intent(out), target :: W(*)
         complex(kind=c_float_complex), intent(out), target :: Q(LDQ,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cheev_c(jobz,uplo,N,c_loc(A),LDA,c_loc(W),T,c_loc(Q),LDQ)
      end subroutine PLASMA_cheev

      subroutine PLASMA_chegv(itype,jobz,uplo,N,A,LDA,B,LDB,W,T,Q,LDQ,info)
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
         complex(kind=c_float_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_float_complex), intent(inout), target :: B(LDB,*)
         real(kind=c_float), intent(out), target :: W(*)
         complex(kind=c_float_complex), intent(out), target :: Q(LDQ,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_chegv_c(itype,jobz,uplo,N,c_loc(A),LDA,c_loc(B),LDB,c_loc(W),T,c_loc(Q),LDQ)
      end subroutine PLASMA_chegv

      subroutine PLASMA_chegst(itype,uplo,N,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_float_complex), intent(in), target :: B(LDB,*)
         complex(kind=c_float_complex), intent(inout), target :: A(LDA,*)
         info = PLASMA_chegst_c(itype,uplo,N,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine PLASMA_chegst

      subroutine PLASMA_chetrd(jobz,uplo,N,A,LDA,D,E,T,Q,LDQ,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDQ
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_float_complex), intent(inout), target :: A(LDA,*)
         real(kind=c_float), intent(out), target :: D(*)
         real(kind=c_float), intent(out), target :: E(*)
         complex(kind=c_float_complex), intent(out), target :: Q(LDQ,*)
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_chetrd_c(jobz,uplo,N,c_loc(A),LDA,c_loc(D),c_loc(E),T,c_loc(Q),LDQ)
      end subroutine PLASMA_chetrd

      function PLASMA_clange(norm,M,N,A,LDA,work)
         use iso_c_binding
         implicit none
         real(kind=c_float) :: PLASMA_clange
         real(kind=c_float), intent(inout), target :: work(*)
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_float_complex), intent(in), target :: A(LDA,*)
         PLASMA_clange = PLASMA_clange_c(norm,M,N,c_loc(A),LDA,c_loc(work))
      end function PLASMA_clange

#if defined(PRECISION_z) || defined(PRECISION_c)
      function PLASMA_clanhe(norm,uplo,N,A,LDA,work)
         use iso_c_binding
         implicit none
         real(kind=c_float) :: PLASMA_clanhe
         real(kind=c_float), intent(inout), target :: work(*)
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_float_complex), intent(in), target :: A(LDA,*)
         PLASMA_clanhe = PLASMA_clanhe_c(norm,uplo,N,c_loc(A),LDA,c_loc(work))
      end function PLASMA_clanhe
#endif

      function PLASMA_clansy(norm,uplo,N,A,LDA,work)
         use iso_c_binding
         implicit none
         real(kind=c_float) :: PLASMA_clansy
         real(kind=c_float), intent(inout), target :: work(*)
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_float_complex), intent(inout), target :: A(LDA,*)
         PLASMA_clansy = PLASMA_clansy_c(norm,uplo,N,c_loc(A),LDA,c_loc(work))
      end function PLASMA_clansy

      subroutine PLASMA_claswp(N,A,LDA,K1,K2,IPIV,INCX,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: INCX
         integer(kind=c_int), intent(in) :: K1
         integer(kind=c_int), intent(in) :: K2
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in), target :: IPIV(*)
         complex(kind=c_float_complex), intent(inout), target :: A(LDA,*)
         info = PLASMA_claswp_c(N,c_loc(A),LDA,K1,K2,c_loc(IPIV),INCX)
      end subroutine PLASMA_claswp

      subroutine PLASMA_clauum(uplo,N,A,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_float_complex), intent(inout), target :: A(LDA,*)
         info = PLASMA_clauum_c(uplo,N,c_loc(A),LDA)
      end subroutine PLASMA_clauum

      subroutine PLASMA_cposv(uplo,N,NRHS,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_float_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_float_complex), intent(inout), target :: B(LDB,*)
         info = PLASMA_cposv_c(uplo,N,NRHS,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine PLASMA_cposv

      subroutine PLASMA_cpotrf(uplo,N,A,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_float_complex), intent(inout), target :: A(LDA,*)
         info = PLASMA_cpotrf_c(uplo,N,c_loc(A),LDA)
      end subroutine PLASMA_cpotrf

      subroutine PLASMA_cpotri(uplo,N,A,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_float_complex), intent(inout), target :: A(LDA,*)
         info = PLASMA_cpotri_c(uplo,N,c_loc(A),LDA)
      end subroutine PLASMA_cpotri

      subroutine PLASMA_cpotrs(uplo,N,NRHS,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_float_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_float_complex), intent(inout), target :: B(LDB,*)
         info = PLASMA_cpotrs_c(uplo,N,NRHS,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine PLASMA_cpotrs

      subroutine PLASMA_csymm(side,uplo,M,N,alpha,A,LDA,B,LDB,beta,C,LDC,info)
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
         complex(kind=c_float_complex), intent(in) :: alpha
         complex(kind=c_float_complex), intent(in) :: beta
         complex(kind=c_float_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_float_complex), intent(in), target :: B(LDB,*)
         complex(kind=c_float_complex), intent(inout), target :: C(LDC,*)
         info = PLASMA_csymm_c(side,uplo,M,N,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
      end subroutine PLASMA_csymm

      subroutine PLASMA_csyrk(uplo,trans,N,K,alpha,A,LDA,beta,C,LDC,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDC
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_float_complex), intent(in) :: alpha
         complex(kind=c_float_complex), intent(in) :: beta
         complex(kind=c_float_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_float_complex), intent(inout), target :: C(LDC,*)
         info = PLASMA_csyrk_c(uplo,trans,N,K,alpha,c_loc(A),LDA,beta,c_loc(C),LDC)
      end subroutine PLASMA_csyrk

      subroutine PLASMA_csyr2k(uplo,trans,N,K,alpha,A,LDA,B,LDB,beta,C,LDC,info)
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
         complex(kind=c_float_complex), intent(in) :: alpha
         complex(kind=c_float_complex), intent(in) :: beta
         complex(kind=c_float_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_float_complex), intent(in), target :: B(LDB,*)
         complex(kind=c_float_complex), intent(inout), target :: C(LDC,*)
         info = PLASMA_csyr2k_c(uplo,trans,N,K,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
      end subroutine PLASMA_csyr2k

      subroutine PLASMA_ctrmm(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB,info)
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
         complex(kind=c_float_complex), intent(in) :: alpha
         complex(kind=c_float_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_float_complex), intent(inout), target :: B(LDB,*)
         info = PLASMA_ctrmm_c(side,uplo,transA,diag,N,NRHS,alpha,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine PLASMA_ctrmm

      subroutine PLASMA_ctrsm(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB,info)
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
         complex(kind=c_float_complex), intent(in) :: alpha
         complex(kind=c_float_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_float_complex), intent(inout), target :: B(LDB,*)
         info = PLASMA_ctrsm_c(side,uplo,transA,diag,N,NRHS,alpha,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine PLASMA_ctrsm

      subroutine PLASMA_ctrsmpl(N,NRHS,A,LDA,L,IPIV,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         complex(kind=c_float_complex), intent(in), target :: A(LDA,*)
         type(c_ptr), value :: L    ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         complex(kind=c_float_complex), intent(inout), target :: B(LDB,*)
         info = PLASMA_ctrsmpl_c(N,NRHS,c_loc(A),LDA,L,IPIV,c_loc(B),LDB)
      end subroutine PLASMA_ctrsmpl

      subroutine PLASMA_ctrtri(uplo,diag,N,A,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_float_complex), intent(inout), target :: A(LDA,*)
         info = PLASMA_ctrtri_c(uplo,diag,N,c_loc(A),LDA)
      end subroutine PLASMA_ctrtri

      subroutine PLASMA_cunglq(M,N,K,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_float_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_float_complex), intent(out), target :: B(LDB,*)
         info = PLASMA_cunglq_c(M,N,K,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_cunglq

      subroutine PLASMA_cungqr(M,N,K,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_float_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_float_complex), intent(out), target :: B(LDB,*)
         info = PLASMA_cungqr_c(M,N,K,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_cungqr

      subroutine PLASMA_cunmlq(side,trans,M,N,K,A,LDA,T,B,LDB,info)
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
         complex(kind=c_float_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_float_complex), intent(inout), target :: B(LDB,*)
         info = PLASMA_cunmlq_c(side,trans,M,N,K,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_cunmlq

      subroutine PLASMA_cunmqr(side,trans,M,N,K,A,LDA,T,B,LDB,info)
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
         complex(kind=c_float_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_float_complex), intent(inout), target :: B(LDB,*)
         info = PLASMA_cunmqr_c(side,trans,M,N,K,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine PLASMA_cunmqr

      subroutine PLASMA_cgecfi(m,n,A,fin,imb,inb,fout,omb,onb,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         complex(kind=c_float_complex), intent(inout), target :: A(*)
         integer(kind=c_int), intent(in) :: fin
         integer(kind=c_int), intent(in) :: fout
         integer(kind=c_int), intent(in) :: imb
         integer(kind=c_int), intent(in) :: inb
         integer(kind=c_int), intent(in) :: omb
         integer(kind=c_int), intent(in) :: onb
         integer(kind=c_int), intent(in) :: m
         integer(kind=c_int), intent(in) :: n
         info = PLASMA_cgecfi_c(m,n,c_loc(A),fin,imb,inb,fout,omb,onb)
      end subroutine PLASMA_cgecfi

      subroutine PLASMA_cgetmi(m,n,A,fin,mb,nb,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         complex(kind=c_float_complex), intent(inout), target :: A(*)
         integer(kind=c_int), intent(in) :: fin
         integer(kind=c_int), intent(in) :: mb
         integer(kind=c_int), intent(in) :: nb
         integer(kind=c_int), intent(in) :: m
         integer(kind=c_int), intent(in) :: n
         info = PLASMA_cgetmi_c(m,n,c_loc(A),fin,mb,nb)
      end subroutine PLASMA_cgetmi

      subroutine PLASMA_cgebrd_Tile(jobu,jobvt,A,D,E,U,VT,T,info)
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
         info = PLASMA_cgebrd_Tile_c(jobu,jobvt,A,c_loc(D),c_loc(E),U,VT,T)
      end subroutine PLASMA_cgebrd_Tile

      subroutine PLASMA_cgelqf_Tile(A,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cgelqf_Tile_c(A,T)
      end subroutine PLASMA_cgelqf_Tile

      subroutine PLASMA_cgelqs_Tile(A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cgelqs_Tile_c(A,T,B)
      end subroutine PLASMA_cgelqs_Tile

      subroutine PLASMA_cgels_Tile(trans,A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cgels_Tile_c(trans,A,T,B)
      end subroutine PLASMA_cgels_Tile

      subroutine PLASMA_cgemm_Tile(transA,transB,alpha,A,B,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: transB
         complex(kind=c_float_complex), intent(in) :: alpha
         complex(kind=c_float_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cgemm_Tile_c(transA,transB,alpha,A,B,beta,C)
      end subroutine PLASMA_cgemm_Tile

      subroutine PLASMA_cgeqrf_Tile(A,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cgeqrf_Tile_c(A,T)
      end subroutine PLASMA_cgeqrf_Tile

      subroutine PLASMA_cgeqrs_Tile(A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cgeqrs_Tile_c(A,T,B)
      end subroutine PLASMA_cgeqrs_Tile

      subroutine PLASMA_cgesv_Tile(A,IPIV,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: IPIV(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cgesv_Tile_c(A,c_loc(IPIV),B)
      end subroutine PLASMA_cgesv_Tile

      subroutine PLASMA_cgesv_incpiv_Tile(A,L,IPIV,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cgesv_incpiv_Tile_c(A,L,IPIV,B)
      end subroutine PLASMA_cgesv_incpiv_Tile

      subroutine PLASMA_cgesvd_Tile(jobu,jobvt,A,S,U,VT,T,info)
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
         info = PLASMA_cgesvd_Tile_c(jobu,jobvt,A,c_loc(S),U,VT,T)
      end subroutine PLASMA_cgesvd_Tile

      subroutine PLASMA_cgetrf_Tile(A,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: IPIV(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cgetrf_Tile_c(A,c_loc(IPIV))
      end subroutine PLASMA_cgetrf_Tile

      subroutine PLASMA_cgetrf_incpiv_Tile(A,L,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cgetrf_incpiv_Tile_c(A,L,IPIV)
      end subroutine PLASMA_cgetrf_incpiv_Tile

      subroutine PLASMA_cgetrs_Tile(trans,A,IPIV,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         integer(kind=c_int), intent(in), target :: IPIV(*)
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cgetrs_Tile_c(trans,A,c_loc(IPIV),B)
      end subroutine PLASMA_cgetrs_Tile

      subroutine PLASMA_cgetrs_incpiv_Tile(A,L,IPIV,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cgetrs_incpiv_Tile_c(A,L,IPIV,B)
      end subroutine PLASMA_cgetrs_incpiv_Tile

#if defined(PRECISION_z) || defined(PRECISION_c)
      subroutine PLASMA_chemm_Tile(side,uplo,alpha,A,B,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_float_complex), intent(in) :: alpha
         complex(kind=c_float_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_chemm_Tile_c(side,uplo,alpha,A,B,beta,C)
      end subroutine PLASMA_chemm_Tile

      subroutine PLASMA_cherk_Tile(uplo,trans,alpha,A,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_float), intent(in) :: alpha
         real(kind=c_float), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cherk_Tile_c(uplo,trans,alpha,A,beta,C)
      end subroutine PLASMA_cherk_Tile

      subroutine PLASMA_cher2k_Tile(uplo,trans,alpha,A,B,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_float_complex), intent(in) :: alpha
         real(kind=c_float), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cher2k_Tile_c(uplo,trans,alpha,A,B,beta,C)
      end subroutine PLASMA_cher2k_Tile
#endif

      subroutine PLASMA_cheev_Tile(jobz,uplo,A,W,T,Q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_float), intent(out), target :: W(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: Q ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cheev_Tile_c(jobz,uplo,A,c_loc(W),T,Q)
      end subroutine PLASMA_cheev_Tile

      subroutine PLASMA_chegv_Tile(itype,jobz,uplo,A,B,W,T,Q,info)
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
         info = PLASMA_chegv_Tile_c(itype,jobz,uplo,A,B,c_loc(W),T,Q)
      end subroutine PLASMA_chegv_Tile

      subroutine PLASMA_chegst_Tile(itype,uplo,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_chegst_Tile_c(itype,uplo,A,B)
      end subroutine PLASMA_chegst_Tile

      subroutine PLASMA_chetrd_Tile(jobz,uplo,A,D,E,T,Q,info)
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
         info = PLASMA_chetrd_Tile_c(jobz,uplo,A,c_loc(D),c_loc(E),T,Q)
      end subroutine PLASMA_chetrd_Tile

      function PLASMA_clange_Tile(norm,A,work)
         use iso_c_binding
         implicit none
         real(kind=c_float) :: PLASMA_clange_Tile
         integer(kind=c_int), intent(in) :: norm
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         real(kind=c_float), intent(inout), target :: work(*)
         PLASMA_clange_Tile = PLASMA_clange_Tile_c(norm,A,c_loc(work))
       end function PLASMA_clange_Tile

#if defined(PRECISION_z) || defined(PRECISION_c)
      function PLASMA_clanhe_Tile(norm,uplo,A,work)
         use iso_c_binding
         implicit none
         real(kind=c_float) :: PLASMA_clanhe_Tile
         real(kind=c_float), intent(inout), target :: work(*)
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         PLASMA_clanhe_Tile = PLASMA_clanhe_Tile_c(norm,uplo,A,c_loc(work))
      end function PLASMA_clanhe_Tile
#endif

      function PLASMA_clansy_Tile(norm,uplo,A,work)
         use iso_c_binding
         implicit none
         real(kind=c_float) :: PLASMA_clansy_Tile
         real(kind=c_float), intent(inout), target :: work(*)
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         PLASMA_clansy_Tile = PLASMA_clansy_Tile_c(norm,uplo,A,c_loc(work))
      end function PLASMA_clansy_Tile

      subroutine PLASMA_claswp_Tile(A,K1,K2,IPIV,INCX,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: INCX
         integer(kind=c_int), intent(in) :: K1
         integer(kind=c_int), intent(in) :: K2
         integer(kind=c_int), intent(in), target :: IPIV(*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_claswp_Tile_c(A,K1,K2,c_loc(IPIV),INCX)
      end subroutine PLASMA_claswp_Tile

      subroutine PLASMA_clauum_Tile(uplo,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_clauum_Tile_c(uplo,A)
      end subroutine PLASMA_clauum_Tile

      subroutine PLASMA_cposv_Tile(uplo,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cposv_Tile_c(uplo,A,B)
      end subroutine PLASMA_cposv_Tile

      subroutine PLASMA_cpotrf_Tile(uplo,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cpotrf_Tile_c(uplo,A)
      end subroutine PLASMA_cpotrf_Tile

      subroutine PLASMA_cpotri_Tile(uplo,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cpotri_Tile_c(uplo,A)
      end subroutine PLASMA_cpotri_Tile

      subroutine PLASMA_cpotrs_Tile(uplo,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cpotrs_Tile_c(uplo,A,B)
      end subroutine PLASMA_cpotrs_Tile

      subroutine PLASMA_csymm_Tile(side,uplo,alpha,A,B,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_float_complex), intent(in) :: alpha
         complex(kind=c_float_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_csymm_Tile_c(side,uplo,alpha,A,B,beta,C)
      end subroutine PLASMA_csymm_Tile

      subroutine PLASMA_csyrk_Tile(uplo,trans,alpha,A,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_float_complex), intent(in) :: alpha
         complex(kind=c_float_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_csyrk_Tile_c(uplo,trans,alpha,A,beta,C)
      end subroutine PLASMA_csyrk_Tile

      subroutine PLASMA_csyr2k_Tile(uplo,trans,alpha,A,B,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_float_complex), intent(in) :: alpha
         complex(kind=c_float_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_csyr2k_Tile_c(uplo,trans,alpha,A,B,beta,C)
      end subroutine PLASMA_csyr2k_Tile

      subroutine PLASMA_ctrmm_Tile(side,uplo,transA,diag,alpha,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_float_complex), intent(in) :: alpha
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ctrmm_Tile_c(side,uplo,transA,diag,alpha,A,B)
      end subroutine PLASMA_ctrmm_Tile

      subroutine PLASMA_ctrsm_Tile(side,uplo,transA,diag,alpha,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_float_complex), intent(in) :: alpha
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ctrsm_Tile_c(side,uplo,transA,diag,alpha,A,B)
      end subroutine PLASMA_ctrsm_Tile

      subroutine PLASMA_ctrsmpl_Tile(A,L,IPIV,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ctrsmpl_Tile_c(A,L,IPIV,B)
      end subroutine PLASMA_ctrsmpl_Tile

      subroutine PLASMA_ctrtri_Tile(uplo,diag,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_ctrtri_Tile_c(uplo,diag,A)
      end subroutine PLASMA_ctrtri_Tile

      subroutine PLASMA_cunglq_Tile(A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cunglq_Tile_c(A,T,B)
      end subroutine PLASMA_cunglq_Tile

      subroutine PLASMA_cungqr_Tile(A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cungqr_Tile_c(A,T,B)
      end subroutine PLASMA_cungqr_Tile

      subroutine PLASMA_cunmlq_Tile(side,trans,A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cunmlq_Tile_c(side,trans,A,T,B)
      end subroutine PLASMA_cunmlq_Tile

      subroutine PLASMA_cunmqr_Tile(side,trans,A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cunmqr_Tile_c(side,trans,A,T,B)
      end subroutine PLASMA_cunmqr_Tile

      subroutine PLASMA_Alloc_Workspace_cgelqf(M,N,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         info = PLASMA_Alloc_Workspace_cgelqf_c(M,N,T)
      end subroutine PLASMA_Alloc_Workspace_cgelqf

      subroutine PLASMA_Alloc_Workspace_cgels(M,N,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         info = PLASMA_Alloc_Workspace_cgels_c(M,N,T)
      end subroutine PLASMA_Alloc_Workspace_cgels

      subroutine PLASMA_Alloc_Workspace_cgeqrf(M,N,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         info = PLASMA_Alloc_Workspace_cgeqrf_c(M,N,T)
      end subroutine PLASMA_Alloc_Workspace_cgeqrf

      subroutine PLASMA_Alloc_Workspace_cgesv_incpiv(N,L,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         type(c_ptr) :: L ! L is **, so pass by reference
         info = PLASMA_Alloc_Workspace_cgesv_incpiv_c(N,L,IPIV)
      end subroutine PLASMA_Alloc_Workspace_cgesv_incpiv

      subroutine PLASMA_Alloc_Workspace_cgetrf_incpiv(M,N,L,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         type(c_ptr) :: L ! L is **, so pass by reference
         info = PLASMA_Alloc_Workspace_cgetrf_incpiv_c(M,N,L,IPIV)
      end subroutine PLASMA_Alloc_Workspace_cgetrf_incpiv

      subroutine PLASMA_Alloc_Workspace_cgeev(N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_cgeev_c(N,descT)
      end subroutine PLASMA_Alloc_Workspace_cgeev

      subroutine PLASMA_Alloc_Workspace_cgebrd(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_cgebrd_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_cgebrd

      subroutine PLASMA_Alloc_Workspace_cgesvd(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_cgesvd_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_cgesvd

      subroutine PLASMA_Alloc_Workspace_cheev(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_cheev_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_cheev

      subroutine PLASMA_Alloc_Workspace_chegv(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_chegv_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_chegv

      subroutine PLASMA_Alloc_Workspace_chetrd(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_chetrd_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_chetrd

      subroutine PLASMA_Alloc_Workspace_cgelqf_Tile(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_cgelqf_Tile_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_cgelqf_Tile

      subroutine PLASMA_Alloc_Workspace_cgels_Tile(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_cgels_Tile_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_cgels_Tile

      subroutine PLASMA_Alloc_Workspace_cgeqrf_Tile(M,N,descT,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         info = PLASMA_Alloc_Workspace_cgeqrf_Tile_c(M,N,descT)
      end subroutine PLASMA_Alloc_Workspace_cgeqrf_Tile

      subroutine PLASMA_Alloc_Workspace_cgesv_incpiv_Tile(N,descL,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         type(c_ptr) :: descL ! descL is **, so pass by reference
         info = PLASMA_Alloc_Workspace_cgesv_incpiv_Tile_c(N,descL,IPIV)
      end subroutine PLASMA_Alloc_Workspace_cgesv_incpiv_Tile

      subroutine PLASMA_Alloc_Workspace_cgetrf_incpiv_Tile(N,descL,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         type(c_ptr) :: descL ! descL is **, so pass by reference
         info = PLASMA_Alloc_Workspace_cgetrf_incpiv_Tile_c(N,descL,IPIV)
      end subroutine PLASMA_Alloc_Workspace_cgetrf_incpiv_Tile

      subroutine PLASMA_cLapack_to_Tile(Af77,LDA,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         complex(kind=c_float_complex), intent(in), target :: Af77(LDA,*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cLapack_to_Tile_c(c_loc(Af77),LDA,A)
      end subroutine PLASMA_cLapack_to_Tile

      subroutine PLASMA_cTile_to_Lapack(A,Af77,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         complex(kind=c_float_complex), intent(out), target :: Af77(LDA,*)
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_cTile_to_Lapack_c(A,c_loc(Af77),LDA)
      end subroutine PLASMA_cTile_to_Lapack

end module plasma_c
