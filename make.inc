###
#
# @file make.inc.gcc.atlas
#
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.4.5
# @author Piotr Luszczek
# @author Mathieu Faverge
# @date 2010-11-15
#
###
#
# OS:       Linux
# Compiler: GCC
# BLAS:     ATLAS
#
###

# PLASMA install directory (./install by default)
#prefix      = /scratch/jvlima/install/plasma-2.4.5

# Compilers, linker/loader, archiver, options
CC          = gcc
FC          = gfortran
LOADER      = $(FC)

ARCH        = ar
ARCHFLAGS   = cr
RANLIB      = ranlib

CFLAGS      = -O2 -DADD_
FFLAGS      = -O2
LDFLAGS     =

# Enable/Disable F90 compilation
PLASMA_F90 = 0

# BLAS library
LIBBLAS     = -L/scratch/jvlima/install/atlas3.9.69/lib -llapack -lcblas -lf77blas -latlas

# CBLAS library
LIBCBLAS    = -L/scratch/jvlima/install/atlas3.9.69/lib -llapack -lcblas -lf77blas -latlas 

# LAPACK and TMG library
LIBLAPACK   = -L/scratch/jvlima/install/atlas3.9.69/lib -llapack -latlas  $(LAPACK_LDFLAGS)
# /scratch/jvlima/install/lapack-3.4.0/libtmglib.a

# LAPACK C Wrapper library
INCCLAPACK  = -I/home/jvlima/install/lapacke/include
LIBCLAPACK  = -L/home/jvlima/install/lapacke/lib -llapacke
