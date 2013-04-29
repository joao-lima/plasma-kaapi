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
NVCC	    = nvcc
FC          = gfortran
LOADER      = $(FC)

ARCH        = ar
ARCHFLAGS   = cr
RANLIB      = ranlib

# CUDA config

CFLAGS	    = -O2 -DADD_ -g -DPLASMA_CUDA $(CUDA_CFLAGS)
FFLAGS      = -O2
LDFLAGS     = $(CUDA_LDFLAGS) -lcublas

#CFLAGS += -DCONFIG_VERBOSE
 
# Enable/Disable F90 compilation
PLASMA_F90 = 0

# Blas Library
LIBBLAS     = -L/Users/jvlima/install/atlas/lib -llapack -lcblas -lf77blas -latlas
# CBlas library
LIBCBLAS    = -L/Users/jvlima/install/atlas/lib -llapack -lcblas -lf77blas -latlas
# lapack and tmg library (lapack is included in acml)
LIBLAPACK   = /Users/jvlima/install/lapack/libtmglib.a /Users/jvlima/install/lapack/liblapack.a 
INCCLAPACK  = -I/Users/jvlima/install/lapacke/include
LIBCLAPACK  = -L/Users/jvlima/install/lapacke/lib -llapacke 

# XKaapi install
XKAAPIDIR = $(HOME)/install/xkaapi/default

