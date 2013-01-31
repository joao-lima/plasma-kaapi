###
#
# @file Makefile
#
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.4.2
# @author Mathieu Faverge
# @date 2010-11-15
#
###

XKAAPIDIR = $(HOME)/install/xkaapi/default
QUARKDIR = $(HOME)/install/xkaapi/default/lib
INCQUARK  = -I$(XKAAPIDIR)/include
LIBQUARK  = -L$(XKAAPIDIR)/lib -lquark -lkaapi

PLASMA_DIR = ..
include ../Makefile.internal

CFLAGS  += -g -DCONFIG_USE_CUDA -I../control/ $(INC) -I$(HOME)/xkaapi/api/quark/core_cublas $(CUDA_CFLAGS)
LDFLAGS +=  -lcore_cublas $(LIB) -L$(HOME)/xkaapi/api/quark/core_cublas $(CUDA_LDFLAGS) -lcublas

#CFLAGS  += -g -I../control/ $(INC)
#LDFLAGS +=  $(LIB) 

ZHDR =  zauxiliary.h zauxiliary.c
CHDR = cauxiliary.h cauxiliary.c
DHDR = dauxiliary.h dauxiliary.c
SHDR = sauxiliary.h sauxiliary.c

#timing.c is in headers because it's a template for every timing routine
ALLHDR= auxiliary.h auxiliary.c timing.h timing.c \
	$(ZHDR) $(CHDR) $(DHDR) $(SHDR)

ZSRC  = time_zgemm.c 		\
	time_zcposv.c 		\
	time_zposv.c 		\
	time_zpotrf.c 		\
	time_zcgesv.c 		\
	time_zgesv.c 		\
	time_zgetrf.c 		\
	time_zgesv_incpiv.c 	\
	time_zgetrf_incpiv.c	\
	time_zgels.c 		\
	time_zgeqrf.c 		\
	time_zgemm_tile.c 	\
	time_zcposv_tile.c 	\
	time_zposv_tile.c 	\
	time_zpotrf_tile.c 	\
	time_zpotri_tile.c	\
	time_zcgesv_tile.c 	\
	time_zgesv_tile.c 	\
	time_zgetrf_tile.c 	\
	time_zgesv_incpiv_tile.c 	\
	time_zgetrf_incpiv_tile.c 	\
	time_zgels_tile.c 	\
	time_zgeqrf_tile.c      \
	time_zgeqrfrh_tile.c    \
	time_zgelqfrh_tile.c    \
        time_zgecfi.c           \
	time_zheev_tile.c	\
	time_zhegv_tile.c	\
	time_zgesvd_tile.c	\
        time_zlapack2tile.c     \
	time_zgetrf_reclap.c	\
	time_zgetrf_rectil.c

CSRC = time_cgemm.c time_cposv.c time_cpotrf.c time_cgesv.c time_cgetrf.c time_cgesv_incpiv.c time_cgetrf_incpiv.c time_cgels.c time_cgeqrf.c time_cgemm_tile.c time_cposv_tile.c time_cpotrf_tile.c time_cpotri_tile.c time_cgesv_tile.c time_cgetrf_tile.c time_cgesv_incpiv_tile.c time_cgetrf_incpiv_tile.c time_cgels_tile.c time_cgeqrf_tile.c time_cgeqrfrh_tile.c time_cgelqfrh_tile.c time_cgecfi.c time_cheev_tile.c time_chegv_tile.c time_cgesvd_tile.c time_clapack2tile.c time_cgetrf_reclap.c time_cgetrf_rectil.c
DSRC = time_dgemm.c time_dsposv.c time_dposv.c time_dpotrf.c time_dsgesv.c time_dgesv.c time_dgetrf.c time_dgesv_incpiv.c time_dgetrf_incpiv.c time_dgels.c time_dgeqrf.c time_dgemm_tile.c time_dsposv_tile.c time_dposv_tile.c time_dpotrf_tile.c time_dpotri_tile.c time_dsgesv_tile.c time_dgesv_tile.c time_dgetrf_tile.c time_dgesv_incpiv_tile.c time_dgetrf_incpiv_tile.c time_dgels_tile.c time_dgeqrf_tile.c time_dgeqrfrh_tile.c time_dgelqfrh_tile.c time_dgecfi.c time_dsyev_tile.c time_dsygv_tile.c time_dgesvd_tile.c time_dlapack2tile.c time_dgetrf_reclap.c time_dgetrf_rectil.c
SSRC = time_sgemm.c time_sposv.c time_spotrf.c time_sgesv.c time_sgetrf.c time_sgesv_incpiv.c time_sgetrf_incpiv.c time_sgels.c time_sgeqrf.c time_sgemm_tile.c time_sposv_tile.c time_spotrf_tile.c time_spotri_tile.c time_sgesv_tile.c time_sgetrf_tile.c time_sgesv_incpiv_tile.c time_sgetrf_incpiv_tile.c time_sgels_tile.c time_sgeqrf_tile.c time_sgeqrfrh_tile.c time_sgelqfrh_tile.c time_sgecfi.c time_ssyev_tile.c time_ssygv_tile.c time_sgesvd_tile.c time_slapack2tile.c time_sgetrf_reclap.c time_sgetrf_rectil.c

ALLSRC = $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJ = $(ALLSRC:.c=.o)
EXE = $(ALLOBJ:.o=)

all: $(EXE)

$(EXE) : $(LIBCOREBLAS) $(LIBPLASMA) $(QUARKDIR)/libquark.a 

clean:
	rm -f *.o *~

cleanall: clean
	rm -f $(EXE)

%.c: $(ALLHDR)

time_%.o : time_%.c $(ALLHDR)
	$(CC) $(CFLAGS) -c $< -o $@

%auxiliary.o : %auxiliary.c $(ALLHDR)
	$(CC) $(CFLAGS) -c $< -o $@

time_s% : time_s%.o sauxiliary.o
	$(LOADER) $^ -o $@ $(LDFLAGS)
time_d% : time_d%.o dauxiliary.o
	$(LOADER) $^ -o $@ $(LDFLAGS)
time_c% : time_c%.o cauxiliary.o
	$(LOADER) $^ -o $@ $(LDFLAGS)
time_z% : time_z%.o zauxiliary.o
	$(LOADER) $^ -o $@ $(LDFLAGS)

time_zlapack2tile time_clapack2tile time_dlapack2tile time_slapack2tile : auxiliary.o
time_zgecfi time_cgecfi time_dgecfi time_sgecfi : auxiliary.o

png: time_dgesv time_zgesv
	./time_dgesv --gnuplot-png | tee timeplot.gpt ; gnuplot timeplot.gpt

.PHONY: all clean
