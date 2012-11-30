###
#
# @file Makefile
#
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.4.2
# @author Bilel Hadri
# @date 2010-11-15
#
###

PLASMA_DIR = ..
include ../Makefile.internal

ZSRC =  example_zgelqf.c  	\
	example_zgels.c 	\
	example_zgeqrs.c	\
	example_zgetrs.c  	\
	example_zpotrf.c  	\
	example_ztrsm.c   	\
        example_zgelqs.c  	\
	example_zgeqrf.c  	\
	example_zgesv.c   	\
	example_zposv.c   	\
	example_zpotrs.c  	\
	example_ztrsmpl.c	\
	example_zunmqr.c 	\
	example_zunmlq.c	\
	\
	example_zposv_f.f	\
	example_zgesv_f.f 	\
	example_zgels_f.f

CSRC = example_cgelqf.c example_cgels.c example_cgeqrs.c example_cgetrs.c example_cpotrf.c example_ctrsm.c example_cgelqs.c example_cgeqrf.c example_cgesv.c example_cposv.c example_cpotrs.c example_ctrsmpl.c example_cunmqr.c example_cunmlq.c example_cposv_f.f example_cgesv_f.f example_cgels_f.f
DSRC = example_dgelqf.c example_dgels.c example_dgeqrs.c example_dgetrs.c example_dpotrf.c example_dtrsm.c example_dgelqs.c example_dgeqrf.c example_dgesv.c example_dposv.c example_dpotrs.c example_dtrsmpl.c example_dormqr.c example_dormlq.c example_dposv_f.f example_dgesv_f.f example_dgels_f.f
SSRC = example_sgelqf.c example_sgels.c example_sgeqrs.c example_sgetrs.c example_spotrf.c example_strsm.c example_sgelqs.c example_sgeqrf.c example_sgesv.c example_sposv.c example_spotrs.c example_strsmpl.c example_sormqr.c example_sormlq.c example_sposv_f.f example_sgesv_f.f example_sgels_f.f

ALLSRC = $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJT= $(ALLSRC:.c=.o)
ALLOBJ = $(ALLOBJT:.f=.o)

EXE = $(ALLOBJ:.o=)

all: $(EXE)

$(EXE) : $(LIBCOREBLAS) $(LIBPLASMA) $(QUARKDIR)/libquark.a 

clean:
	rm -f $(ALLOBJ) 

cleanall: clean
	rm -f $(EXE)

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

.f.o:
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

example_% : example_%.o
	$(LOADER) $(LDFLAGS) $< -o $@ $(LIB)

.PHONY: clean cleanall
