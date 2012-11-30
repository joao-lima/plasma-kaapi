###
#
# @file Makefile
#
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.4.2
# @author Bilel Hadri
# @author Hatem Ltaief
# @date 2010-11-15
#
###

PLASMA_DIR = ..
include ../Makefile.internal

ZHDR =  testing_zmain.h
CHDR = testing_cmain.h
DHDR = testing_dmain.h
SHDR = testing_smain.h
ALLHDR= $(ZHDR) $(CHDR) $(DHDR) $(SHDR)

ZSRC  = testing_zmain.c 	\
	testing_zgesv.c 	\
	testing_zgesvd.c 	\
	testing_zgels.c 	\
	testing_zposv.c 	\
	testing_zgemm.c 	\
	testing_zheev.c 	\
	testing_zhegst.c 	\
	testing_zhegv.c 	\
	testing_zsyrk.c 	\
	testing_zsyr2k.c 	\
	testing_zsymm.c 	\
	testing_ztrsm.c 	\
	testing_ztrmm.c 	\
	testing_zlange.c	\
	testing_zherk.c		\
	testing_zher2k.c	\
	testing_zhemm.c		\
	testing_zgecfi.c	\
	testing_zgetmi.c	\
	testing_zcposv.c 	\
	testing_zcgesv.c 	\
	testing_zcungesv.c 

CSRC = testing_cmain.c testing_cgesv.c testing_cgesvd.c testing_cgels.c testing_cposv.c testing_cgemm.c testing_cheev.c testing_chegst.c testing_chegv.c testing_csyrk.c testing_csyr2k.c testing_csymm.c testing_ctrsm.c testing_ctrmm.c testing_clange.c testing_cherk.c testing_cher2k.c testing_chemm.c testing_cgecfi.c testing_cgetmi.c
DSRC = testing_dmain.c testing_dgesv.c testing_dgesvd.c testing_dgels.c testing_dposv.c testing_dgemm.c testing_dsyev.c testing_dsygst.c testing_dsygv.c testing_dsyrk.c testing_dsyr2k.c testing_dsymm.c testing_dtrsm.c testing_dtrmm.c testing_dlange.c testing_dgecfi.c testing_dgetmi.c testing_dsposv.c testing_dsgesv.c testing_dsungesv.c
SSRC = testing_smain.c testing_sgesv.c testing_sgesvd.c testing_sgels.c testing_sposv.c testing_sgemm.c testing_ssyev.c testing_ssygst.c testing_ssygv.c testing_ssyrk.c testing_ssyr2k.c testing_ssymm.c testing_strsm.c testing_strmm.c testing_slange.c testing_sgecfi.c testing_sgetmi.c

ZOBJ = $(ZSRC:.c=.o)
COBJ = $(CSRC:.c=.o)
DOBJ = $(DSRC:.c=.o)
SOBJ = $(SSRC:.c=.o)

ALLSRC = $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJ = $(ALLSRC:.c=.o)

EXE = ztesting ctesting dtesting stesting

all: $(EXE)

$(EXE) : $(LIBCOREBLAS) $(LIBPLASMA) $(QUARKDIR)/libquark.a 

clean:
	rm -f $(ALLOBJ) *~

cleanall: clean
	rm -f $(EXE)

testing_z%.o: testing_z%.c testing_zmain.h
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

testing_c%.o: testing_c%.c testing_cmain.h
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

testing_d%.o: testing_d%.c testing_dmain.h
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

testing_s%.o: testing_s%.c testing_smain.h
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

.f.o:
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

ztesting : testing_zmain.h $(ZOBJ)
	$(LOADER) $(LDFLAGS) $(ZOBJ) -o $@ $(LIB)

ctesting : testing_cmain.h $(COBJ)
	$(LOADER) $(LDFLAGS) $(COBJ) -o $@ $(LIB)

dtesting : testing_dmain.h $(DOBJ)
	$(LOADER) $(LDFLAGS) $(DOBJ) -o $@ $(LIB)

stesting : testing_smain.h $(SOBJ)
	$(LOADER) $(LDFLAGS) $(SOBJ) -o $@ $(LIB)

.PHONY: clean cleanall
