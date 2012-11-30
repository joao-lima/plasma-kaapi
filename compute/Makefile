###
#
# @file Makefile
#
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.4.2
# @author Jakub Kurzak
# @date 2010-11-15
#
###

# Overwritten in make.inc
PLASMA_DIR = ..
include ../Makefile.internal

INC := -I$(PLASMA_DIR)/control/ $(INC)

ZHDR =  zgecfi2.h
CHDR = cgecfi2.h
DHDR = dgecfi2.h
SHDR = sgecfi2.h

ZSRC =  zgels.c  zgelqf.c zgelqs.c zgeqrf.c zgeqrs.c 	\
        zgesv_incpiv.c  zgetrf_incpiv.c zgetrs_incpiv.c	\
        zgesv.c  zgetrf.c zgetrs.c zlaswp.c		\
        zposv.c  zpotrf.c zpotrs.c             		\
	ztrsm.c  zgemm.c  zsymm.c zsyrk.c ztrsmpl.c  	\
	zpotri.c ztrtri.c zlauum.c ztrmm.c zsyr2k.c	\
	zlange.c zlansy.c				\
	zunmqr.c zunmlq.c zungqr.c zunglq.c 		\
	zhemm.c  zherk.c  zlanhe.c zher2k.c   		\
	zcposv.c zcgesv.c zcungesv.c			\
        zgecfi2.c zgecfi.c zgetmi.c                     \
	zplrnt.c zplghe.c zplgsy.c			\
	zhegv.c zheev.c zhegst.c zhetrd.c		\
	zgesvd.c zgebrd.c	 			\
        \
	pzgelqf.c pzgeqrf.c pzpotrf.c pztrtri.c pzlauum.c  	\
        pzgelqfrh.c pzgeqrfrh.c                                 \
        pztrsm.c  pzgemm.c  pzsymm.c  pzsyrk.c  pztrsmpl.c 	\
	pzlange.c pzlansy.c pzsyr2k.c				\
	pzlacpy.c pzaxpy.c  pzgetrf_incpiv.c pztrmm.c		\
	pzunmqr.c pzunmlq.c pzungqr.c pzunglq.c                 \
	pzunmqrrh.c pzungqrrh.c				        \
	pzunmlqrh.c pzunglqrh.c				        \
        pzlanhe.c pzhemm.c pzherk.c pzher2k.c			\
	pzlag2c.c pzlaset.c pzlaset2.c                          \
        pzpack.c  pzshift.c pzgetmi2.c				\
	pzplrnt.c pzplghe.c pzplgsy.c                           \
        pzgetrf_reclap.c pzgetrf_rectil.c pzlaswp.c		\
	pzbarrier.c                                             \
	pzgerbb.c pzherbt.c pzhegst.c				\
	pzgbrdb.c pzhbrdt.c

CSRC = cgels.c cgelqf.c cgelqs.c cgeqrf.c cgeqrs.c cgesv_incpiv.c cgetrf_incpiv.c cgetrs_incpiv.c cgesv.c cgetrf.c cgetrs.c claswp.c cposv.c cpotrf.c cpotrs.c ctrsm.c cgemm.c csymm.c csyrk.c ctrsmpl.c cpotri.c ctrtri.c clauum.c ctrmm.c csyr2k.c clange.c clansy.c cunmqr.c cunmlq.c cungqr.c cunglq.c chemm.c cherk.c clanhe.c cher2k.c cgecfi2.c cgecfi.c cgetmi.c cplrnt.c cplghe.c cplgsy.c chegv.c cheev.c chegst.c chetrd.c cgesvd.c cgebrd.c pcgelqf.c pcgeqrf.c pcpotrf.c pctrtri.c pclauum.c pcgelqfrh.c pcgeqrfrh.c pctrsm.c pcgemm.c pcsymm.c pcsyrk.c pctrsmpl.c pclange.c pclansy.c pcsyr2k.c pclacpy.c pcaxpy.c pcgetrf_incpiv.c pctrmm.c pcunmqr.c pcunmlq.c pcungqr.c pcunglq.c pcunmqrrh.c pcungqrrh.c pcunmlqrh.c pcunglqrh.c pclanhe.c pchemm.c pcherk.c pcher2k.c pclaset.c pclaset2.c pcpack.c pcshift.c pcgetmi2.c pcplrnt.c pcplghe.c pcplgsy.c pcgetrf_reclap.c pcgetrf_rectil.c pclaswp.c pcbarrier.c pcgerbb.c pcherbt.c pchegst.c pcgbrdb.c pchbrdt.c
DSRC = dgels.c dgelqf.c dgelqs.c dgeqrf.c dgeqrs.c dgesv_incpiv.c dgetrf_incpiv.c dgetrs_incpiv.c dgesv.c dgetrf.c dgetrs.c dlaswp.c dposv.c dpotrf.c dpotrs.c dtrsm.c dgemm.c dsymm.c dsyrk.c dtrsmpl.c dpotri.c dtrtri.c dlauum.c dtrmm.c dsyr2k.c dlange.c dlansy.c dormqr.c dormlq.c dorgqr.c dorglq.c dsposv.c dsgesv.c dsungesv.c dgecfi2.c dgecfi.c dgetmi.c dplrnt.c dplgsy.c dsygv.c dsyev.c dsygst.c dsytrd.c dgesvd.c dgebrd.c pdgelqf.c pdgeqrf.c pdpotrf.c pdtrtri.c pdlauum.c pdgelqfrh.c pdgeqrfrh.c pdtrsm.c pdgemm.c pdsymm.c pdsyrk.c pdtrsmpl.c pdlange.c pdlansy.c pdsyr2k.c pdlacpy.c pdaxpy.c pdgetrf_incpiv.c pdtrmm.c pdormqr.c pdormlq.c pdorgqr.c pdorglq.c pdormqrrh.c pdorgqrrh.c pdormlqrh.c pdorglqrh.c pdlag2s.c pdlaset.c pdlaset2.c pdpack.c pdshift.c pdgetmi2.c pdplrnt.c pdplgsy.c pdgetrf_reclap.c pdgetrf_rectil.c pdlaswp.c pdbarrier.c pdgerbb.c pdsyrbt.c pdsygst.c pdgbrdb.c pdsbrdt.c
SSRC = sgels.c sgelqf.c sgelqs.c sgeqrf.c sgeqrs.c sgesv_incpiv.c sgetrf_incpiv.c sgetrs_incpiv.c sgesv.c sgetrf.c sgetrs.c slaswp.c sposv.c spotrf.c spotrs.c strsm.c sgemm.c ssymm.c ssyrk.c strsmpl.c spotri.c strtri.c slauum.c strmm.c ssyr2k.c slange.c slansy.c sormqr.c sormlq.c sorgqr.c sorglq.c sgecfi2.c sgecfi.c sgetmi.c splrnt.c splgsy.c ssygv.c ssyev.c ssygst.c ssytrd.c sgesvd.c sgebrd.c psgelqf.c psgeqrf.c pspotrf.c pstrtri.c pslauum.c psgelqfrh.c psgeqrfrh.c pstrsm.c psgemm.c pssymm.c pssyrk.c pstrsmpl.c pslange.c pslansy.c pssyr2k.c pslacpy.c psaxpy.c psgetrf_incpiv.c pstrmm.c psormqr.c psormlq.c psorgqr.c psorglq.c psormqrrh.c psorgqrrh.c psormlqrh.c psorglqrh.c pslaset.c pslaset2.c pspack.c psshift.c psgetmi2.c psplrnt.c psplgsy.c psgetrf_reclap.c psgetrf_rectil.c pslaswp.c psbarrier.c psgerbb.c pssyrbt.c pssygst.c psgbrdb.c pssbrdt.c

ALLSRC = $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJ = $(ALLSRC:.c=.o)

all: $(LIBPLASMA)

$(LIBPLASMA): $(HEADERS) $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

clean:
	rm -f *.o *~

cleanall: clean
	rm -f *.a

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

zgecfi2.c zgecfi.c : $(ZHDR)
cgecfi2.c cgecfi.c : $(CHDR)
dgecfi2.c dgecfi.c : $(DHDR)
sgecfi2.c sgecfi.c : $(SHDR)

.PHONY: clean cleanall
