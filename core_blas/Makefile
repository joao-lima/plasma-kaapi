###
#
# @file Makefile
#
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.4.2
# @date 2010-11-15
#
###

PLASMA_DIR = ..
include ../Makefile.internal

# Add the headers from src directory
INC := -I../control/ $(INC)

ZSRC =  core_zgeqrt.c core_ztsqrt.c core_ztsmqr.c core_ztsrfb.c \
        core_zgelqt.c core_ztslqt.c core_ztsmlq.c               \
                      core_zttqrt.c core_zttmqr.c core_zttrfb.c \
                      core_zttlqt.c core_zttmlq.c               \
	core_zgetrf_incpiv.c core_ztstrf.c core_zgessm.c core_zssssm.c \
	core_zpotrf.c core_ztrtri.c core_zlauum.c               \
	core_zgemm.c  core_zsymm.c  core_zsyrk.c core_ztrsm.c   \
	core_ztrmm.c  core_zaxpy.c  core_zsyr2k.c               \
	core_zlacpy.c core_zlange.c core_zlansy.c               \
	core_dzasum.c core_zunmqr.c core_zunmlq.c               \
	core_zhemm.c  core_zherk.c  core_zlanhe.c core_zher2k.c \
	core_zlag2c.c                                           \
        core_zswpab.c core_zgetrip.c core_zshift.c		\
	core_zplrnt.c core_zplghe.c core_zplgsy.c		\
	core_zlaswp.c core_zgetrf.c 				\
	core_zgetrf_rectil.c core_zgetrf_reclap.c               \
	core_zlaset.c core_zlaset2.c 				\
	core_zherfb.c core_zhegst.c 				\
	core_zlarfx_tbrd.c core_ztrdalg.c core_zbrdalg.c	\
	core_zhbelr.c core_zhbrce.c core_zhblrx.c              	\
	core_zgbelr.c core_zgbrce.c core_zgblrx.c 		\
	core_ztsmqr_hetra1.c core_ztsmqr_corner.c		\
	core_ztsmlq_hetra1.c core_ztsmlq_corner.c

CSRC = core_cgeqrt.c core_ctsqrt.c core_ctsmqr.c core_ctsrfb.c core_cgelqt.c core_ctslqt.c core_ctsmlq.c core_cttqrt.c core_cttmqr.c core_cttrfb.c core_cttlqt.c core_cttmlq.c core_cgetrf_incpiv.c core_ctstrf.c core_cgessm.c core_cssssm.c core_cpotrf.c core_ctrtri.c core_clauum.c core_cgemm.c core_csymm.c core_csyrk.c core_ctrsm.c core_ctrmm.c core_caxpy.c core_csyr2k.c core_clacpy.c core_clange.c core_clansy.c core_scasum.c core_cunmqr.c core_cunmlq.c core_chemm.c core_cherk.c core_clanhe.c core_cher2k.c core_cswpab.c core_cgetrip.c core_cshift.c core_cplrnt.c core_cplghe.c core_cplgsy.c core_claswp.c core_cgetrf.c core_cgetrf_rectil.c core_cgetrf_reclap.c core_claset.c core_claset2.c core_cherfb.c core_chegst.c core_clarfx_tbrd.c core_ctrdalg.c core_cbrdalg.c core_chbelr.c core_chbrce.c core_chblrx.c core_cgbelr.c core_cgbrce.c core_cgblrx.c core_ctsmqr_hetra1.c core_ctsmqr_corner.c core_ctsmlq_hetra1.c core_ctsmlq_corner.c
DSRC = core_dgeqrt.c core_dtsqrt.c core_dtsmqr.c core_dtsrfb.c core_dgelqt.c core_dtslqt.c core_dtsmlq.c core_dttqrt.c core_dttmqr.c core_dttrfb.c core_dttlqt.c core_dttmlq.c core_dgetrf_incpiv.c core_dtstrf.c core_dgessm.c core_dssssm.c core_dpotrf.c core_dtrtri.c core_dlauum.c core_dgemm.c core_dsymm.c core_dsyrk.c core_dtrsm.c core_dtrmm.c core_daxpy.c core_dsyr2k.c core_dlacpy.c core_dlange.c core_dlansy.c core_dasum.c core_dormqr.c core_dormlq.c core_dlag2s.c core_dswpab.c core_dgetrip.c core_dshift.c core_dplrnt.c core_dplgsy.c core_dlaswp.c core_dgetrf.c core_dgetrf_rectil.c core_dgetrf_reclap.c core_dlaset.c core_dlaset2.c core_dsyrfb.c core_dsygst.c core_dlarfx_tbrd.c core_dtrdalg.c core_dbrdalg.c core_dhbelr.c core_dhbrce.c core_dhblrx.c core_dgbelr.c core_dgbrce.c core_dgblrx.c core_dtsmqr_sytra1.c core_dtsmqr_corner.c core_dtsmlq_sytra1.c core_dtsmlq_corner.c
SSRC = core_sgeqrt.c core_stsqrt.c core_stsmqr.c core_stsrfb.c core_sgelqt.c core_stslqt.c core_stsmlq.c core_sttqrt.c core_sttmqr.c core_sttrfb.c core_sttlqt.c core_sttmlq.c core_sgetrf_incpiv.c core_ststrf.c core_sgessm.c core_sssssm.c core_spotrf.c core_strtri.c core_slauum.c core_sgemm.c core_ssymm.c core_ssyrk.c core_strsm.c core_strmm.c core_saxpy.c core_ssyr2k.c core_slacpy.c core_slange.c core_slansy.c core_sasum.c core_sormqr.c core_sormlq.c core_sswpab.c core_sgetrip.c core_sshift.c core_splrnt.c core_splgsy.c core_slaswp.c core_sgetrf.c core_sgetrf_rectil.c core_sgetrf_reclap.c core_slaset.c core_slaset2.c core_ssyrfb.c core_ssygst.c core_slarfx_tbrd.c core_strdalg.c core_sbrdalg.c core_shbelr.c core_shbrce.c core_shblrx.c core_sgbelr.c core_sgbrce.c core_sgblrx.c core_stsmqr_sytra1.c core_stsmqr_corner.c core_stsmlq_sytra1.c core_stsmlq_corner.c

ALLSRC = core_alloc.c global.c	\
	 $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJ = $(ALLSRC:.c=.o)

all: $(LIBCOREBLAS)

$(LIBCOREBLAS): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

clean:
	rm -f *.o *~

cleanall: clean
	rm -f *.a

.c.o :
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

.PHONY: all clean cleanall
