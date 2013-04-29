
PLASMA_DIR = ..
include ../Makefile.internal

# Add the headers from src directory
INC := -I../control/ $(INC)

#CFLAGS = -g -DKAAPI_USE_CUDA -I.. -I$(HOME)/install/xkaapi/default/include \
#	  $(CBLAS_CFLAGS) \
#	  $(LAPACKE_CFLAGS) \
#	  -I$(HOME)/plasma_2.4.2/include \
#	  -I$(HOME)/plasma_2.4.2 \
#	  $(CUDA_CFLAGS) 
##CFLAGS += -DCONFIG_VERBOSE

SRC = \
    core_dpotrf.c \
    core_dgemm.c \
    core_sgemm.c \
    core_dsyrk.c \
    core_dtrsm.c \
    core_dssssm.c \
    core_dgessm.c \
    core_dtsmqr.c \
    core_dtsrfb.c \
    core_dormqr.c \
    core_dtstrf.c \
    core_dlarfb_gemm.c \
    core_init.c


ALLOBJ := $(SRC:.c=.o) magmablas.o core_dtstrf_kernel.o

all: $(LIBCOREKBLAS)

magmablas.o: magmablas.cu
	$(NVCC) --machine 64 -arch=sm_20 $(CFLAGS) -c $< -o $@

core_dtstrf_kernel.o: core_dtstrf_kernel.cu
	$(NVCC) --machine 64 -arch=sm_20 $(CFLAGS) -c $< -o $@

$(LIBCOREKBLAS): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

clean:
	rm -f *.o *~

cleanall: clean
	rm -f *.a

.c.o :
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

.PHONY: all clean cleanall
