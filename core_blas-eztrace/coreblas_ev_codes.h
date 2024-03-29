/**
 *
 * @file coreblas_ev_codes.h
 *
 *  PLASMA core_blas tracing kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 **/
#ifndef EVCODES_COREBLAS_H
#define EVCODES_COREBLAS_H

#define COREBLAS_EVENTS_ID     0x03
#define COREBLAS_PREFIX       (COREBLAS_EVENTS_ID << NB_BITS_EVENTS)

#define FUT_COREBLAS_STOP     (COREBLAS_PREFIX | 0x0000)
#define FUT_COREBLAS_TASK     (COREBLAS_PREFIX | 0x1000)
#define FUT_COREBLAS_TASKW    (COREBLAS_PREFIX | 0x1001)

/* Level 3 Blas */
#define FUT_COREBLAS_GEMM     (COREBLAS_PREFIX | 0x0001)
#define FUT_COREBLAS_HERK     (COREBLAS_PREFIX | 0x0002)
#define FUT_COREBLAS_SYRK     (COREBLAS_PREFIX | 0x0003)
#define FUT_COREBLAS_HEMM     (COREBLAS_PREFIX | 0x0005)
#define FUT_COREBLAS_SYMM     (COREBLAS_PREFIX | 0x0006)
#define FUT_COREBLAS_TRMM     (COREBLAS_PREFIX | 0x0007)
#define FUT_COREBLAS_TRSM     (COREBLAS_PREFIX | 0x0008)
#define FUT_COREBLAS_HER2K    (COREBLAS_PREFIX | 0x0009)
#define FUT_COREBLAS_SYR2K    (COREBLAS_PREFIX | 0x000a)

/* Level 2 Blas */
#define FUT_COREBLAS_GEMV     (COREBLAS_PREFIX | 0x0010)
#define FUT_COREBLAS_GBMV     (COREBLAS_PREFIX | 0x0011)
#define FUT_COREBLAS_HEMV     (COREBLAS_PREFIX | 0x0012)
#define FUT_COREBLAS_HBMV     (COREBLAS_PREFIX | 0x0013)
#define FUT_COREBLAS_HPMV     (COREBLAS_PREFIX | 0x0014)
#define FUT_COREBLAS_SYMV     (COREBLAS_PREFIX | 0x0015)
#define FUT_COREBLAS_SBMV     (COREBLAS_PREFIX | 0x0016)
#define FUT_COREBLAS_SPMV     (COREBLAS_PREFIX | 0x0017)
#define FUT_COREBLAS_TRMV     (COREBLAS_PREFIX | 0x0018)
#define FUT_COREBLAS_TBMV     (COREBLAS_PREFIX | 0x0019)
#define FUT_COREBLAS_TPMV     (COREBLAS_PREFIX | 0x001a)
#define FUT_COREBLAS_TRSV     (COREBLAS_PREFIX | 0x001b)
#define FUT_COREBLAS_TBSV     (COREBLAS_PREFIX | 0x001c)
#define FUT_COREBLAS_TPSV     (COREBLAS_PREFIX | 0x001d)
#define FUT_COREBLAS_GER      (COREBLAS_PREFIX | 0x001e)
#define FUT_COREBLAS_GERU     (COREBLAS_PREFIX | 0x001f)
#define FUT_COREBLAS_GERC     (COREBLAS_PREFIX | 0x0020)
#define FUT_COREBLAS_HER      (COREBLAS_PREFIX | 0x0021)
#define FUT_COREBLAS_HPR      (COREBLAS_PREFIX | 0x0022)
#define FUT_COREBLAS_HER2     (COREBLAS_PREFIX | 0x0023)
#define FUT_COREBLAS_HPR2     (COREBLAS_PREFIX | 0x0024)
#define FUT_COREBLAS_SYR      (COREBLAS_PREFIX | 0x0025)
#define FUT_COREBLAS_SPR      (COREBLAS_PREFIX | 0x0026)
#define FUT_COREBLAS_SYR2     (COREBLAS_PREFIX | 0x0027)
#define FUT_COREBLAS_SPR2     (COREBLAS_PREFIX | 0x0028)

/* Level 1 BLAS */
#define FUT_COREBLAS_ROTG     (COREBLAS_PREFIX | 0x0030)
#define FUT_COREBLAS_ROTMG    (COREBLAS_PREFIX | 0x0031)
#define FUT_COREBLAS_ROT      (COREBLAS_PREFIX | 0x0032)
#define FUT_COREBLAS_ROTM     (COREBLAS_PREFIX | 0x0033)
#define FUT_COREBLAS_SWAP     (COREBLAS_PREFIX | 0x0034)
#define FUT_COREBLAS_SCAL     (COREBLAS_PREFIX | 0x0035)
#define FUT_COREBLAS_COPY     (COREBLAS_PREFIX | 0x0036)
#define FUT_COREBLAS_AXPY     (COREBLAS_PREFIX | 0x0037)
#define FUT_COREBLAS_DOT      (COREBLAS_PREFIX | 0x0038)
#define FUT_COREBLAS_DOTU     (COREBLAS_PREFIX | 0x0039)
#define FUT_COREBLAS_DOTC     (COREBLAS_PREFIX | 0x003a)
#define FUT_COREBLAS_xDOT     (COREBLAS_PREFIX | 0x003b)
#define FUT_COREBLAS_NRM2     (COREBLAS_PREFIX | 0x003c)
#define FUT_COREBLAS_ASUM     (COREBLAS_PREFIX | 0x003d)
#define FUT_COREBLAS_AMAX     (COREBLAS_PREFIX | 0x003e)

/* Lapack */
#define FUT_COREBLAS_LACPY    (COREBLAS_PREFIX | 0x0050)
#define FUT_COREBLAS_LANGE    (COREBLAS_PREFIX | 0x0051)
#define FUT_COREBLAS_LANHE    (COREBLAS_PREFIX | 0x0052)
#define FUT_COREBLAS_LANSY    (COREBLAS_PREFIX | 0x0053)
#define FUT_COREBLAS_LARFB    (COREBLAS_PREFIX | 0x0054)
#define FUT_COREBLAS_LARFT    (COREBLAS_PREFIX | 0x0055)
#define FUT_COREBLAS_LASWP    (COREBLAS_PREFIX | 0x0056)
#define FUT_COREBLAS_LAUUM    (COREBLAS_PREFIX | 0x0057)
#define FUT_COREBLAS_POTRF    (COREBLAS_PREFIX | 0x0058)
#define FUT_COREBLAS_TRTRI    (COREBLAS_PREFIX | 0x0059)
#define FUT_COREBLAS_LASET    (COREBLAS_PREFIX | 0x0060)


/* PLASMA coreblas */
#define FUT_COREBLAS_GELQT    (COREBLAS_PREFIX | 0x0101)
#define FUT_COREBLAS_GEQRT    (COREBLAS_PREFIX | 0x0102)
#define FUT_COREBLAS_GESSM    (COREBLAS_PREFIX | 0x0103)
#define FUT_COREBLAS_GETRF    (COREBLAS_PREFIX | 0x0105)
#define FUT_COREBLAS_GETRO    (COREBLAS_PREFIX | 0x0106)
#define FUT_COREBLAS_SSSSM    (COREBLAS_PREFIX | 0x0107)
#define FUT_COREBLAS_TITRO    (COREBLAS_PREFIX | 0x0108)
#define FUT_COREBLAS_TRBMM    (COREBLAS_PREFIX | 0x0109)
#define FUT_COREBLAS_TRGMM    (COREBLAS_PREFIX | 0x010a)
#define FUT_COREBLAS_TSLQT    (COREBLAS_PREFIX | 0x010b)
#define FUT_COREBLAS_TSMLQ    (COREBLAS_PREFIX | 0x010c)
#define FUT_COREBLAS_TSMQR    (COREBLAS_PREFIX | 0x010d)
#define FUT_COREBLAS_TSQRT    (COREBLAS_PREFIX | 0x010e)
#define FUT_COREBLAS_TSRFB    (COREBLAS_PREFIX | 0x010f)
#define FUT_COREBLAS_TSTRF    (COREBLAS_PREFIX | 0x0110)
#define FUT_COREBLAS_TTLQT    (COREBLAS_PREFIX | 0x0111)
#define FUT_COREBLAS_TTMLQ    (COREBLAS_PREFIX | 0x0112)
#define FUT_COREBLAS_TTMQR    (COREBLAS_PREFIX | 0x0113)
#define FUT_COREBLAS_TTQRT    (COREBLAS_PREFIX | 0x0114)
#define FUT_COREBLAS_TTRFB    (COREBLAS_PREFIX | 0x0115)
#define FUT_COREBLAS_UNMLQ    (COREBLAS_PREFIX | 0x0116)
#define FUT_COREBLAS_UNMQR    (COREBLAS_PREFIX | 0x0117)
#define FUT_COREBLAS_GETRIP   (COREBLAS_PREFIX | 0x0118)
#define FUT_COREBLAS_PLGHE    (COREBLAS_PREFIX | 0x0119)
#define FUT_COREBLAS_PLGSY    (COREBLAS_PREFIX | 0x011a)
#define FUT_COREBLAS_SHIFT    (COREBLAS_PREFIX | 0x011b)
#define FUT_COREBLAS_SHIFTW   (COREBLAS_PREFIX | 0x011c)
#define FUT_COREBLAS_SWPAB    (COREBLAS_PREFIX | 0x011d)
#define FUT_COREBLAS_PLRNT    (COREBLAS_PREFIX | 0x011e)

#define FUT_COREBLAS_BRDALG   (COREBLAS_PREFIX | 0x0120)
#define FUT_COREBLAS_TRDALG   (COREBLAS_PREFIX | 0x0121)
#define FUT_COREBLAS_HEGST    (COREBLAS_PREFIX | 0x0122)
#define FUT_COREBLAS_SYGST    (COREBLAS_PREFIX | 0x0123)
#define FUT_COREBLAS_HERFB    (COREBLAS_PREFIX | 0x0124)
#define FUT_COREBLAS_SYRFB    (COREBLAS_PREFIX | 0x0125)

#define COREBLAS_MASK_EVENTS   0x0fff
#define COREBLAS_NBMAX_EVENTS  0x0126

#endif /* COREBLAS_CODES_H */
