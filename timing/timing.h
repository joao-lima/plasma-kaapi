#ifndef TIMING_H
#define TIMING_H

typedef double real_Double_t;

enum iparam_timing {
  TIMING_CHECK,
  TIMING_WARMUP,
  TIMING_NITER,
  TIMING_N,
  TIMING_M,
  TIMING_NB,
  /*TIMING_MB,*/
  TIMING_IB,
  TIMING_NRHS,
  TIMING_THRDNBR,
  TIMING_THRDNBR_SUBGRP,
  TIMING_SCHEDULER,
  TIMING_AUTOTUNING,
  TIMING_INPUTFMT,
  TIMING_OUTPUTFMT,
  TIMING_INBPARAM
};

enum dparam_timing {
  TIMING_TIME,
  TIMING_ANORM,
  TIMING_BNORM,
  TIMING_XNORM,
  TIMING_RES,
  TIMING_DNBPARAM
};

#endif /* TIMING_H */
