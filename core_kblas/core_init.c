
#include <stdio.h>

#include "core_kblas.h"

void core_cublas_init(void)
{
  const char* arch = getenv("QUARK_ARCH_GPU_ONLY");
  const char* prio = getenv("QUARK_PRIORITY");
  uint8_t quark_arch = QUARK_ARCH_DEFAULT;
  uint8_t quark_prio = 2;
  
  if(arch != 0)
  {
    uint8_t arch_set = atoi(arch);
    if(arch_set != 0)
    {
      quark_arch = QUARK_ARCH_GPU_ONLY;
    }
  }
  
  if(prio != 0)
  {
    quark_prio = atoi(prio);
  }
  
  QUARK_Task_Set_Function_Params("CORE_sgemm_quark", CORE_sgemm_quark, NULL, CORE_sgemm_quark_cublas,
                                 quark_arch,
                                 quark_prio);
  
  QUARK_Task_Set_Function_Params("CORE_dgemm_quark", CORE_dgemm_quark, NULL, CORE_dgemm_quark_cublas,
                                 quark_arch,
                                 quark_prio);
  
  QUARK_Task_Set_Function_Params("CORE_dgemm_f2_quark", CORE_dgemm_f2_quark, NULL, CORE_dgemm_f2_quark_cublas,
                                 quark_arch,
                                 quark_prio);
  
  QUARK_Task_Set_Function_Params("CORE_dtrsm_quark", CORE_dtrsm_quark, NULL, CORE_dtrsm_quark_cublas,
                                 quark_arch,
                                 quark_prio);
  
  QUARK_Task_Set_Function_Params("CORE_dsyrk_quark", CORE_dsyrk_quark, NULL, CORE_dsyrk_quark_cublas,
                                 quark_arch,
                                 quark_prio);
  
//  QUARK_Task_Set_Function_Params("CORE_dpotrf_quark", CORE_dpotrf_quark, CORE_dpotrf_quark_cpu, NULL,
  QUARK_Task_Set_Function_Params("CORE_dpotrf_quark", CORE_dpotrf_quark, NULL, NULL,
                                 QUARK_ARCH_CPU_ONLY,
                                 quark_prio);
  
  /* LU */  
  QUARK_Task_Set_Function_Params("CORE_dgetrf_incpiv_quark", CORE_dgetrf_incpiv_quark, NULL, NULL,
                                 QUARK_ARCH_CPU_ONLY,
                                 quark_prio);

  /* LU */
  QUARK_Task_Set_Function_Params("CORE_dtstrf_quark", CORE_dtstrf_quark, NULL, CORE_dtstrf_quark_cublas,
                                 QUARK_ARCH_CPU_ONLY,
                                 quark_prio);

  /* LU */
  QUARK_Task_Set_Function_Params("CORE_dssssm_quark", CORE_dssssm_quark, NULL, CORE_dssssm_quark_cublas,
                                 quark_arch,
                                 quark_prio);
  
  /* LU */
  QUARK_Task_Set_Function_Params("CORE_dgessm_quark", CORE_dgessm_quark, NULL, CORE_dgessm_quark_cublas,
                                 quark_arch,
                                 quark_prio);

  QUARK_Task_Set_Function_Params("CORE_dgeqrt_quark", CORE_dgeqrt_quark, NULL, NULL,
                                 QUARK_ARCH_CPU_ONLY,
                                 quark_prio);

  /* QR */
  QUARK_Task_Set_Function_Params("CORE_dormqr_quark", CORE_dormqr_quark, NULL, CORE_dormqr_quark_cublas,
                                 quark_arch,
                                 quark_prio);

  /* QR */
  QUARK_Task_Set_Function_Params("CORE_dtsqrt_quark", CORE_dtsqrt_quark, NULL, NULL,
                                 QUARK_ARCH_CPU_ONLY,
                                 quark_prio);

  /* QR */
  QUARK_Task_Set_Function_Params("CORE_dtsmqr_quark", CORE_dtsmqr_quark, NULL, CORE_dtsmqr_quark_cublas,
                                 quark_arch,
                                 quark_prio);
  
  /* Tile CAQR */
  QUARK_Task_Set_Function_Params("CORE_dttqrt_quark", CORE_dttqrt_quark, NULL, NULL,
                                 QUARK_ARCH_CPU_ONLY,
                                 quark_prio);
  
  /* Tile CAQR */
  QUARK_Task_Set_Function_Params("CORE_dttmqr_quark", CORE_dttmqr_quark, NULL, NULL,
                                 QUARK_ARCH_CPU_ONLY,
                                 quark_prio);
}
