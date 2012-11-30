/**
 *
 * @file eztrace_convert_coreblas.c
 *
 *  PLASMA core_blas tracing kernels
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * This file provides the functions to generate the trace 
 * in function of the events.
 *
 * @version 2.4.2
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 * Easy way to add new kernel:
 *   1 - Add them to the file coreblas_ev_codes.h
 *   2 - First replacement for HANDLE:
 * #define FUT_COREBLAS_\([0-9A-Z]*\)[ ]*(COREBLAS_PREFIX | 0x[0-9a-f]*) -> HANDLE(\,(downcase \1))
 *   3 - second replacement:
 * #define FUT_COREBLAS_\([0-9A-Z]*\)[ ]*(COREBLAS_PREFIX | 0x[0-9a-f]*) ->   addEntityValue("\,(downcase \1)", COREBLAS_STATE, "\,(downcase \1)", GTG_PINK    );
 *   4 - third and last replacement:
 * 
 *
 **/
#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <GTG.h>
#include <ev_codes.h>
#include <eztrace_list.h>
#include <eztrace_convert.h>
#include "coreblas_ev_codes.h"

#ifndef min
#define min( a, b ) ( (a) < (b) ? (a) : (b) )
#endif
#ifndef max
#define max( a, b ) ( (a) > (b) ? (a) : (b) )
#endif

#define COREBLAS_STATE      "ST_Thread"
#define COREBLAS_TASK_NAME  "Submitted Tasks counter"
#define COREBLAS_TASK_ALIAS "STasks"

#define COREBLAS_TASKR_NAME  "Global Ready Tasks counter"
#define COREBLAS_TASKR_ALIAS "GRTasks"

#define COREBLAS_TASKWR_NAME  "Local Ready Tasks counter"
#define COREBLAS_TASKWR_ALIAS "LRTasks"

#define COREBLAS_THREADS_MAX 4096

/*
 * Statistics
 */
typedef struct coreblas_stats_s {
    int  nb;
    double sum;
    double min;
    double max;
} coreblas_stats_t;

static coreblas_stats_t *statsarray = NULL;

typedef struct coreblas_thrdstate_s {
    unsigned int tid;
    int  active;
    double lasttime;
} coreblas_thrdstate_t;

static coreblas_thrdstate_t *thrdstate = NULL;
static int nbtrhd = 0;

#include "coreblas_string.c"

/* 
 * Case with priority and request need to be handle correctly in GTG
 * before to be included here
 */
#ifdef TRACE_BY_SEQUENCE

#define MAX_SEQUENCE 100
typedef struct sequence_s {
    uint64_t id;
    char *name;
} sequence_t;

sequence_t *seqtab;
gtg_color_t colors[20];

void sequenceInit(){
    seqtab = (sequence_t*)malloc(MAX_SEQUENCE * sizeof(sequence_t));
    memset(seqtab, 0, MAX_SEQUENCE * sizeof(sequence_t));

    colors[ 0] = GTG_RED;       
    colors[ 1] = GTG_GREEN;     
    colors[ 2] = GTG_BLUE;      
    colors[ 3] = GTG_WHITE;     
    colors[ 4] = GTG_TEAL;      
    colors[ 5] = GTG_DARKGREY;  
    colors[ 6] = GTG_YELLOW;    
    colors[ 7] = GTG_PURPLE;    
    colors[ 8] = GTG_LIGHTBROWN;
    colors[ 9] = GTG_DARKBLUE;  
    colors[10] = GTG_PINK;      
    colors[11] = GTG_DARKPINK;  
    colors[12] = GTG_SEABLUE;   
    colors[13] = GTG_KAKI;      
    colors[14] = GTG_REDBLOOD;  
    colors[15] = GTG_BROWN;     
    colors[16] = GTG_GRENAT;    
    colors[17] = GTG_ORANGE;    
    colors[18] = GTG_MAUVE;     
    colors[19] = GTG_LIGHTPINK; 

}

void sequenceDestroy(){
    int i=0;
    while( i < MAX_SEQUENCE && seqtab[i].id != 0)
    {
        free(seqtab[i].name);
        i++;
    }
    free(seqtab);
}

char *getSequence(uint64_t seq)
{
    int i=0;

    while ( (i < MAX_SEQUENCE)
            && (seqtab[i].id != 0)
            && (seqtab[i].id != seq) )
        i++;
    
    if (i < MAX_SEQUENCE)
    {
        if ( seqtab[i].id == seq )
        {
            return seqtab[i].name;
        }
        else
        {
            seqtab[i].id = seq;
            asprintf(&(seqtab[i].name), "Sequence%03d", i);
            addEntityValue(seqtab[i].name, COREBLAS_STATE, seqtab[i].name, colors[i%20] );
            return seqtab[i].name;
        }
    } else {
        fprintf(stderr, "WARNING: Too many sequences, you need to increase the limit and recompile\n");
        return "SequenceOutOfRange";
    }
}
#define HANDLE(func)							\
    void handle_coreblas_##func##_start (struct fxt_ev_64 *ev)          \
    {									\
        FUNC_NAME;                                                      \
        INIT_THREAD_ID(_threadstr);                                     \
        if ( GET_NBPARAMS(ev) > 0 ) {                                   \
            CHANGE() setState (CURRENT, COREBLAS_STATE, _threadstr, getSequence(GET_PARAM(ev, 2))); \
        } else {                                                        \
            CHANGE() setState (CURRENT, COREBLAS_STATE, _threadstr, #func); \
        }                                                               \
        free(_threadstr);                                               \
    }
#else
#define HANDLE(func)							\
    void handle_coreblas_##func##_start (struct fxt_ev_64 *ev)          \
    {									\
        FUNC_NAME;                                                      \
        INIT_THREAD_ID(_threadstr);                                     \
        if ( GET_NBPARAMS(ev) > 0 ) {                                   \
            CHANGE() setState (CURRENT, COREBLAS_STATE, _threadstr, #func); \
        } else {                                                        \
            CHANGE() setState (CURRENT, COREBLAS_STATE, _threadstr, #func); \
        }                                                               \
        free(_threadstr);                                               \
    }
#endif

void handle_coreblas_task (struct fxt_ev_64 *ev)
{
    FUNC_NAME;
    INIT_PROCESS_ID(process_id);
    assert( GET_NBPARAMS(ev) == 1 );
    int value = (int)GET_PARAM(ev, 1);
    CHANGE() addVar (CURRENT, COREBLAS_TASK_ALIAS, process_id, (varPrec)value);
    free(process_id);
}

void handle_coreblas_taskw (struct fxt_ev_64 *ev)
{
    FUNC_NAME;
    assert( GET_NBPARAMS(ev) == 2 );
    INIT_PROCESS_ID(process_id);
    INIT_SPECIFIC_THREAD_ID(thread_id, CUR_ID, (unsigned int)GET_PARAM(ev, 1));
    int value = (int)GET_PARAM(ev, 2);
    CHANGE() addVar (CURRENT, COREBLAS_TASKR_ALIAS,  process_id, (varPrec)value);
    CHANGE() addVar (CURRENT, COREBLAS_TASKWR_ALIAS, thread_id,  (varPrec)value);
    free(thread_id);
    free(process_id);
}

/* Level 3 Blas */
HANDLE(gemm )
HANDLE(herk )
HANDLE(syrk )
HANDLE(hemm )
HANDLE(symm )
HANDLE(trmm )
HANDLE(trsm )
HANDLE(her2k)
HANDLE(syr2k)

/* Level 2 Blas */
HANDLE(gemv )
HANDLE(gbmv )
HANDLE(hemv )
HANDLE(hbmv )
HANDLE(hpmv )
HANDLE(symv )
HANDLE(sbmv )
HANDLE(spmv )
HANDLE(trmv )
HANDLE(tbmv )
HANDLE(tpmv )
HANDLE(trsv )
HANDLE(tbsv )
HANDLE(tpsv )
HANDLE(ger  )
HANDLE(geru )
HANDLE(gerc )
HANDLE(her  )
HANDLE(hpr  )
HANDLE(her2 )
HANDLE(hpr2 )
HANDLE(syr  )
HANDLE(spr  )
HANDLE(syr2 )
HANDLE(spr2 )

/* Level 1 BLAS */
HANDLE(rotg )
HANDLE(rotmg)
HANDLE(rot  )
HANDLE(rotm )
HANDLE(swap )
HANDLE(scal )
HANDLE(copy )
HANDLE(axpy )
HANDLE(dot  )
HANDLE(dotu )
HANDLE(dotc )
HANDLE(xdot )
HANDLE(nrm2 )
HANDLE(asum )
HANDLE(amax )

/* lapack */
HANDLE(lacpy)
HANDLE(lange)
HANDLE(lanhe)
HANDLE(lansy)
HANDLE(larfb)
HANDLE(larft)
HANDLE(laswp)
HANDLE(lauum)
HANDLE(potrf)
HANDLE(trtri)
HANDLE(laset)

/* plasma coreblas */
HANDLE(gelqt)
HANDLE(geqrt)
HANDLE(gessm)
HANDLE(getrf)
HANDLE(getro)
HANDLE(ssssm)
HANDLE(titro)
HANDLE(trbmm)
HANDLE(trgmm)
HANDLE(tslqt)
HANDLE(tsmlq)
HANDLE(tsmqr)
HANDLE(tsqrt)
HANDLE(tsrfb)
HANDLE(tstrf)
HANDLE(ttlqt)
HANDLE(ttmlq)
HANDLE(ttmqr)
HANDLE(ttqrt)
HANDLE(ttrfb)
HANDLE(unmlq)
HANDLE(unmqr)
HANDLE(getrip)
HANDLE(plghe)
HANDLE(plgsy)
HANDLE(shift)
HANDLE(shiftw)
HANDLE(swpab)
HANDLE(plrnt)

HANDLE(brdalg)
HANDLE(trdalg)
HANDLE(hegst)
HANDLE(sygst)
HANDLE(herfb)
HANDLE(syrfb)

void
handle_coreblas_stop ()
{
    FUNC_NAME;
    INIT_THREAD_ID(_threadstr);
    CHANGE() setState (CURRENT, COREBLAS_STATE, _threadstr, "wait");
    free(_threadstr);
}

int
eztrace_convert_coreblas_init()
{
    addVarType( COREBLAS_TASK_ALIAS,   COREBLAS_TASK_NAME,   "CT_Process" );
    addVarType( COREBLAS_TASKR_ALIAS,  COREBLAS_TASKR_NAME,  "CT_Process" );
    addVarType( COREBLAS_TASKWR_ALIAS, COREBLAS_TASKWR_NAME, "CT_Thread"  );

    /* Level 3 Blas */
    addEntityValue("gemm" , COREBLAS_STATE, "gemm" , GTG_YELLOW );
    addEntityValue("herk" , COREBLAS_STATE, "herk" , GTG_WHITE  );
    addEntityValue("syrk" , COREBLAS_STATE, "syrk" , GTG_WHITE  );
    addEntityValue("hemm" , COREBLAS_STATE, "hemm" , GTG_DARKPINK);
    addEntityValue("symm" , COREBLAS_STATE, "symm" , GTG_DARKPINK);
    addEntityValue("trmm" , COREBLAS_STATE, "trmm" , GTG_PURPLE );
    addEntityValue("trsm" , COREBLAS_STATE, "trsm" , GTG_RED    );
    addEntityValue("her2k", COREBLAS_STATE, "her2k", GTG_PINK   );
    addEntityValue("syr2k", COREBLAS_STATE, "syr2k", GTG_PINK   );

    /* Level 2 Blas */
    addEntityValue("gemv" , COREBLAS_STATE, "gemv" , GTG_TEAL   );
    addEntityValue("gbmv" , COREBLAS_STATE, "gbmv" , GTG_TEAL   );
    addEntityValue("hemv" , COREBLAS_STATE, "hemv" , GTG_TEAL   );
    addEntityValue("hbmv" , COREBLAS_STATE, "hbmv" , GTG_TEAL   );
    addEntityValue("hpmv" , COREBLAS_STATE, "hpmv" , GTG_TEAL   );
    addEntityValue("symv" , COREBLAS_STATE, "symv" , GTG_TEAL   );
    addEntityValue("sbmv" , COREBLAS_STATE, "sbmv" , GTG_TEAL   );
    addEntityValue("spmv" , COREBLAS_STATE, "spmv" , GTG_TEAL   );
    addEntityValue("trmv" , COREBLAS_STATE, "trmv" , GTG_TEAL   );
    addEntityValue("tbmv" , COREBLAS_STATE, "tbmv" , GTG_TEAL   );
    addEntityValue("tpmv" , COREBLAS_STATE, "tpmv" , GTG_TEAL   );
    addEntityValue("trsv" , COREBLAS_STATE, "trsv" , GTG_ORANGE );
    addEntityValue("tbsv" , COREBLAS_STATE, "tbsv" , GTG_ORANGE );
    addEntityValue("tpsv" , COREBLAS_STATE, "tpsv" , GTG_ORANGE );
    addEntityValue("ger"  , COREBLAS_STATE, "ger"  , GTG_SEABLUE   );
    addEntityValue("geru" , COREBLAS_STATE, "geru" , GTG_SEABLUE   );
    addEntityValue("gerc" , COREBLAS_STATE, "gerc" , GTG_SEABLUE   );
    addEntityValue("her"  , COREBLAS_STATE, "her"  , GTG_SEABLUE   );
    addEntityValue("hpr"  , COREBLAS_STATE, "hpr"  , GTG_SEABLUE   );
    addEntityValue("her2" , COREBLAS_STATE, "her2" , GTG_SEABLUE   );
    addEntityValue("hpr2" , COREBLAS_STATE, "hpr2" , GTG_SEABLUE   );
    addEntityValue("syr"  , COREBLAS_STATE, "syr"  , GTG_SEABLUE   );
    addEntityValue("spr"  , COREBLAS_STATE, "spr"  , GTG_SEABLUE   );
    addEntityValue("syr2" , COREBLAS_STATE, "syr2" , GTG_SEABLUE   );
    addEntityValue("spr2" , COREBLAS_STATE, "spr2" , GTG_SEABLUE   );

    /* Level 1 BLAS */
    addEntityValue("rotg" , COREBLAS_STATE, "rotg" , GTG_PURPLE   );
    addEntityValue("rotmg", COREBLAS_STATE, "rotmg", GTG_PURPLE   );
    addEntityValue("rot"  , COREBLAS_STATE, "rot"  , GTG_PURPLE   );
    addEntityValue("rotm" , COREBLAS_STATE, "rotm" , GTG_PURPLE   );
    addEntityValue("swap" , COREBLAS_STATE, "swap" , GTG_ORANGE   );
    addEntityValue("scal" , COREBLAS_STATE, "scal" , GTG_ORANGE   );
    addEntityValue("copy" , COREBLAS_STATE, "copy" , GTG_ORANGE   );
    addEntityValue("axpy" , COREBLAS_STATE, "axpy" , GTG_ORANGE   );
    addEntityValue("dot"  , COREBLAS_STATE, "dot"  , GTG_LIGHTPINK);
    addEntityValue("dotu" , COREBLAS_STATE, "dotu" , GTG_LIGHTPINK);
    addEntityValue("dotc" , COREBLAS_STATE, "dotc" , GTG_LIGHTPINK);
    addEntityValue("xdot" , COREBLAS_STATE, "xdot" , GTG_LIGHTPINK);
    addEntityValue("nrm2" , COREBLAS_STATE, "nrm2" , GTG_LIGHTPINK);
    addEntityValue("asum" , COREBLAS_STATE, "asum" , GTG_LIGHTPINK);
    addEntityValue("amax" , COREBLAS_STATE, "amax" , GTG_LIGHTPINK);

    /* Lapack */
    addEntityValue("lacpy", COREBLAS_STATE, "lacpy", GTG_LIGHTPINK );
    addEntityValue("lange", COREBLAS_STATE, "lange", GTG_LIGHTPINK );
    addEntityValue("lanhe", COREBLAS_STATE, "lanhe", GTG_LIGHTPINK );
    addEntityValue("lansy", COREBLAS_STATE, "lansy", GTG_LIGHTPINK );
    addEntityValue("larfb", COREBLAS_STATE, "larfb", GTG_YELLOW    );
    addEntityValue("larft", COREBLAS_STATE, "larft", GTG_RED       );
    addEntityValue("laswp", COREBLAS_STATE, "laswp", GTG_ORANGE    );
    addEntityValue("lauum", COREBLAS_STATE, "lauum", GTG_LIGHTPINK );
    addEntityValue("potrf", COREBLAS_STATE, "potrf", GTG_GREEN     );
    addEntityValue("trtri", COREBLAS_STATE, "trtri", GTG_LIGHTPINK );
    addEntityValue("laset", COREBLAS_STATE, "laset", GTG_LIGHTPINK );

    /* PLASMA coreblas */
    addEntityValue("gelqt", COREBLAS_STATE, "gelqt", GTG_GREEN     );
    addEntityValue("geqrt", COREBLAS_STATE, "geqrt", GTG_GREEN     );
    addEntityValue("gessm", COREBLAS_STATE, "gessm", GTG_BLUE      );
    addEntityValue("getrf", COREBLAS_STATE, "getrf", GTG_GREEN     );
    addEntityValue("getro", COREBLAS_STATE, "getro", GTG_ORANGE    );
    addEntityValue("ssssm", COREBLAS_STATE, "ssssm", GTG_YELLOW    );
    addEntityValue("titro", COREBLAS_STATE, "titro", GTG_LIGHTPINK );
    addEntityValue("trbmm", COREBLAS_STATE, "trbmm", GTG_BLUE      );
    addEntityValue("trgmm", COREBLAS_STATE, "trgmm", GTG_BLUE      );
    addEntityValue("tslqt", COREBLAS_STATE, "tslqt", GTG_RED       );
    addEntityValue("tsmlq", COREBLAS_STATE, "tsmlq", GTG_YELLOW    );
    addEntityValue("tsmqr", COREBLAS_STATE, "tsmqr", GTG_YELLOW    );
    addEntityValue("tsqrt", COREBLAS_STATE, "tsqrt", GTG_RED       );
    addEntityValue("tsrfb", COREBLAS_STATE, "tsrfb", GTG_BLUE      );
    addEntityValue("tstrf", COREBLAS_STATE, "tstrf", GTG_BLUE      );
    addEntityValue("ttlqt", COREBLAS_STATE, "ttlqt", GTG_REDBLOOD  );
    addEntityValue("ttmlq", COREBLAS_STATE, "ttmlq", GTG_ORANGE    );
    addEntityValue("ttmqr", COREBLAS_STATE, "ttmqr", GTG_ORANGE    );
    addEntityValue("ttqrt", COREBLAS_STATE, "ttqrt", GTG_REDBLOOD  );
    addEntityValue("ttrfb", COREBLAS_STATE, "ttrfb", GTG_SEABLUE   );
    addEntityValue("unmlq", COREBLAS_STATE, "unmlq", GTG_YELLOW    );
    addEntityValue("unmqr", COREBLAS_STATE, "unmqr", GTG_YELLOW    );
    addEntityValue("getrip",COREBLAS_STATE, "getrip",GTG_LIGHTPINK );
    addEntityValue("plghe", COREBLAS_STATE, "plghe", GTG_LIGHTPINK );
    addEntityValue("plgsy", COREBLAS_STATE, "plgsy", GTG_LIGHTPINK );
    addEntityValue("shift", COREBLAS_STATE, "shift", GTG_LIGHTPINK );
    addEntityValue("shiftw",COREBLAS_STATE, "shiftw",GTG_LIGHTPINK );
    addEntityValue("swpab", COREBLAS_STATE, "swpab", GTG_LIGHTPINK );
    addEntityValue("plrnt", COREBLAS_STATE, "plrnt", GTG_LIGHTPINK );

    addEntityValue("brdalg", COREBLAS_STATE, "brdalg", GTG_LIGHTPINK );
    addEntityValue("trdalg", COREBLAS_STATE, "trdalg", GTG_LIGHTPINK );
    addEntityValue("hegst",  COREBLAS_STATE, "hegst",  GTG_LIGHTPINK );
    addEntityValue("sygst",  COREBLAS_STATE, "sygst",  GTG_LIGHTPINK );
    addEntityValue("herfb",  COREBLAS_STATE, "herfb",  GTG_LIGHTPINK );
    addEntityValue("syrfb",  COREBLAS_STATE, "syrfb",  GTG_LIGHTPINK );

    /* plasma coreblas */
    addEntityValue ("wait", COREBLAS_STATE, "wait", GTG_BLACK );

#ifdef TRACE_BY_SEQUENCE
    sequenceInit();
#endif
    return 0;
}

int 
handle_coreblas_events(struct fxt_ev_64 *ev)
{

    switch (ev->code)
        {
            /* Level 3 Blas */
        case FUT_COREBLAS_GEMM  : handle_coreblas_gemm_start(ev);  break;
        case FUT_COREBLAS_HERK  : handle_coreblas_herk_start(ev);  break;
        case FUT_COREBLAS_SYRK  : handle_coreblas_syrk_start(ev);  break;
        case FUT_COREBLAS_HEMM  : handle_coreblas_hemm_start(ev);  break;
        case FUT_COREBLAS_SYMM  : handle_coreblas_symm_start(ev);  break;
        case FUT_COREBLAS_TRMM  : handle_coreblas_trmm_start(ev);  break;
        case FUT_COREBLAS_TRSM  : handle_coreblas_trsm_start(ev);  break;
        case FUT_COREBLAS_HER2K : handle_coreblas_her2k_start(ev); break;
        case FUT_COREBLAS_SYR2K : handle_coreblas_syr2k_start(ev); break;

            /* Level 2 Blas */
        case FUT_COREBLAS_GEMV  : handle_coreblas_gemv_start(ev);  break;
        case FUT_COREBLAS_GBMV  : handle_coreblas_gbmv_start(ev);  break;
        case FUT_COREBLAS_HEMV  : handle_coreblas_hemv_start(ev);  break;
        case FUT_COREBLAS_HBMV  : handle_coreblas_hbmv_start(ev);  break;
        case FUT_COREBLAS_HPMV  : handle_coreblas_hpmv_start(ev);  break;
        case FUT_COREBLAS_SYMV  : handle_coreblas_symv_start(ev);  break;
        case FUT_COREBLAS_SBMV  : handle_coreblas_sbmv_start(ev);  break;
        case FUT_COREBLAS_SPMV  : handle_coreblas_spmv_start(ev);  break;
        case FUT_COREBLAS_TRMV  : handle_coreblas_trmv_start(ev);  break;
        case FUT_COREBLAS_TBMV  : handle_coreblas_tbmv_start(ev);  break;
        case FUT_COREBLAS_TPMV  : handle_coreblas_tpmv_start(ev);  break;
        case FUT_COREBLAS_TRSV  : handle_coreblas_trsv_start(ev);  break;
        case FUT_COREBLAS_TBSV  : handle_coreblas_tbsv_start(ev);  break;
        case FUT_COREBLAS_TPSV  : handle_coreblas_tpsv_start(ev);  break;
        case FUT_COREBLAS_GER   : handle_coreblas_ger_start(ev);   break;
        case FUT_COREBLAS_GERU  : handle_coreblas_geru_start(ev);  break;
        case FUT_COREBLAS_GERC  : handle_coreblas_gerc_start(ev);  break;
        case FUT_COREBLAS_HER   : handle_coreblas_her_start(ev);   break;
        case FUT_COREBLAS_HPR   : handle_coreblas_hpr_start(ev);   break;
        case FUT_COREBLAS_HER2  : handle_coreblas_her2_start(ev);  break;
        case FUT_COREBLAS_HPR2  : handle_coreblas_hpr2_start(ev);  break;
        case FUT_COREBLAS_SYR   : handle_coreblas_syr_start(ev);   break;
        case FUT_COREBLAS_SPR   : handle_coreblas_spr_start(ev);   break;
        case FUT_COREBLAS_SYR2  : handle_coreblas_syr2_start(ev);  break;
        case FUT_COREBLAS_SPR2  : handle_coreblas_spr2_start(ev);  break;

            /* Level 1 BLAS */
        case FUT_COREBLAS_ROTG  : handle_coreblas_rotg_start(ev);  break;
        case FUT_COREBLAS_ROTMG : handle_coreblas_rotmg_start(ev); break;
        case FUT_COREBLAS_ROT   : handle_coreblas_rot_start(ev);   break;
        case FUT_COREBLAS_ROTM  : handle_coreblas_rotm_start(ev);  break;
        case FUT_COREBLAS_SWAP  : handle_coreblas_swap_start(ev);  break;
        case FUT_COREBLAS_SCAL  : handle_coreblas_scal_start(ev);  break;
        case FUT_COREBLAS_COPY  : handle_coreblas_copy_start(ev);  break;
        case FUT_COREBLAS_AXPY  : handle_coreblas_axpy_start(ev);  break;
        case FUT_COREBLAS_DOT   : handle_coreblas_dot_start(ev);   break;
        case FUT_COREBLAS_DOTU  : handle_coreblas_dotu_start(ev);  break;
        case FUT_COREBLAS_DOTC  : handle_coreblas_dotc_start(ev);  break;
        case FUT_COREBLAS_xDOT  : handle_coreblas_xdot_start(ev);  break;
        case FUT_COREBLAS_NRM2  : handle_coreblas_nrm2_start(ev);  break;
        case FUT_COREBLAS_ASUM  : handle_coreblas_asum_start(ev);  break;
        case FUT_COREBLAS_AMAX  : handle_coreblas_amax_start(ev);  break;

            /* Lapack */
        case FUT_COREBLAS_LACPY : handle_coreblas_lacpy_start(ev);  break;
        case FUT_COREBLAS_LANGE : handle_coreblas_lange_start(ev);  break;
        case FUT_COREBLAS_LANHE : handle_coreblas_lanhe_start(ev);  break;
        case FUT_COREBLAS_LANSY : handle_coreblas_lansy_start(ev);  break;
        case FUT_COREBLAS_LARFB : handle_coreblas_larfb_start(ev);  break;
        case FUT_COREBLAS_LARFT : handle_coreblas_larft_start(ev);  break;
        case FUT_COREBLAS_LASWP : handle_coreblas_laswp_start(ev);  break;
        case FUT_COREBLAS_LAUUM : handle_coreblas_lauum_start(ev);  break;
        case FUT_COREBLAS_POTRF : handle_coreblas_potrf_start(ev);  break;
        case FUT_COREBLAS_TRTRI : handle_coreblas_trtri_start(ev);  break;
        case FUT_COREBLAS_LASET : handle_coreblas_laset_start(ev);  break;

            /* PLASMA coreblas */
        case FUT_COREBLAS_GELQT : handle_coreblas_gelqt_start(ev);  break;
        case FUT_COREBLAS_GEQRT : handle_coreblas_geqrt_start(ev);  break;
        case FUT_COREBLAS_GESSM : handle_coreblas_gessm_start(ev);  break;
        case FUT_COREBLAS_GETRF : handle_coreblas_getrf_start(ev);  break;
        case FUT_COREBLAS_GETRO : handle_coreblas_getro_start(ev);  break;
        case FUT_COREBLAS_SSSSM : handle_coreblas_ssssm_start(ev);  break;
        case FUT_COREBLAS_TITRO : handle_coreblas_titro_start(ev);  break;
        case FUT_COREBLAS_TRBMM : handle_coreblas_trbmm_start(ev);  break;
        case FUT_COREBLAS_TRGMM : handle_coreblas_trgmm_start(ev);  break;
        case FUT_COREBLAS_TSLQT : handle_coreblas_tslqt_start(ev);  break;
        case FUT_COREBLAS_TSMLQ : handle_coreblas_tsmlq_start(ev);  break;
        case FUT_COREBLAS_TSMQR : handle_coreblas_tsmqr_start(ev);  break;
        case FUT_COREBLAS_TSQRT : handle_coreblas_tsqrt_start(ev);  break;
        case FUT_COREBLAS_TSRFB : handle_coreblas_tsrfb_start(ev);  break;
        case FUT_COREBLAS_TSTRF : handle_coreblas_tstrf_start(ev);  break;
        case FUT_COREBLAS_TTLQT : handle_coreblas_ttlqt_start(ev);  break;
        case FUT_COREBLAS_TTMLQ : handle_coreblas_ttmlq_start(ev);  break;
        case FUT_COREBLAS_TTMQR : handle_coreblas_ttmqr_start(ev);  break;
        case FUT_COREBLAS_TTQRT : handle_coreblas_ttqrt_start(ev);  break;
        case FUT_COREBLAS_TTRFB : handle_coreblas_ttrfb_start(ev);  break;
        case FUT_COREBLAS_UNMLQ : handle_coreblas_unmlq_start(ev);  break;
        case FUT_COREBLAS_UNMQR : handle_coreblas_unmqr_start(ev);  break;
        case FUT_COREBLAS_GETRIP: handle_coreblas_getrip_start(ev); break;
        case FUT_COREBLAS_PLGHE : handle_coreblas_plghe_start(ev);  break;
        case FUT_COREBLAS_PLGSY : handle_coreblas_plgsy_start(ev);  break;
        case FUT_COREBLAS_SHIFT : handle_coreblas_shift_start(ev);  break;
        case FUT_COREBLAS_SHIFTW: handle_coreblas_shiftw_start(ev); break;
        case FUT_COREBLAS_SWPAB : handle_coreblas_swpab_start(ev);  break;
        case FUT_COREBLAS_PLRNT : handle_coreblas_plrnt_start(ev);  break;

        case FUT_COREBLAS_BRDALG : handle_coreblas_brdalg_start(ev);  break;
        case FUT_COREBLAS_TRDALG : handle_coreblas_trdalg_start(ev);  break;
        case FUT_COREBLAS_HEGST  : handle_coreblas_hegst_start(ev);  break;
        case FUT_COREBLAS_SYGST  : handle_coreblas_sygst_start(ev);  break;
        case FUT_COREBLAS_HERFB  : handle_coreblas_herfb_start(ev);  break;
        case FUT_COREBLAS_SYRFB  : handle_coreblas_syrfb_start(ev);  break;

        case FUT_COREBLAS_STOP  : handle_coreblas_stop(ev);  break;
        case FUT_COREBLAS_TASK  : handle_coreblas_task(ev);  break;
        case FUT_COREBLAS_TASKW : handle_coreblas_taskw(ev); break;

        default:
            return 0;
        }
  
    return 1;
  
}

void
eztrace_convert_coreblas_finalize()
{

}


int 
handle_coreblas_stats(struct fxt_ev_64 *ev)
{
    int i;
    double time;

    if ( statsarray == NULL ) {
        statsarray = (coreblas_stats_t *)malloc(COREBLAS_NBMAX_EVENTS * sizeof(coreblas_stats_t));
        memset(statsarray, 0, COREBLAS_NBMAX_EVENTS * sizeof(coreblas_stats_t));

        thrdstate = (coreblas_thrdstate_t*)malloc(COREBLAS_THREADS_MAX * sizeof(coreblas_thrdstate_t));
        memset( thrdstate, 0, COREBLAS_THREADS_MAX * sizeof(coreblas_thrdstate_t));
    }   

    switch (ev->code)
        {
        case FUT_COREBLAS_STOP  : 
          {
              for (i=0; i<nbtrhd; i++) {
                  if ( thrdstate[i].tid == (unsigned int)CUR_THREAD_ID) {
                      if ( thrdstate[i].active == 0 ) {
                          fprintf(stderr, "WARNING: The end of a state appears before the beginning\n");
                          return 0;
                      }
                      
                      time = ( CURRENT - thrdstate[i].lasttime );

                      if( statsarray[ thrdstate[i].active ].nb == 0 ) {
                          statsarray[ thrdstate[i].active ].sum = 0.;
                          statsarray[ thrdstate[i].active ].max = 0.;
                          statsarray[ thrdstate[i].active ].min = 999999999999.;                          
                      }
                      statsarray[ thrdstate[i].active ].nb++;
                      statsarray[ thrdstate[i].active ].sum += time;
                      statsarray[ thrdstate[i].active ].max = max( statsarray[ thrdstate[i].active ].max, time );
                      statsarray[ thrdstate[i].active ].min = min( statsarray[ thrdstate[i].active ].min, time );

                      thrdstate[i].active = 0;
                      thrdstate[i].lasttime = 0;
                      return 1;
                  }
              }
              return 0;
          }
          break;

        case FUT_COREBLAS_TASK  : 
          break;
        case FUT_COREBLAS_TASKW : 
          break;

        default: /* All the different states */
            if ( ( (ev->code) & COREBLAS_PREFIX) ) {
                for (i=0; i<nbtrhd; i++) {
                    if ( thrdstate[i].tid == (unsigned int)CUR_THREAD_ID) {
                        if ( thrdstate[i].active != 0 ) {
                            fprintf(stderr, "WARNING: thread %d change to state %d before to stop previous state %d\n",
                                    (int)CUR_THREAD_ID, thrdstate[i].active,  (int)( (ev->code) & COREBLAS_MASK_EVENTS));
                        }
                        
                        thrdstate[i].active = (ev->code) & COREBLAS_MASK_EVENTS;
                        thrdstate[i].lasttime = CURRENT;
                        return 1;
                    }
                }

                /* Thread not found, we add it */
                if ( nbtrhd < COREBLAS_THREADS_MAX ) {
                    thrdstate[nbtrhd].tid = (unsigned int)CUR_THREAD_ID;
                    thrdstate[i].active = ev->code & COREBLAS_MASK_EVENTS;
                    thrdstate[nbtrhd].lasttime = CURRENT;
                    nbtrhd++;
                    return 1;
                }
            }
            return 0;
        }
  
    return 1;
}

/* 
 * Print the results of statistics.
 */
void print_coreblas_stats() {
    int i;

    coreblas_stats_strings_init();

    printf ( "\nCoreblas Module:\n");
    printf (   "-----------\n");

    for(i=0; i<COREBLAS_NBMAX_EVENTS; i++) {
        if ( statsarray[ i ].nb > 0 ) {
            printf ( "%s : %d calls\n"
                     "\tAverage time: %.3f ms\n"
                     "\tMaximun time: %.3f ms\n"
                     "\tMinimun time: %.3f ms\n",
                     coreblas_stats_strings[ i ], statsarray[ i ].nb,
                     statsarray[ i ].sum / (double)(statsarray[ i ].nb), statsarray[ i ].max, statsarray[ i ].min);
        }
    }
}

struct eztrace_convert_module coreblas_module;

void libinit(void) __attribute__ ((constructor));
void libinit(void)
{
  /* Specify the initialization function.
   * This function will be called once all the plugins are loaded
   * and the trace is started.
   * This function usually declared StateTypes, LinkTypes, etc.
   */
  coreblas_module.init = eztrace_convert_coreblas_init;

  /* Specify the function to call for handling an event
   */
  coreblas_module.handle = handle_coreblas_events;

  /* Specify the function to call for handling an event when eztrace_stats is called
   */
  coreblas_module.handle_stats = handle_coreblas_stats;

  /* Print the results of statistics
   */
  coreblas_module.print_stats = print_coreblas_stats;

  /* Specify the module prefix */
  coreblas_module.module_prefix = COREBLAS_EVENTS_ID;

  asprintf(&coreblas_module.name, "coreblas");
  asprintf(&coreblas_module.description, "Module for kernels used in PLASMA (BLAS, LAPACK and coreblas)");

  coreblas_module.token.data = &coreblas_module;

  /* Register the module to eztrace_convert */
  eztrace_convert_register_module(&coreblas_module);

  printf("module Coreblas loaded\n");
}

void libfinalize(void) __attribute__ ((destructor));
void libfinalize(void)
{
#ifdef TRACE_BY_SEQUENCE
    sequenceDestroy();
#endif
}
