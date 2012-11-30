/**
 *
 * @file global.h
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Jakub Kurzak
 * @author Piotr Luszczek
 * @date 2010-11-15
 *
 **/

/***************************************************************************//**
 *  PLASMA internals of interest to PLASMA core developers, but not necessarily
 *  of interest to PLASMA community contributors.
 **/
#ifndef _PLASMA_GLOBAL_H_
#define _PLASMA_GLOBAL_H_

#include <plasma.h>

#include <string.h>

#if defined( _WIN32 ) || defined( _WIN64 )
#include "plasmawinthread.h"
#else
#include <pthread.h>
#endif

/***************************************************************************//**
 *  Configuration
 **/
// maximum contexts
#define CONTEXTS_MAX         256
// maximum cores per context
#define CONTEXT_THREADS_MAX  256
// size of parallel functions arguments buffer
#define ARGS_BUFF_SIZE       512
// cache line size
#define CACHE_LINE_SIZE      128
// standard page size
#define STANDARD_PAGE_SIZE  4096

/***************************************************************************//**
 *  Action commands
 **/
#define PLASMA_ACT_STAND_BY     0
#define PLASMA_ACT_PARALLEL     1
#define PLASMA_ACT_DYNAMIC      2
#define PLASMA_ACT_FINALIZE     3

/***************************************************************************//**
 *  Numerical operations
 **/
#define PLASMA_FUNC_SGELS    1
#define PLASMA_FUNC_SPOSV    2
#define PLASMA_FUNC_SGESV    3
#define PLASMA_FUNC_DGELS    4
#define PLASMA_FUNC_DPOSV    5
#define PLASMA_FUNC_DGESV    6
#define PLASMA_FUNC_CGELS    7
#define PLASMA_FUNC_CPOSV    8
#define PLASMA_FUNC_CGESV    9
#define PLASMA_FUNC_ZGELS   10
#define PLASMA_FUNC_ZPOSV   11
#define PLASMA_FUNC_ZGESV   12
#define PLASMA_FUNC_ZCGESV  13
#define PLASMA_FUNC_DSGESV  14
#define PLASMA_FUNC_ZCPOSV  15
#define PLASMA_FUNC_DSPOSV  16
#define PLASMA_FUNC_DSGELS  17
#define PLASMA_FUNC_ZCGELS  18
#define PLASMA_FUNC_SGEMM   19
#define PLASMA_FUNC_DGEMM   20
#define PLASMA_FUNC_CGEMM   21
#define PLASMA_FUNC_ZGEMM   22
#define PLASMA_FUNC_SSYMM   23
#define PLASMA_FUNC_DSYMM   24
#define PLASMA_FUNC_CSYMM   25
#define PLASMA_FUNC_ZSYMM   26
#define PLASMA_FUNC_CHERK   27
#define PLASMA_FUNC_ZHERK   28
#define PLASMA_FUNC_SSYRK   29
#define PLASMA_FUNC_DSYRK   30
#define PLASMA_FUNC_CSYRK   31
#define PLASMA_FUNC_ZSYRK   32
#define PLASMA_FUNC_CHEMM   33
#define PLASMA_FUNC_ZHEMM   34
#define PLASMA_FUNC_ZHEEV   35
#define PLASMA_FUNC_CHEEV   36
#define PLASMA_FUNC_DSYEV   37
#define PLASMA_FUNC_SSYEV   38
#define PLASMA_FUNC_ZHEGST  39
#define PLASMA_FUNC_CHEGST  40
#define PLASMA_FUNC_DSYGST  41
#define PLASMA_FUNC_SSYGST  42
#define PLASMA_FUNC_ZHEGV   43
#define PLASMA_FUNC_CHEGV   44
#define PLASMA_FUNC_DSYGV   45
#define PLASMA_FUNC_SSYGV   46
#define PLASMA_FUNC_ZHETRD  47
#define PLASMA_FUNC_CHETRD  48
#define PLASMA_FUNC_DSYTRD  49
#define PLASMA_FUNC_SSYTRD  50
#define PLASMA_FUNC_ZGESVD  51
#define PLASMA_FUNC_CGESVD  52
#define PLASMA_FUNC_DGESVD  53
#define PLASMA_FUNC_SGESVD  54
#define PLASMA_FUNC_ZGEEV   55
#define PLASMA_FUNC_CGEEV   56
#define PLASMA_FUNC_DGEEV   57
#define PLASMA_FUNC_SGEEV   58
#define PLASMA_FUNC_ZGEHRD  59
#define PLASMA_FUNC_CGEHRD  60
#define PLASMA_FUNC_DGEHRD  61
#define PLASMA_FUNC_SGEHRD  62
#define PLASMA_FUNC_ZGEBRD  63
#define PLASMA_FUNC_CGEBRD  64
#define PLASMA_FUNC_DGEBRD  65
#define PLASMA_FUNC_SGEBRD  66

/***************************************************************************//**
 *  Parallel function call - packing of arguments
 **/
#define plasma_pack_args_1( \
    type1, arg1) \
{ \
    type1 var1 = (arg1); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_1", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
}

#define plasma_pack_args_2( \
    type1, arg1, \
    type2, arg2) \
{ \
    type1 var1 = (arg1); \
    type2 var2 = (arg2); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) + \
        sizeof(type2) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_2", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
    memcpy(plasma_ptr, &var2, sizeof(type2)); plasma_ptr += sizeof(type2); \
}

#define plasma_pack_args_3( \
    type1, arg1, \
    type2, arg2, \
    type3, arg3) \
{ \
    type1 var1 = (arg1); \
    type2 var2 = (arg2); \
    type3 var3 = (arg3); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) + \
        sizeof(type2) + \
        sizeof(type3) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_3", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
    memcpy(plasma_ptr, &var2, sizeof(type2)); plasma_ptr += sizeof(type2); \
    memcpy(plasma_ptr, &var3, sizeof(type3)); plasma_ptr += sizeof(type3); \
}

#define plasma_pack_args_4( \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4) \
{ \
    type1 var1 = (arg1); \
    type2 var2 = (arg2); \
    type3 var3 = (arg3); \
    type4 var4 = (arg4); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) + \
        sizeof(type2) + \
        sizeof(type3) + \
        sizeof(type4) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_4", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
    memcpy(plasma_ptr, &var2, sizeof(type2)); plasma_ptr += sizeof(type2); \
    memcpy(plasma_ptr, &var3, sizeof(type3)); plasma_ptr += sizeof(type3); \
    memcpy(plasma_ptr, &var4, sizeof(type4)); plasma_ptr += sizeof(type4); \
}

#define plasma_pack_args_5( \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5) \
{ \
    type1 var1 = (arg1); \
    type2 var2 = (arg2); \
    type3 var3 = (arg3); \
    type4 var4 = (arg4); \
    type5 var5 = (arg5); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) + \
        sizeof(type2) + \
        sizeof(type3) + \
        sizeof(type4) + \
        sizeof(type5) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_5", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
    memcpy(plasma_ptr, &var2, sizeof(type2)); plasma_ptr += sizeof(type2); \
    memcpy(plasma_ptr, &var3, sizeof(type3)); plasma_ptr += sizeof(type3); \
    memcpy(plasma_ptr, &var4, sizeof(type4)); plasma_ptr += sizeof(type4); \
    memcpy(plasma_ptr, &var5, sizeof(type5)); plasma_ptr += sizeof(type5); \
}

#define plasma_pack_args_6( \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6) \
{ \
    type1 var1 = (arg1); \
    type2 var2 = (arg2); \
    type3 var3 = (arg3); \
    type4 var4 = (arg4); \
    type5 var5 = (arg5); \
    type6 var6 = (arg6); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) + \
        sizeof(type2) + \
        sizeof(type3) + \
        sizeof(type4) + \
        sizeof(type5) + \
        sizeof(type6) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_6", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
    memcpy(plasma_ptr, &var2, sizeof(type2)); plasma_ptr += sizeof(type2); \
    memcpy(plasma_ptr, &var3, sizeof(type3)); plasma_ptr += sizeof(type3); \
    memcpy(plasma_ptr, &var4, sizeof(type4)); plasma_ptr += sizeof(type4); \
    memcpy(plasma_ptr, &var5, sizeof(type5)); plasma_ptr += sizeof(type5); \
    memcpy(plasma_ptr, &var6, sizeof(type6)); plasma_ptr += sizeof(type6); \
}

#define plasma_pack_args_7( \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7) \
{ \
    type1 var1 = (arg1); \
    type2 var2 = (arg2); \
    type3 var3 = (arg3); \
    type4 var4 = (arg4); \
    type5 var5 = (arg5); \
    type6 var6 = (arg6); \
    type7 var7 = (arg7); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) + \
        sizeof(type2) + \
        sizeof(type3) + \
        sizeof(type4) + \
        sizeof(type5) + \
        sizeof(type6) + \
        sizeof(type7) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_7", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
    memcpy(plasma_ptr, &var2, sizeof(type2)); plasma_ptr += sizeof(type2); \
    memcpy(plasma_ptr, &var3, sizeof(type3)); plasma_ptr += sizeof(type3); \
    memcpy(plasma_ptr, &var4, sizeof(type4)); plasma_ptr += sizeof(type4); \
    memcpy(plasma_ptr, &var5, sizeof(type5)); plasma_ptr += sizeof(type5); \
    memcpy(plasma_ptr, &var6, sizeof(type6)); plasma_ptr += sizeof(type6); \
    memcpy(plasma_ptr, &var7, sizeof(type7)); plasma_ptr += sizeof(type7); \
}

#define plasma_pack_args_8( \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8) \
{ \
    type1 var1 = (arg1); \
    type2 var2 = (arg2); \
    type3 var3 = (arg3); \
    type4 var4 = (arg4); \
    type5 var5 = (arg5); \
    type6 var6 = (arg6); \
    type7 var7 = (arg7); \
    type8 var8 = (arg8); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) + \
        sizeof(type2) + \
        sizeof(type3) + \
        sizeof(type4) + \
        sizeof(type5) + \
        sizeof(type6) + \
        sizeof(type7) + \
        sizeof(type8) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_8", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
    memcpy(plasma_ptr, &var2, sizeof(type2)); plasma_ptr += sizeof(type2); \
    memcpy(plasma_ptr, &var3, sizeof(type3)); plasma_ptr += sizeof(type3); \
    memcpy(plasma_ptr, &var4, sizeof(type4)); plasma_ptr += sizeof(type4); \
    memcpy(plasma_ptr, &var5, sizeof(type5)); plasma_ptr += sizeof(type5); \
    memcpy(plasma_ptr, &var6, sizeof(type6)); plasma_ptr += sizeof(type6); \
    memcpy(plasma_ptr, &var7, sizeof(type7)); plasma_ptr += sizeof(type7); \
    memcpy(plasma_ptr, &var8, sizeof(type8)); plasma_ptr += sizeof(type8); \
}

#define plasma_pack_args_9( \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9) \
{ \
    type1 var1 = (arg1); \
    type2 var2 = (arg2); \
    type3 var3 = (arg3); \
    type4 var4 = (arg4); \
    type5 var5 = (arg5); \
    type6 var6 = (arg6); \
    type7 var7 = (arg7); \
    type8 var8 = (arg8); \
    type9 var9 = (arg9); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) + \
        sizeof(type2) + \
        sizeof(type3) + \
        sizeof(type4) + \
        sizeof(type5) + \
        sizeof(type6) + \
        sizeof(type7) + \
        sizeof(type8) + \
        sizeof(type9) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_9", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
    memcpy(plasma_ptr, &var2, sizeof(type2)); plasma_ptr += sizeof(type2); \
    memcpy(plasma_ptr, &var3, sizeof(type3)); plasma_ptr += sizeof(type3); \
    memcpy(plasma_ptr, &var4, sizeof(type4)); plasma_ptr += sizeof(type4); \
    memcpy(plasma_ptr, &var5, sizeof(type5)); plasma_ptr += sizeof(type5); \
    memcpy(plasma_ptr, &var6, sizeof(type6)); plasma_ptr += sizeof(type6); \
    memcpy(plasma_ptr, &var7, sizeof(type7)); plasma_ptr += sizeof(type7); \
    memcpy(plasma_ptr, &var8, sizeof(type8)); plasma_ptr += sizeof(type8); \
    memcpy(plasma_ptr, &var9, sizeof(type9)); plasma_ptr += sizeof(type9); \
}

#define plasma_pack_args_10( \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10) \
{ \
    type1 var1 = (arg1); \
    type2 var2 = (arg2); \
    type3 var3 = (arg3); \
    type4 var4 = (arg4); \
    type5 var5 = (arg5); \
    type6 var6 = (arg6); \
    type7 var7 = (arg7); \
    type8 var8 = (arg8); \
    type9 var9 = (arg9); \
    type10 var10 = (arg10); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) + \
        sizeof(type2) + \
        sizeof(type3) + \
        sizeof(type4) + \
        sizeof(type5) + \
        sizeof(type6) + \
        sizeof(type7) + \
        sizeof(type8) + \
        sizeof(type9) + \
        sizeof(type10) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_9", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
    memcpy(plasma_ptr, &var2, sizeof(type2)); plasma_ptr += sizeof(type2); \
    memcpy(plasma_ptr, &var3, sizeof(type3)); plasma_ptr += sizeof(type3); \
    memcpy(plasma_ptr, &var4, sizeof(type4)); plasma_ptr += sizeof(type4); \
    memcpy(plasma_ptr, &var5, sizeof(type5)); plasma_ptr += sizeof(type5); \
    memcpy(plasma_ptr, &var6, sizeof(type6)); plasma_ptr += sizeof(type6); \
    memcpy(plasma_ptr, &var7, sizeof(type7)); plasma_ptr += sizeof(type7); \
    memcpy(plasma_ptr, &var8, sizeof(type8)); plasma_ptr += sizeof(type8); \
    memcpy(plasma_ptr, &var9, sizeof(type9)); plasma_ptr += sizeof(type9); \
    memcpy(plasma_ptr, &var10, sizeof(type10)); plasma_ptr += sizeof(type10); \
}

/***************************************************************************//**
 *  Sync after dynamically scheduled section
 **/
#define plasma_dynamic_sync() \
{ \
    if (plasma->dynamic_section) { \
        QUARK_Waitall(plasma->quark); \
        plasma_barrier(plasma); \
        plasma->dynamic_section = PLASMA_FALSE; \
    } \
}

/***************************************************************************//**
 *  Parallel SPMD function call - thread control
 **/
#define plasma_static_call(parallel_function) \
{ \
    if (plasma->dynamic_section) \
        plasma_dynamic_sync(); \
    pthread_mutex_lock(&plasma->action_mutex); \
    plasma->action = PLASMA_ACT_PARALLEL; \
    plasma->parallel_func_ptr = &parallel_function; \
    pthread_mutex_unlock(&plasma->action_mutex); \
    pthread_cond_broadcast(&plasma->action_condt); \
    plasma_barrier(plasma); \
    plasma->action = PLASMA_ACT_STAND_BY; \
    parallel_function(plasma); \
    plasma_barrier(plasma); \
}

/***************************************************************************//**
 *  Start dynamically scheduled section
 **/
#define plasma_dynamic_spawn() \
{ \
    if (!plasma->dynamic_section) { \
        plasma->dynamic_section = PLASMA_TRUE; \
        pthread_mutex_lock(&plasma->action_mutex); \
        plasma->action = PLASMA_ACT_DYNAMIC; \
        pthread_mutex_unlock(&plasma->action_mutex); \
        pthread_cond_broadcast(&plasma->action_condt); \
        plasma_barrier(plasma); \
        plasma->action = PLASMA_ACT_STAND_BY; \
    } \
}

/***************************************************************************//**
 *  Parallel call for functions with static versions only
 **/
#define plasma_static_call_1( \
           parallel_function, \
    type1, arg1) \
    plasma_pack_args_1( \
        type1, (arg1)) \
    plasma_static_call(parallel_function) \

#define plasma_static_call_2( \
           parallel_function, \
    type1, arg1, \
    type2, arg2) \
    plasma_pack_args_2( \
        type1, (arg1), \
        type2, (arg2)) \
    plasma_static_call(parallel_function) \

#define plasma_static_call_3( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3) \
    plasma_pack_args_3( \
        type1, (arg1), \
        type2, (arg2), \
        type3, (arg3)) \
    plasma_static_call(parallel_function) \

#define plasma_static_call_4( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4) \
    plasma_pack_args_4( \
        type1, (arg1), \
        type2, (arg2), \
        type3, (arg3), \
        type4, (arg4)) \
    plasma_static_call(parallel_function) \

#define plasma_static_call_5( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5) \
    plasma_pack_args_5( \
        type1, (arg1), \
        type2, (arg2), \
        type3, (arg3), \
        type4, (arg4), \
        type5, (arg5)) \
    plasma_static_call(parallel_function) \

#define plasma_static_call_6( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6) \
    plasma_pack_args_6( \
        type1, (arg1), \
        type2, (arg2), \
        type3, (arg3), \
        type4, (arg4), \
        type5, (arg5), \
        type6, (arg6)) \
    plasma_static_call(parallel_function) \

#define plasma_static_call_7( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7) \
    plasma_pack_args_7( \
        type1, (arg1), \
        type2, (arg2), \
        type3, (arg3), \
        type4, (arg4), \
        type5, (arg5), \
        type6, (arg6), \
        type7, (arg7)) \
    plasma_static_call(parallel_function) \

#define plasma_static_call_8( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8) \
    plasma_pack_args_8( \
        type1, (arg1), \
        type2, (arg2), \
        type3, (arg3), \
        type4, (arg4), \
        type5, (arg5), \
        type6, (arg6), \
        type7, (arg7), \
        type8, (arg8)) \
    plasma_static_call(parallel_function) \

#define plasma_static_call_9( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9) \
    plasma_pack_args_9( \
        type1, (arg1), \
        type2, (arg2), \
        type3, (arg3), \
        type4, (arg4), \
        type5, (arg5), \
        type6, (arg6), \
        type7, (arg7), \
        type8, (arg8), \
        type9, (arg9)) \
    plasma_static_call(parallel_function) \

#define plasma_static_call_10( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10) \
    plasma_pack_args_10( \
        type1, (arg1), \
        type2, (arg2), \
        type3, (arg3), \
        type4, (arg4), \
        type5, (arg5), \
        type6, (arg6), \
        type7, (arg7), \
        type8, (arg8), \
        type9, (arg9), \
        type10, (arg10)) \
    plasma_static_call(parallel_function) \

/***************************************************************************//**
 *  Parallel call for functions with both static and dynamic versions
 **/
#define plasma_parallel_call_1( \
           parallel_function, \
    type1, arg1) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_1( \
            type1, (arg1)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1); \
    }

#define plasma_parallel_call_2( \
           parallel_function, \
    type1, arg1, \
    type2, arg2) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_2( \
            type1, (arg1), \
            type2, (arg2)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1, \
            arg2); \
    }

#define plasma_parallel_call_3( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_3( \
            type1, (arg1), \
            type2, (arg2), \
            type3, (arg3)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1, \
            arg2, \
            arg3); \
    }

#define plasma_parallel_call_4( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_4( \
            type1, (arg1), \
            type2, (arg2), \
            type3, (arg3), \
            type4, (arg4)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1, \
            arg2, \
            arg3, \
            arg4); \
    }

#define plasma_parallel_call_5( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_5( \
            type1, (arg1), \
            type2, (arg2), \
            type3, (arg3), \
            type4, (arg4), \
            type5, (arg5)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1, \
            arg2, \
            arg3, \
            arg4, \
            arg5); \
    }

#define plasma_parallel_call_6( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_6( \
            type1, (arg1), \
            type2, (arg2), \
            type3, (arg3), \
            type4, (arg4), \
            type5, (arg5), \
            type6, (arg6)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1, \
            arg2, \
            arg3, \
            arg4, \
            arg5, \
            arg6); \
    }

#define plasma_parallel_call_7( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_7( \
            type1, (arg1), \
            type2, (arg2), \
            type3, (arg3), \
            type4, (arg4), \
            type5, (arg5), \
            type6, (arg6), \
            type7, (arg7)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1, \
            arg2, \
            arg3, \
            arg4, \
            arg5, \
            arg6, \
            arg7); \
    }

#define plasma_parallel_call_8( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_8( \
            type1, (arg1), \
            type2, (arg2), \
            type3, (arg3), \
            type4, (arg4), \
            type5, (arg5), \
            type6, (arg6), \
            type7, (arg7), \
            type8, (arg8)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1, \
            arg2, \
            arg3, \
            arg4, \
            arg5, \
            arg6, \
            arg7, \
            arg8); \
    }

#define plasma_parallel_call_9( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_9( \
            type1, (arg1), \
            type2, (arg2), \
            type3, (arg3), \
            type4, (arg4), \
            type5, (arg5), \
            type6, (arg6), \
            type7, (arg7), \
            type8, (arg8), \
            type9, (arg9)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1, \
            arg2, \
            arg3, \
            arg4, \
            arg5, \
            arg6, \
            arg7, \
            arg8, \
            arg9); \
    }

#define plasma_parallel_call_10( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_10( \
            type1, (arg1), \
            type2, (arg2), \
            type3, (arg3), \
            type4, (arg4), \
            type5, (arg5), \
            type6, (arg6), \
            type7, (arg7), \
            type8, (arg8), \
            type9, (arg9), \
            type10, (arg10)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1, \
            arg2, \
            arg3, \
            arg4, \
            arg5, \
            arg6, \
            arg7, \
            arg8, \
            arg9, \
            arg10); \
    }

/***************************************************************************//**
 *  Parallel call for functions with dynamic versions only
 **/
#define plasma_dynamic_call_1( \
           parallel_function, \
    type1, arg1) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1); \

#define plasma_dynamic_call_2( \
           parallel_function, \
    type1, arg1, \
    type2, arg2) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1, \
        arg2); \

#define plasma_dynamic_call_3( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1, \
        arg2, \
        arg3); \

#define plasma_dynamic_call_4( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1, \
        arg2, \
        arg3, \
        arg4); \

#define plasma_dynamic_call_5( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1, \
        arg2, \
        arg3, \
        arg4, \
        arg5); \

#define plasma_dynamic_call_6( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1, \
        arg2, \
        arg3, \
        arg4, \
        arg5, \
        arg6); \

#define plasma_dynamic_call_7( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1, \
        arg2, \
        arg3, \
        arg4, \
        arg5, \
        arg6, \
        arg7); \

#define plasma_dynamic_call_8( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1, \
        arg2, \
        arg3, \
        arg4, \
        arg5, \
        arg6, \
        arg7, \
        arg8); \

#define plasma_dynamic_call_9( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1, \
        arg2, \
        arg3, \
        arg4, \
        arg5, \
        arg6, \
        arg7, \
        arg8, \
        arg9);

#define plasma_dynamic_call_10( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1, \
        arg2, \
        arg3, \
        arg4, \
        arg5, \
        arg6, \
        arg7, \
        arg8, \
        arg9, \
        arg10);

/***************************************************************************//**
 *  Parallel function call - unpacking of arguments
 **/
#define plasma_unpack_args_1( \
    arg1) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
}

#define plasma_unpack_args_2( \
    arg1, \
    arg2) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
}

#define plasma_unpack_args_3( \
    arg1, \
    arg2, \
    arg3) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
    memcpy(&arg3, plasma_ptr, sizeof(arg3)); plasma_ptr += sizeof(arg3); \
}

#define plasma_unpack_args_4( \
    arg1, \
    arg2, \
    arg3, \
    arg4) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
    memcpy(&arg3, plasma_ptr, sizeof(arg3)); plasma_ptr += sizeof(arg3); \
    memcpy(&arg4, plasma_ptr, sizeof(arg4)); plasma_ptr += sizeof(arg4); \
}

#define plasma_unpack_args_5( \
    arg1, \
    arg2, \
    arg3, \
    arg4, \
    arg5) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
    memcpy(&arg3, plasma_ptr, sizeof(arg3)); plasma_ptr += sizeof(arg3); \
    memcpy(&arg4, plasma_ptr, sizeof(arg4)); plasma_ptr += sizeof(arg4); \
    memcpy(&arg5, plasma_ptr, sizeof(arg5)); plasma_ptr += sizeof(arg5); \
}

#define plasma_unpack_args_6( \
    arg1, \
    arg2, \
    arg3, \
    arg4, \
    arg5, \
    arg6) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
    memcpy(&arg3, plasma_ptr, sizeof(arg3)); plasma_ptr += sizeof(arg3); \
    memcpy(&arg4, plasma_ptr, sizeof(arg4)); plasma_ptr += sizeof(arg4); \
    memcpy(&arg5, plasma_ptr, sizeof(arg5)); plasma_ptr += sizeof(arg5); \
    memcpy(&arg6, plasma_ptr, sizeof(arg6)); plasma_ptr += sizeof(arg6); \
}

#define plasma_unpack_args_7( \
    arg1, \
    arg2, \
    arg3, \
    arg4, \
    arg5, \
    arg6, \
    arg7) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
    memcpy(&arg3, plasma_ptr, sizeof(arg3)); plasma_ptr += sizeof(arg3); \
    memcpy(&arg4, plasma_ptr, sizeof(arg4)); plasma_ptr += sizeof(arg4); \
    memcpy(&arg5, plasma_ptr, sizeof(arg5)); plasma_ptr += sizeof(arg5); \
    memcpy(&arg6, plasma_ptr, sizeof(arg6)); plasma_ptr += sizeof(arg6); \
    memcpy(&arg7, plasma_ptr, sizeof(arg7)); plasma_ptr += sizeof(arg7); \
}

#define plasma_unpack_args_8( \
    arg1, \
    arg2, \
    arg3, \
    arg4, \
    arg5, \
    arg6, \
    arg7, \
    arg8) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
    memcpy(&arg3, plasma_ptr, sizeof(arg3)); plasma_ptr += sizeof(arg3); \
    memcpy(&arg4, plasma_ptr, sizeof(arg4)); plasma_ptr += sizeof(arg4); \
    memcpy(&arg5, plasma_ptr, sizeof(arg5)); plasma_ptr += sizeof(arg5); \
    memcpy(&arg6, plasma_ptr, sizeof(arg6)); plasma_ptr += sizeof(arg6); \
    memcpy(&arg7, plasma_ptr, sizeof(arg7)); plasma_ptr += sizeof(arg7); \
    memcpy(&arg8, plasma_ptr, sizeof(arg8)); plasma_ptr += sizeof(arg8); \
}

#define plasma_unpack_args_9( \
    arg1, \
    arg2, \
    arg3, \
    arg4, \
    arg5, \
    arg6, \
    arg7, \
    arg8, \
    arg9) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
    memcpy(&arg3, plasma_ptr, sizeof(arg3)); plasma_ptr += sizeof(arg3); \
    memcpy(&arg4, plasma_ptr, sizeof(arg4)); plasma_ptr += sizeof(arg4); \
    memcpy(&arg5, plasma_ptr, sizeof(arg5)); plasma_ptr += sizeof(arg5); \
    memcpy(&arg6, plasma_ptr, sizeof(arg6)); plasma_ptr += sizeof(arg6); \
    memcpy(&arg7, plasma_ptr, sizeof(arg7)); plasma_ptr += sizeof(arg7); \
    memcpy(&arg8, plasma_ptr, sizeof(arg8)); plasma_ptr += sizeof(arg8); \
    memcpy(&arg9, plasma_ptr, sizeof(arg9)); plasma_ptr += sizeof(arg9); \
}

#define plasma_unpack_args_10( \
    arg1, \
    arg2, \
    arg3, \
    arg4, \
    arg5, \
    arg6, \
    arg7, \
    arg8, \
    arg9, \
    arg10) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
    memcpy(&arg3, plasma_ptr, sizeof(arg3)); plasma_ptr += sizeof(arg3); \
    memcpy(&arg4, plasma_ptr, sizeof(arg4)); plasma_ptr += sizeof(arg4); \
    memcpy(&arg5, plasma_ptr, sizeof(arg5)); plasma_ptr += sizeof(arg5); \
    memcpy(&arg6, plasma_ptr, sizeof(arg6)); plasma_ptr += sizeof(arg6); \
    memcpy(&arg7, plasma_ptr, sizeof(arg7)); plasma_ptr += sizeof(arg7); \
    memcpy(&arg8, plasma_ptr, sizeof(arg8)); plasma_ptr += sizeof(arg8); \
    memcpy(&arg9, plasma_ptr, sizeof(arg9)); plasma_ptr += sizeof(arg9); \
    memcpy(&arg10, plasma_ptr, sizeof(arg10)); plasma_ptr += sizeof(arg10); \
}
#endif
