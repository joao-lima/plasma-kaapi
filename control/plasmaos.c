/**
 *
 * @file plasmaos.c
 *
 *  This file handles the mapping from pthreads calls to windows threads
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Piotr Luszczek
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 **/

#if defined(linux) || defined(__linux) || defined(__linux__)
#define PLASMA_OS_LINUX 1
#define _GNU_SOURCE
#include <unistd.h>
#include <sched.h>
#elif defined( _WIN32 ) || defined( _WIN64 )
#define PLASMA_OS_WINDOWS 1
#include <Windows.h>
#elif (defined __APPLE__) || (defined macintosh) || (defined __MACOSX__)
#define PLASMA_OS_MACOS 1
#include <sys/param.h>
#include <sys/sysctl.h>
#include <mach/mach_init.h>
#include <mach/thread_policy.h>
kern_return_t thread_policy_set(thread_act_t thread, thread_policy_flavor_t flavor,
                                thread_policy_t policy_info, mach_msg_type_number_t count);
#elif (defined _AIX)
#define PLASMA_OS_AIX 1
#else
#error "Cannot find the runing system or system not supported. Please define try to PLASMA_OS_[LINUX|MACOS|AIX|WINDOWS]"
#endif

#if defined(PLASMA_HWLOC) && (defined PLASMA_AFFINITY_DISABLE)
#undef PLASMA_HWLOC
#endif

#include <errno.h>
#include <stdlib.h>
#include "common.h"

#ifdef __cplusplus
extern "C" {
#endif

static pthread_mutex_t  mutextopo = PTHREAD_MUTEX_INITIALIZER;
static volatile int sys_corenbr = 1;
static volatile int topo_initialized = 0;

  /*
   * Topology functions
   */
#ifdef PLASMA_HWLOC
#include "plasmaos-hwloc.c"
#else 

void plasma_topology_init(){
    pthread_mutex_lock(&mutextopo);
    if ( !topo_initialized ) {
#if (defined PLASMA_OS_LINUX) || (defined PLASMA_OS_AIX)

        sys_corenbr = sysconf(_SC_NPROCESSORS_ONLN);

#elif (defined PLASMA_OS_MACOS)

        int mib[4];
        int cpu;
        size_t len = sizeof(cpu);
        
        /* set the mib for hw.ncpu */
        mib[0] = CTL_HW;
        mib[1] = HW_AVAILCPU;

        /* get the number of CPUs from the system */
        sysctl(mib, 2, &cpu, &len, NULL, 0);
        if( cpu < 1 ) {
            mib[1] = HW_NCPU;
            sysctl( mib, 2, &cpu, &len, NULL, 0 );
        }
        if( cpu < 1 ) {
            cpu = 1;
        }
        sys_corenbr = cpu;
#elif (defined PLASMA_OS_WINDOWS)
        SYSTEM_INFO sysinfo;
        GetSystemInfo(&sysinfo);
        sys_corenbr = sysinfo.dwNumberOfProcessors;
#endif
    }
    pthread_mutex_unlock(&mutextopo);
}

void plasma_topology_finalize(){
    plasma_unsetaffinity();
}

/**
 This routine will set affinity for the calling thread that has rank 'rank'.
 Ranks start with 0.

 If there are multiple instances of PLASMA then affinity will be wrong: all ranks 0
 will be pinned to core 0.

 Also, affinity is not resotred when PLASMA_Finalize() is called.
 */
int plasma_setaffinity(int rank) {
#ifndef PLASMA_AFFINITY_DISABLE
#if (defined PLASMA_OS_LINUX)
    {
        cpu_set_t set;
        CPU_ZERO( &set );
        CPU_SET( rank, &set );
        
#if (defined HAVE_OLD_SCHED_SETAFFINITY)
        if( sched_setaffinity( 0, &set ) < 0 )
#else /* HAVE_OLD_SCHED_SETAFFINITY */
        if( sched_setaffinity( 0, sizeof(set), &set) < 0 )
#endif /* HAVE_OLD_SCHED_SETAFFINITY */
            {
                return PLASMA_ERR_UNEXPECTED;
            }

        return PLASMA_SUCCESS;
    }
#elif (defined PLASMA_OS_MACOS)
    {
        thread_affinity_policy_data_t ap;
        int                           ret;
        
        ap.affinity_tag = 1; /* non-null affinity tag */
        ret = thread_policy_set( mach_thread_self(),
                                 THREAD_AFFINITY_POLICY,
                                 (integer_t*) &ap,
                                 THREAD_AFFINITY_POLICY_COUNT
            );
        if(ret != 0) 
            return PLASMA_ERR_UNEXPECTED;

        return PLASMA_SUCCESS;
    }
#elif (defined PLASMA_OS_WINDOWS)
    {
        DWORD mask = 1 << rank;

        if( SetThreadAffinityMask(GetCurrentThread(), mask) == 0)
            return PLASMA_ERR_UNEXPECTED;
        
        return PLASMA_SUCCESS;
    }
#elif (defined PLASMA_OS_AIX)
    {
        tid_t self_ktid = thread_self ();
        bindprocessor(BINDTHREAD, self_ktid, rank);
        return PLASMA_SUCCESS;
    }
#else
    return PLASMA_ERR_NOT_SUPPORTED;
#endif
#endif /* PLASMA_AFFINITY_DISABLE */
}

/**
 This routine will set affinity for the calling thread that has rank 'rank'.
 Ranks start with 0.

 If there are multiple instances of PLASMA then affinity will be wrong: all ranks 0
 will be pinned to core 0.

 Also, affinity is not resotred when PLASMA_Finalize() is called.
 */
int plasma_unsetaffinity(int rank) {
#ifndef PLASMA_AFFINITY_DISABLE
#if (defined PLASMA_OS_LINUX)
    {
        int i;
        cpu_set_t set;
        CPU_ZERO( &set );

        for(i=0; i<sys_corenbr; i++)
            CPU_SET( i, &set );

#if (defined HAVE_OLD_SCHED_SETAFFINITY)
        if( sched_setaffinity( 0, &set ) < 0 )
#else /* HAVE_OLD_SCHED_SETAFFINITY */
        if( sched_setaffinity( 0, sizeof(set), &set) < 0 )
#endif /* HAVE_OLD_SCHED_SETAFFINITY */
            {
                plasma_warning("plasma_unsetaffinity", "Could not unbind thread");
                return PLASMA_ERR_UNEXPECTED;
            }

        return PLASMA_SUCCESS;
    }
#elif (defined PLASMA_OS_MACOS)
    {
        /* TODO: check how to unbind the main thread if necessary for OpenMP */
        /* thread_affinity_policy_data_t ap; */
        /* int                           ret; */
        
        /* ap.affinity_tag = 1; /\* non-null affinity tag *\/ */
        /* ret = thread_policy_set( mach_thread_self(), */
        /*                          THREAD_AFFINITY_POLICY, */
        /*                          (integer_t*) &ap, */
        /*                          THREAD_AFFINITY_POLICY_COUNT */
        /*     ); */
        /* if(ret != 0) { */
        /*     plasma_warning("plasma_unsetaffinity", "Could not unbind thread"); */
        /*     return PLASMA_ERR_UNEXPECTED; */
        /* } */

        return PLASMA_SUCCESS;
    }
#elif (defined PLASMA_OS_WINDOWS)
    {
        DWORD mask = 0;

        for(i=0; i<sys_corenbr; i++)
            mask |= 1 << i;

        if( SetThreadAffinityMask(GetCurrentThread(), mask) == 0) {
            plasma_warning("plasma_unsetaffinity", "Could not unbind thread");
            return PLASMA_ERR_UNEXPECTED;
        }
        return PLASMA_SUCCESS;
    }
#elif (defined PLASMA_OS_AIX)
    {
        /* TODO: check how to unbind the main thread if necessary for OpenMP */
        /* tid_t self_ktid = thread_self (); */
        /* bindprocessor(BINDTHREAD, self_ktid, rank); */
        return PLASMA_SUCCESS;
    }
#else
    return PLASMA_ERR_NOT_SUPPORTED;
#endif
#endif /* PLASMA_AFFINITY_DISABLE */
}
#endif /* PLASMA_HWLOC */

/** ****************************************************************************
   A thread can unlock the CPU if it has nothing to do to let 
   another thread of less priority running for example for I/O.
 */
int plasma_yield() {
#if (defined PLASMA_OS_LINUX) || (defined PLASMA_OS_MACOS) || (defined PLASMA_OS_AIX)
    return sched_yield();
#elif PLASMA_OS_WINDOWS
    return SleepEx(0,0);
#else
    return PLASMA_ERR_NOT_SUPPORTED;
#endif
}
    
#ifdef PLASMA_OS_WINDOWS
#define PLASMA_GETENV(var, str) {                    \
        int len = 512;                               \
        int ret;                                     \
        str = (char*)malloc(len * sizeof(char));     \
        ret = GetEnvironmentVariable(var, str, len); \
        if (ret == 0) {                              \
            free(str);                               \
            str = NULL;                              \
        }                                            \
    }

#define PLASMA_CLEANENV(str) if (str != NULL) free(str);

#else /* Other OS systems */

#define PLASMA_GETENV(var, str)  envstr = getenv(var);
#define PLASMA_CLEANENV(str)

#endif

/** ****************************************************************************
 * Check for an integer in an environment variable, returning the
 * integer value or a provided default value 
*/
int plasma_get_numthreads()
{
    char    *envstr  = NULL;
    char    *endptr;
    long int thrdnbr = -1;
    extern int errno;
    
    /* Env variable does not exist, we search the system number of core */
    PLASMA_GETENV("PLASMA_NUM_THREADS", envstr);
    if ( envstr == NULL ) {
        thrdnbr = sys_corenbr;
    } else {
        /* Convert to long, checking for errors */
        thrdnbr = strtol(envstr, &endptr, 10);
        if ((errno == ERANGE) || ((thrdnbr==0) && (endptr==envstr))) {
            PLASMA_CLEANENV(envstr);
            return -1;
        }
    } 
    PLASMA_CLEANENV(envstr);
    return (int)thrdnbr;
}

int plasma_get_numthreads_numa()
{
    char    *envstr  = NULL;
    char    *endptr;
    long int thrdnbr = -1;
    extern int errno;
    
    /* Env variable does not exist, we search the system number of core */
    PLASMA_GETENV("PLASMA_NUM_THREADS_NUMA", envstr);
    if ( envstr != NULL ) {
        /* Convert to long, checking for errors */
        thrdnbr = strtol(envstr, &endptr, 10);
        if ((errno == ERANGE) || ((thrdnbr==0) && (endptr==envstr))) {
            PLASMA_CLEANENV(envstr);
            return -1;
        }
    } else {
#ifdef PLASMA_HWLOC
      thrdnbr = plasma_getnuma_size();
#else
      thrdnbr = 1;
#endif
    } 

    PLASMA_CLEANENV(envstr);
    return (int)thrdnbr;
}

int plasma_get_affthreads(int *coresbind) {
    char *envstr = NULL;
    int i;

    /* Env variable does not exist, we search the system number of core */
    PLASMA_GETENV("PLASMA_AFF_THREADS", envstr);
    if ( envstr == NULL) {
        for (i = 0; i < CONTEXT_THREADS_MAX; i++)
            coresbind[i] = i % sys_corenbr;
    }
    else {
        char *endptr;
        int   wrap = 0;
        int   nbr  = 0;
        long int val;

        /* We use the content of the PLASMA_AFF_THREADS env. variable */
        for (i = 0; i < CONTEXT_THREADS_MAX; i++) {
            if (!wrap) {
                val = strtol(envstr, &endptr, 10);
                if (endptr != envstr) {
                    coresbind[i] = (int)val;
                    envstr = endptr;
                }
                else {
                    /* there must be at least one entry */
                    if (i < 1) {
                        plasma_error("plasma_get_affthreads", "PLASMA_AFF_THREADS should have at least one entry => everything will be bind on core 0");
                        coresbind[i] = 0;
                        i++;
                    }
                    
                    /* there is no more values in the string                                 */
                    /* the next threads are binded with a round robin policy over this array */
                    wrap = 1;
                    nbr = i;

                    coresbind[i] = coresbind[0];
                }
            }
            else {
                coresbind[i] = coresbind[i % nbr];
            }
        }
    }
    PLASMA_CLEANENV(envstr);
    return PLASMA_SUCCESS;
}

#ifdef __cplusplus
}
#endif

#ifndef PLASMA_HAS_COMPLEX_H

#ifdef __cplusplus
extern "C" {
#endif
float
cabsf(float _Complex z) {
  float *zp, x, y;
  zp = (float *)&z;

  /* first find out the large component */
  if (zp[0] > zp[1]) {
    x = zp[0];
    y = zp[1];
  } else {
    x = zp[1];
    y = zp[0];
  }

  return fabsf(x) * sqrtf(1.0f + y / x);
}

double
cabs(double _Complex z) {
  double *zp, x, y;
  zp = (double *)&z;

  /* first find out the large component */
  if (zp[0] > zp[1]) {
    x = zp[0];
    y = zp[1];
  } else {
    x = zp[1];
    y = zp[0];
  }

  return fabs(x) * sqrt(1.0f + y / x);
}

double
cimag(PLASMA_Complex64_t z) {
  return ((double *)&z)[1];
}
double
creal(PLASMA_Complex64_t z) {
  return ((double *)&z)[0];
}

PLASMA_Complex64_t
conj(PLASMA_Complex64_t z) {
  double *zp, *vp;
  PLASMA_Complex64_t v;

  zp = (double *)&z;
  vp = (double *)&v;
  vp[0] = zp[0];
  vp[1] = -zp[1];
  return v;
}

#ifdef __cplusplus
}
#endif

#endif /* PLASMA_HAS_COMPLEX */
