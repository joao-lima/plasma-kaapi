/**
 *
 * @file time_main.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author ???
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 **/

#if defined( _WIN32 ) || defined( _WIN64 )
#define int64_t __int64
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef PLASMA_EZTRACE
#include <eztrace.h>
#endif

#if defined( _WIN32 ) || defined( _WIN64 )
#include <windows.h>
#include <time.h>
#include <sys/timeb.h>
#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif

struct timezone
{
    int  tz_minuteswest; /* minutes W of Greenwich */
    int  tz_dsttime;     /* type of dst correction */
};

int gettimeofday(struct timeval* tv, struct timezone* tz)
{
    FILETIME         ft;
    unsigned __int64 tmpres = 0;
    static int       tzflag;

    if (NULL != tv)
        {
            GetSystemTimeAsFileTime(&ft);
            tmpres |=  ft.dwHighDateTime;
            tmpres <<= 32;
            tmpres |=  ft.dwLowDateTime;

            /*converting file time to unix epoch*/
            tmpres /= 10;  /*convert into microseconds*/
            tmpres -= DELTA_EPOCH_IN_MICROSECS;

            tv->tv_sec  = (long)(tmpres / 1000000UL);
            tv->tv_usec = (long)(tmpres % 1000000UL);
        }
    if (NULL != tz)
        {
            if (!tzflag)
                {
                    _tzset();
                    tzflag++;
                }
            tz->tz_minuteswest = _timezone / 60;
            tz->tz_dsttime     = _daylight;
        }
    return 0;
}

#else  /* Non-Windows */
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <cblas.h>
#include <lapacke.h>
#include <plasma.h>
#include <core_blas.h>

#if defined(PLASMA_CUDA)
#include <cuda_runtime_api.h>
#endif

#if defined(PLASMA_KAAPI)
#include "kaapi.h"
#include <core_kblas.h>
#endif

#include <plasma_tmg.h>
#include "flops.h"
#include "timing.h"
#include "auxiliary.h"

static int RunTest(int *iparam, _PREC *dparam, real_Double_t *t_);

real_Double_t cWtime(void);

int ISEED[4] = {0,0,0,1};   /* initial seed for zlarnv() */

/*
 * struct timeval {time_t tv_sec; suseconds_t tv_usec;};
 */
real_Double_t cWtime(void)
{
    struct timeval tp;
    gettimeofday( &tp, NULL );
    return tp.tv_sec + 1e-6 * tp.tv_usec;
}

static int
Test(int64_t n, int *iparam) {
    int           i, j, iter;
    int thrdnbr, niter, nrhs, m;
    real_Double_t *t, gflops;
    _PREC         eps = _LAMCH( 'e' );
    _PREC         dparam[TIMING_DNBPARAM];
    real_Double_t fmuls, fadds, fp_per_mul, fp_per_add;
    double        sumgf, sumgf2, sumt, sd;
    char         *s;
    char         *env[] = {
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "GOTO_NUM_THREADS",
        "ACML_NUM_THREADS",
        "ATLAS_NUM_THREADS",
        "BLAS_NUM_THREADS", ""
    };
    int gnuplot = 0;

    memset( &dparam, 0, TIMING_DNBPARAM * sizeof(_PREC) );

    thrdnbr = iparam[TIMING_THRDNBR];
    niter   = iparam[TIMING_NITER];
    nrhs    = iparam[TIMING_NRHS];
    m       = iparam[TIMING_M];
    
       if (n < 0 || thrdnbr < 0) {
#if 0
        if  (iparam[TIMING_CHECK] )
            printf( "#   N NRHS threads seconds   Gflop/s Deviation    ||Ax-b||       ||A||       ||x||       ||b||"
                    "         eps ||Ax-b||/N/eps/(||A||||x||+||b||)\n" );
        else
            printf( "#N,NB,IB,NRHS,threads,seconds,Gflop/s,Deviation\n" );
#endif
//            printf( "#   N NRHS threads seconds   Gflop/s Deviation\n" );

        if (gnuplot) {
            printf( "set title '%d_NUM_THREADS: ", thrdnbr );
            for (i = 0; env[i][0]; ++i) {
                s = getenv( env[i] );

                if (i) printf( " " ); /* separating space */

                for (j = 0; j < 5 && env[i][j] && env[i][j] != '_'; ++j)
                    printf( "%c", env[i][j] );

                if (s)
                    printf( "=%s", s );
                else
                    printf( "->%s", "?" );
            }
            printf( "'\n" );
            printf( "%s\n%s\n%s\n%s\n%s%s%s\n",
                    "set xlabel 'Matrix size'",
                    "set ylabel 'Gflop/s'",
                    "set key bottom",
                    gnuplot > 1 ? "set terminal png giant\nset output 'timeplot.png'" : "",
                    "plot '-' using 1:5 title '", _NAME, "' with linespoints" );
        }

        return 0;
    }

#if 0
    printf( "%d,%d,%d,%d,%d,", iparam[TIMING_N], iparam[TIMING_NB], iparam[TIMING_IB], iparam[TIMING_NRHS], iparam[TIMING_THRDNBR] );
    fflush( stdout );
#endif

    t = (real_Double_t*)malloc(niter*sizeof(real_Double_t));

    if (sizeof(_TYPE) == sizeof(_PREC)) {
        fp_per_mul = 1;
        fp_per_add = 1;
    } else {
        fp_per_mul = 6;
        fp_per_add = 2;
    }

    fadds = (double)(_FADDS);
    fmuls = (double)(_FMULS);
    gflops = 0.0;

    if ( iparam[TIMING_WARMUP] ) {
        RunTest( iparam, dparam, &(t[0]));
    }

    sumgf  = 0.0;
    sumgf2 = 0.0;
    sumt   = 0.0;
    
    for (iter = 0; iter < niter; iter++)
    {
#ifdef PLASMA_EZTRACE
        if( iter == 0 ) {
            eztrace_start();
            RunTest( iparam, dparam, &(t[iter]));
            eztrace_stop();
        }
        else
#endif
            RunTest( iparam, dparam, &(t[iter]));
        
        gflops = 1e-9 * (fmuls * fp_per_mul + fadds * fp_per_add) / t[iter];
        sumt   += t[iter];
        sumgf  += gflops;
        sumgf2 += gflops*gflops;
    }

    gflops = sumgf/niter;
    sd = sqrt((sumgf2 - (sumgf*sumgf)/niter)/niter);

    if ( iparam[TIMING_CHECK] )
    {
        printf( "%s,%.10f,%.10f,%.10f,%8.5e,%8.5e,%8.5e,%8.5e,%8.5e,%8.5e\n",
               _NAME,
                sumt/niter, gflops, sd, dparam[TIMING_RES], dparam[TIMING_ANORM], dparam[TIMING_XNORM], dparam[TIMING_BNORM], eps, 
                dparam[TIMING_RES] / n / eps / (dparam[TIMING_ANORM] * dparam[TIMING_XNORM] + dparam[TIMING_BNORM] ));
    }
    else
    {
      const char* sched = getenv("KAAPI_SCHED");
      const char* cpuset = getenv("KAAPI_CPUSET");
      const char* gpuset = getenv("KAAPI_GPUSET");
      const char* ldpath = getenv("LD_LIBRARY_PATH");
      const char* priority = getenv("QUARK_PRIORITY");
      const char* cuda_window = getenv("KAAPI_CUDA_WINDOW_SIZE");
      double alpha= 1.f;
      double beta = 1.f;
      int steal= 0;
      int calibrate= 0;
      int setarch=0;
      int prio= 0;

      if(getenv("KAAPI_PERFMODEL_CALIBRATE") != 0)
      {
        calibrate = 1;
      }

      if(getenv("QUARK_ARCH_GPU_ONLY") != 0)
      {
        setarch = 1;
      }
      
      if(priority != 0)
      {
        prio = atoi(priority);
      }
      
      if(getenv("KAAPI_HEFT_STEAL") != 0)
      {
        steal = 1;
      }
      if(getenv("KAAPI_GDUAL_STEAL") != 0)
      {
        steal = 1;
      }      
      
      if(getenv("KAAPI_HEFT_ALPHA") != 0)
      {
        const char* str = getenv("KAAPI_HEFT_ALPHA");
        alpha = atof(str);
      }
      if(getenv("KAAPI_HEFT_BETA") != 0)
      {
        const char* str = getenv("KAAPI_HEFT_BETA");
        beta = atof(str);
      }
      if(getenv("KAAPI_GDUAL_ALPHA") != 0)
      {
        const char* str = getenv("KAAPI_GDUAL_ALPHA");
        alpha = atof(str);
      }
      if(getenv("KAAPI_GDUAL_BETA") != 0)
      {
        const char* str = getenv("KAAPI_GDUAL_BETA");
        beta = atof(str);        
      }
//      fprintf(stdout, "name setarch priority cpu gpu CUDA_WINDOW scheduler alpha beta calibrate steal(heft/gdual) N NB IB NRHS THREADS(PLASMA) time Gflops SD  cpuset gpuset LD_LIBRARY_PATH\n");
        fprintf(stdout, "%s;%d;%d;%d;%d;%s;%s;%.4f;%.4f;%d;%d;%d;%d;%d;%d;%d;%.10f;%.10f;%.10f;%s;%s;%s\n",
                _NAME, setarch, prio, kaapi_getconcurrency_cpu(), kaapi_getconcurrency_gpu(), cuda_window, sched, alpha, beta, calibrate, steal,
                iparam[TIMING_N], iparam[TIMING_NB], iparam[TIMING_IB], iparam[TIMING_NRHS], iparam[TIMING_THRDNBR] ,sumt/niter, gflops, sd, cpuset, gpuset, ldpath );
    }
    fflush( stdout );
    free(t);

    return 0;
}

static int
startswith(const char *s, const char *prefix) {
    size_t n = strlen( prefix );
    if (strncmp( s, prefix, n ))
        return 0;
    return 1;
}

static int
get_range(char *range, int *start_p, int *stop_p, int *step_p) {
    char *s, *s1, buf[21];
    int colon_count, copy_len, nbuf=20, n;
    int start=1000, stop=10000, step=1000;

    colon_count = 0;
    for (s = strchr( range, ':'); s; s = strchr( s+1, ':'))
        colon_count++;

    if (colon_count == 0) { /* No colon in range. */
        if (sscanf( range, "%d", &start ) < 1 || start < 1)
            return -1;
        step = start / 10;
        if (step < 1) step = 1;
        stop = start + 10 * step;

    } else if (colon_count == 1) { /* One colon in range.*/
        /* First, get the second number (after colon): the stop value. */
        s = strchr( range, ':' );
        if (sscanf( s+1, "%d", &stop ) < 1 || stop < 1)
            return -1;

        /* Next, get the first number (before colon): the start value. */
        n = s - range;
        copy_len = n > nbuf ? nbuf : n;
        strncpy( buf, range, copy_len );
        buf[copy_len] = 0;
        if (sscanf( buf, "%d", &start ) < 1 || start > stop || start < 1)
            return -1;

        /* Let's have 10 steps or less. */
        step = (stop - start) / 10;
        if (step < 1)
            step = 1;
    } else if (colon_count == 2) { /* Two colons in range. */
        /* First, get the first number (before the first colon): the start value. */
        s = strchr( range, ':' );
        n = s - range;
        copy_len = n > nbuf ? nbuf : n;
        strncpy( buf, range, copy_len );
        buf[copy_len] = 0;
        if (sscanf( buf, "%d", &start ) < 1 || start < 1)
            return -1;

        /* Next, get the second number (after the first colon): the stop value. */
        s1 = strchr( s+1, ':' );
        n = s1 - (s + 1);
        copy_len = n > nbuf ? nbuf : n;
        strncpy( buf, s+1, copy_len );
        buf[copy_len] = 0;
        if (sscanf( buf, "%d", &stop ) < 1 || stop < start)
            return -1;

        /* Finally, get the third number (after the second colon): the step value. */
        if (sscanf( s1+1, "%d", &step ) < 1 || step < 1)
            return -1;
    } else

        return -1;

    *start_p = start;
    *stop_p = stop;
    *step_p = step;

    return 0;
}

static void
show_help(char *prog_name) {
    printf( "Usage:\n%s [options]\n\n", prog_name );
    printf( "Options are:\n" );
    printf( "  --threads=C    Number of threads (default: 1)\n" );
    printf( "  --n_range=R    Range of N values: Start:Stop:Step (default: 500:5000:500)\n" );
    //    printf( "  --gnuplot      produce output suitable for gnuplot" );
    printf( "  --[no]check    Check result (default: nocheck)\n" );
    printf( "  --[no]warmup   Perform a warmup run to pre-load libraries (default: warmup)\n");
    printf( "  --niter=N      Number of iterations (default: 1)\n");
    printf( "  --nb=N         Nb size. Not used if autotuning is activated (default: 128)\n");
    printf( "  --ib=N         IB size. Not used if autotuning is activated (default: 32)\n");
    printf( "  --[no]dyn      Activate Dynamic scheduling (default: nodyn)\n");
    printf( "  --[no]atun     Activate autotuning (default: noatun)\n");
    printf( "  --ifmt         Input format. 0: CM, 1: CCRB, 2: CRRB, 3: RCRB, 4: RRRB, 5: RM (default: 0)\n");
    printf( "  --ofmt         Output format. 0: CM, 1: CCRB, 2: CRRB, 3: RCRB, 4: RRRB, 5: RM (default: 1)\n");
    printf( "  --thrdbypb     Number of threads per subproblem for inplace transformation (default: 1)\n");
}
static void
get_thread_count(int *thrdnbr) {
#if defined WIN32 || defined WIN64
    sscanf( getenv( "NUMBER_OF_PROCESSORS" ), "%d", thrdnbr );
#else
    *thrdnbr = sysconf(_SC_NPROCESSORS_ONLN);
#endif
}

int
main(int argc, char *argv[]) {
    int i, m;
    int start =  500;
    int stop  = 5000;
    int step  =  500;
    int iparam[TIMING_INBPARAM];

    iparam[TIMING_CHECK         ] = 0;
    iparam[TIMING_WARMUP        ] = 1;
    iparam[TIMING_NITER         ] = 1;
    iparam[TIMING_M             ] = -1;
    iparam[TIMING_N             ] = 500;
    iparam[TIMING_NB            ] = 128;
    iparam[TIMING_IB            ] = 32;
    iparam[TIMING_NRHS          ] = 1;
    iparam[TIMING_THRDNBR       ] = 1;
    iparam[TIMING_THRDNBR_SUBGRP] = 1;
    iparam[TIMING_SCHEDULER     ] = 0;
    iparam[TIMING_AUTOTUNING    ] = 1;
    iparam[TIMING_INPUTFMT      ] = 0;
    iparam[TIMING_OUTPUTFMT     ] = 0;

    get_thread_count( &(iparam[TIMING_THRDNBR]) );

    for (i = 1; i < argc && argv[i]; ++i) {
        if (startswith( argv[i], "--help" )) {
            show_help( argv[0] );
            return EXIT_SUCCESS;
        } else if (startswith( argv[i], "--n_range=" )) {
            get_range( strchr( argv[i], '=' ) + 1, &start, &stop, &step );
        } else if (startswith( argv[i], "--threads=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[TIMING_THRDNBR]) );
        /* } else if (startswith( argv[i], "--gnuplot-png" )) { */
        /*     gnuplot = 2; */
        /* } else if (startswith( argv[i], "--gnuplot" )) { */
        /*     gnuplot = 1; */
        } else if (startswith( argv[i], "--check" )) {
            iparam[TIMING_CHECK] = 1;
        } else if (startswith( argv[i], "--nocheck" )) {
            iparam[TIMING_CHECK] = 0;
        } else if (startswith( argv[i], "--warmup" )) {
            iparam[TIMING_WARMUP] = 1;
        } else if (startswith( argv[i], "--nowarmup" )) {
            iparam[TIMING_WARMUP] = 0;
        } else if (startswith( argv[i], "--dyn" )) {
            iparam[TIMING_SCHEDULER] = 1;
        } else if (startswith( argv[i], "--nodyn" )) {
            iparam[TIMING_SCHEDULER] = 0;
        } else if (startswith( argv[i], "--atun" )) {
            iparam[TIMING_AUTOTUNING] = 1;
        } else if (startswith( argv[i], "--noatun" )) {
            iparam[TIMING_AUTOTUNING] = 0;
        } else if (startswith( argv[i], "--m=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[TIMING_M]) );
        } else if (startswith( argv[i], "--nb=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[TIMING_NB]) );
        } else if (startswith( argv[i], "--nrhs=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[TIMING_NRHS]) );
        } else if (startswith( argv[i], "--ib=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[TIMING_IB]) );
        } else if (startswith( argv[i], "--ifmt=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[TIMING_INPUTFMT]) );
        } else if (startswith( argv[i], "--ofmt=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[TIMING_OUTPUTFMT]) );
        } else if (startswith( argv[i], "--thrdbypb=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[TIMING_THRDNBR_SUBGRP]) );
        } else if (startswith( argv[i], "--niter=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &iparam[TIMING_NITER] );
        } else {
            fprintf( stderr, "Unknown option: %s\n", argv[i] );
        }
    }

    m = iparam[TIMING_M];

    if (step < 1) step = 1;

    Test( -1, iparam ); /* print header */
    for (i = start; i <= stop; i += step)
    {        
        if ( m == -1 )
            iparam[TIMING_M] = i;
        iparam[TIMING_N] = i;
        Test( i, iparam );
    }

    /* if (gnuplot) { */
    /*         printf( "%s\n%s\n", */
    /*                 "e", */
    /*                 gnuplot > 1 ? "" : "pause 10" ); */
    /* } */

    return EXIT_SUCCESS;
}
