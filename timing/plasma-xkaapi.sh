#!/bin/bash

export LD_LIBRARY_PATH=$HOME/install/xkaapi/default/lib:$LD_LIBRARY_PATH
version="$(date +%s)"

dorun="yes"
ALLITER=10

#verif="--check"
#debug=1

#export KAAPI_RECORD_TRACE=1
#export KAAPI_RECORD_MASK="COMPUTE,IDLE"

# export KAAPI_DUMP_GRAPH=1
# export KAAPI_DOT_NOVERSION_LINK=1
# export KAAPI_DOT_NODATA_LINK=1
# export KAAPI_DOT_NOLABEL_INFO=1
# export KAAPI_DOT_NOACTIVATION_LINK=1
# export KAAPI_DOT_NOLABEL=1

export KAAPI_DISPLAY_PERF=1

#export KAAPI_PUSH_AFFINITY="writer"
#export KAAPI_PUSH_AFFINITY="heft"
export KAAPI_HEFT_CALIBRATE=1
#    export KAAPI_STEAL_AFFINITY="writer"
#    export KAAPI_PUSH_AFFINITY="locality"
#    export KAAPI_STEAL_AFFINITY="locality"

#ib=1024
#ib=512
#ib=256
ib=128
#ib=64

options=""

if [ -n "$ib" ]
then
  options="$options --ib=$ib"
fi

function run_test {
#  export KAAPI_CPUSET="0:11"
  export KAAPI_CPUSET="0,5,6,11"
  export KAAPI_GPUSET="0~1"
#  export KAAPI_CPUSET="0,5,6,11"
#  export KAAPI_GPUSET="0~1,1~2,2~3,3~4,4~7,5~8,6~9,7~10"

  export KAAPI_WINDOW_SIZE=2

  #testing="time_dpotrf"
  #testing="time_dpotrf_tile"
  #testing="time_dgetrf_tile"
  #testing="time_dgetrf_incpiv_tile"
  #testing="time_dgetrf"
  #testing="time_dgetrf_incpiv"
#  testing="time_dgemm_tile"
  #testing="time_sgemm_tile"
  #testing="time_dpotrf_tile time_dgemm_tile"
#  testing="time_dgeqrf_tile"
  testing="time_dgeqrfrh_tile"

  #cache_policy="lru_double lru"
  cache_policy="lru_double"
  #cache_policy="lru"

  #range="20480:20480:20480"
  #range="32768:32768:32768"
  #range="40960:40960:40960"
#  range="10240:10240:10240"
  #range="10240:10240:20480"
  range="8192:8192:8192"
  #range="2048:8192:2048"
#  range="2048:2048:2048"
  #nblocks="2048"
#  nblocks="1024"
  nblocks="512"
  #nblocks=" $(seq 400 100 1200)"
  niter=1

  for test in $testing
  do
    for cache in $cache_policy
    do
      for nb in $nblocks
      do
	echo "(ncpu=$ncpu,ngpu=$ngpu,niter=$niter) KAAPI_GPU_CACHE_POLICY=$cache ./$test --threads=1 --dyn --n_range=$range --nb=$nb \
		  --niter=$niter --nowarmup --ifmt=0 $options $verif"

	if [ -n "$debug" ]
	then
	  echo "debug"
	  KAAPI_STACKSIZE_MASTER=536870912 \
		    gdb ./$test 

	else
    #     KAAPI_STACKSIZE_MASTER=1073741824 
	  for i in `seq 1 $niter`
	  do
	     KAAPI_STACKSIZE_MASTER=536870912 \
	     KAAPI_GPU_CACHE_POLICY=$cache \
		      ./$test --threads=1 --dyn --n_range=$range --niter=1 \
		      --nowarmup --ifmt=0  --nb=$nb $options $verif 
	   done
	fi
      done
    done
  done
}

#run_test
#exit 0

function generic_plasma {
  ncpu="$1"
  cpuset="$2"
  ngpu="$3"
  gpuset="$4"
  testing="$5"
  cache_policy="$6"
  range="$7"
  nblocks="$8"
  niter="$9"
  affinity="${10}"

  export KAAPI_CPUSET="$cpuset"
  export KAAPI_GPUSET="$gpuset"
  export KAAPI_CUDA_WINDOW_SIZE=2

  for aff in $affinity
  do
    for test in $testing
    do
      for cache in $cache_policy
      do
	for nb in $nblocks
	do
	  output="$HOME/res/plasma-${test}-${aff}-${ncpu}cpu${ngpu}gpu-${version}.csv"
#	output="$HOME/res/plasma-${test}-blocktest-${ncpu}cpu${ngpu}gpu-${version}.csv"

	  echo "(ncpu=$ncpu,ngpu=$ngpu,niter=$niter,$aff) ./$test --threads=1 --dyn --n_range=$range --nb=$nb --niter=$niter --nowarmup --ifmt=0 $options $verif $output"

	  if [ -n "$debug" ]
	  then
	    echo "debug"
	    KAAPI_STACKSIZE_MASTER=536870912 \
		      gdb ./$test 

	  else
      #     KAAPI_STACKSIZE_MASTER=1073741824 
	    for i in `seq 1 $niter`
	    do
	      if [ -n "$dorun" ]
	      then
		 KAAPI_STACKSIZE_MASTER=536870912 \
		 KAAPI_GPU_CACHE_POLICY=$cache \
		 KAAPI_PUSH_AFFINITY=$aff \
			  ./$test --threads=1 --dyn --n_range=$range --niter=1 \
			  --nowarmup --ifmt=0  --nb=$nb $options $verif  &>> $output
	       fi
	     done
	  fi
	done
      done
    done
  done
}

function hybrid_weak {
    #testing="time_dpotrf"
    #testing="time_dpotrf_tile"
    #testing="time_dgetrf_tile"
    #testing="time_dgetrf_incpiv_tile"
    #testing="time_dgetrf"
    #testing="time_dgetrf_incpiv"
#    testing="time_dgemm_tile"
#    testing="time_dgeqrf_tile"
    #testing="time_sgemm_tile"
    #testing="time_dpotrf_tile time_dgemm_tile"
    testing="
    time_dpotrf_tile time_dgetrf_incpiv_tile time_dgeqrf_tile time_dgeqrfrh_tile time_dgemm_tile
    "

    #cache_policy="lru_double lru"
    cache_policy="lru_double"
    #cache_policy="lru"

    affinity="default writer heft"

    #range="20480:20480:20480"
    #range="32768:32768:32768"
    #range="40960:40960:40960"
    range="10240:10240:10240"
    #range="10240:10240:20480"
    #range="8192:8192:8192"
    #range="2048:8192:2048"
    #range="2048:2048:2048"
    #nblocks="2048"
    #nblocks="1024"
    nblocks="512"
#    nblocks=" $(seq 400 100 1200)"
    niter=$ALLITER

    ncpu=4
    cpuset="0,5,6,11"
    export KAAPI_NCPU=$ncpu
    ngpu=8
    export KAAPI_NGPU=$ngpu
    gpuset="0~1,1~2,2~3,3~4,4~7,5~8,6~9,7~10"
    generic_plasma "$ncpu" "$cpuset" "$ngpu" "$gpuset" "$testing" "$cache_policy" "$range" "$nblocks" "$niter" "$affinity"

    ncpu=5
    cpuset="0,5,6,9,11"
    export KAAPI_NCPU=$ncpu
    ngpu=7
    export KAAPI_NGPU=$ngpu
    gpuset="0~1,1~2,2~3,3~4,4~7,5~8,7~10"
    generic_plasma "$ncpu" "$cpuset" "$ngpu" "$gpuset" "$testing" "$cache_policy" "$range" "$nblocks" "$niter" "$affinity"

    ncpu=6
    cpuset="0,1,5,6,10,11"
    export KAAPI_NCPU=$ncpu
    ngpu=6
    export KAAPI_NGPU=$ngpu
    gpuset="1~2,2~3,3~4,4~7,5~8,6~9"
    generic_plasma "$ncpu" "$cpuset" "$ngpu" "$gpuset" "$testing" "$cache_policy" "$range" "$nblocks" "$niter" "$affinity"

    ncpu=7
    cpuset="0,1,5,6,9,10,11"
    export KAAPI_NCPU=$ncpu
    ngpu=5
    export KAAPI_NGPU=$ngpu
    gpuset="1~2,2~3,3~4,4~7,5~8"
    generic_plasma "$ncpu" "$cpuset" "$ngpu" "$gpuset" "$testing" "$cache_policy" "$range" "$nblocks" "$niter" "$affinity"

    ncpu=8
    cpuset="0,2,4,5,6,8,10,11"
    export KAAPI_NCPU=$ncpu
    ngpu=4
    export KAAPI_NGPU=$ngpu
    gpuset="1~1,2~3,4~7,6~9"
    generic_plasma "$ncpu" "$cpuset" "$ngpu" "$gpuset" "$testing" "$cache_policy" "$range" "$nblocks" "$niter" "$affinity"

    ncpu=9
    cpuset="0,2,4,5,6,8,9,10,11"
    export KAAPI_NCPU=$ncpu
    ngpu=3
    export KAAPI_NGPU=$ngpu
    gpuset="1~1,2~3,4~7"
    generic_plasma "$ncpu" "$cpuset" "$ngpu" "$gpuset" "$testing" "$cache_policy" "$range" "$nblocks" "$niter" "$affinity"

    ncpu=10
    cpuset="0,2,3,4,5,6,8,9,10,11"
    export KAAPI_NCPU=$ncpu
    ngpu=2
    export KAAPI_NGPU=$ngpu
    gpuset="1~1,4~7"
    generic_plasma "$ncpu" "$cpuset" "$ngpu" "$gpuset" "$testing" "$cache_policy" "$range" "$nblocks" "$niter" "$affinity"

    ncpu=11
    cpuset="0,2,3,4,5,6,7,8,9,10,11"
    export KAAPI_NCPU=$ncpu
    ngpu=1
    export KAAPI_NGPU=$ngpu
    gpuset="1~1"
    generic_plasma "$ncpu" "$cpuset" "$ngpu" "$gpuset" "$testing" "$cache_policy" "$range" "$nblocks" "$niter" "$affinity"
}

function hybrid_strong {
    #testing="time_dpotrf"
    #testing="time_dpotrf_tile"
    #testing="time_dgetrf_tile"
    #testing="time_dgetrf_incpiv_tile"
    #testing="time_dgetrf"
    #testing="time_dgetrf_incpiv"
#    testing="time_dgemm_tile"
#    testing="time_dgeqrf_tile"
    #testing="time_sgemm_tile"
    testing="
    time_dpotrf_tile time_dgetrf_incpiv_tile time_dgeqrf_tile time_dgeqrfrh_tile time_dgemm_tile
    "

    #cache_policy="lru_double lru"
    cache_policy="lru_double"
    #cache_policy="lru"

    #range="20480:20480:20480"
    #range="32768:32768:32768"
    #range="40960:40960:40960"
    #range="10240:10240:10240"
    #range="10240:10240:20480"
    #range="8192:8192:8192"
    #range="2048:8192:2048"
    #range="2048:2048:2048"
    #nblocks="2048"
    #nblocks="1024"
#    nblocks="512 1024"
    nblocks="512"
#    nblocks=" $(seq 400 100 1200)"
    niter=$ALLITER
#    ninputs="$(seq 4096 4096 20480)"
#    ninputs="$(seq 4096 4096 16384)"
    ninputs="10240"

    ncpu=4
    cpuset="0,5,6,11"
    ngpu=8
    gpuset="0~1,1~2,2~3,3~4,4~7,5~8,6~9,7~10"

    for n in $ninputs
    do
      range="$n:$n:$n"
      generic_plasma "$ncpu" "$cpuset" "$ngpu" "$gpuset" "$testing" "$cache_policy" "$range" "$nblocks" "$niter"
    done
}

hybrid_weak
#hybrid_strong
