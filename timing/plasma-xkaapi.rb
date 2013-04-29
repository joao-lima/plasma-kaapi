#!/usr/bin/env ruby

require 'optparse'
require 'pp'

# 1 to run, 0 to check output
$kaapi_run_test = 0
# Warmup
$kaapi_warmup= 0
# log file
$kaapi_log = STDOUT
# relative to HOME
$kaapi_install_dir = ENV['HOME'] + "/install/xkaapi/sc2013"
# results
$kaapi_result_dir = ENV['HOME'] + "/" + "res"
# output prefix 
$kaapi_output_prefix = ""
# version 
$kaapi_output_version = `date +%s`.chomp
# iterations
$niter = 2
# default environmental variables
$global_variables = [
    [ 'KAAPI_STACKSIZE_MASTER', '536870912' ],
#    [ 'KAAPI_CUDA_PEER',		'1' ],
    [ 'KAAPI_DISPLAY_PERF',      '1'        ]
  ]

def kaapi_get_map_cpuset(ncpu, ngpu)
  if ncpu == 4
    cpuset="0,5,6,11"
    gpuset="0~1,1~2,2~3,3~4,4~7,5~8,6~9,7~10"
  end
    
  if ncpu == 5
    cpuset="0,5,6,9,11"
    gpuset="0~1,1~2,2~3,3~4,4~7,5~8,7~10"
  end
  
  if ncpu == 6
    cpuset="0,1,5,6,10,11"
    gpuset="1~2,2~3,3~4,4~7,5~8,6~9"
  end
  
  if ncpu == 7
    cpuset="0,1,5,6,9,10,11"
    gpuset="1~2,2~3,3~4,4~7,5~8"
  end

  if ncpu == 8
    cpuset="0,2,4,5,6,8,10,11"
    gpuset="1~1,2~3,4~7,6~9"
  end

  if ncpu == 9
    cpuset="0,2,4,5,6,8,9,10,11"
    gpuset="1~1,2~3,4~7"
  end

  if ncpu == 10
    cpuset="0,2,3,4,5,6,8,9,10,11"
    gpuset="1~1,4~7"
  end

  if ncpu == 11
    cpuset="0,2,3,4,5,6,7,8,9,10,11"
    gpuset="1~1"
  end

  if ncpu == 12
    cpuset="0:11"
    gpuset=""
  end

  return cpuset, gpuset
end

def run_kaapi_program(program, n, nb, ib, ncpu, ngpu, win_size, cache, affinity, niter, foutput, ftrace)
  cmd = ""
  
  cmd = "LD_LIBRARY_PATH=" + $kaapi_install_dir + "/" + "lib"
  if ENV['LD_LIBRARY_PATH'] != nil
    cmd = cmd + ":" + ENV['LD_LIBRARY_PATH']
  end
  
  cmd = cmd + " " 
  
  $global_variables.each do |v|
    cmd = cmd +  v[0] + "=" +  v[1] + " "
  end
  
  cpuset, gpuset = kaapi_get_map_cpuset(ncpu, ngpu)
  
  cmd = cmd +
      "KAAPI_CPUSET="+ cpuset + " " +
      "KAAPI_GPUSET=" + gpuset + " " +
      "KAAPI_PUSH_AFFINITY=" + affinity + " " +
      "KAAPI_CUDA_WINDOW_SIZE=" + win_size.to_s  + " "
  
  if affinity == "heft" || affinity == "gdual"
    cmd = cmd + "KAAPI_PERFMODEL_CALIBRATE=1 KAAPI_PERFMODEL_TRANSFER=idgraf "
  end
  
  # PLASMA specific call
  range = n.to_s + ":" + n.to_s + ":" + n.to_s
  cmd = cmd + "./" + program['file'] + " --threads=1 --dyn " +
      "--n_range=" + range + " --niter=1 --nowarmup --ifmt=0 " +
      "--nb=" + nb.to_s + " --ib=" + ib.to_s
  
#  print cmd, "\n"
#  print $kaapi_result_dir + "/" + output + "\n"
#  print $kaapi_result_dir + "/" + trace + "\n"
#  $kaapi_log.print cmd + "\n"
  
  if affinity == "heft"
    if $kaapi_run_test == 1
      if $kaapi_warmup == 1
	result = `#{cmd} &> /dev/null`
      end
    end
  end
  niter.times do |i|
    $kaapi_log.print "(",$kaapi_output_version,",",ncpu,",",ngpu,") ","(",i,"/",niter,") ",program['name']," size=",n," block=",nb," ib=",ib," affinity=",affinity," prefix=",$kaapi_output_prefix,"\n"
    if $kaapi_run_test == 1
      result = `#{cmd}`
      foutput.print result.scan(/^#{program['name']}.*$/), "\n"
      ftrace.print result
      #puts result.sub!(/^(#{testname}.*)$/,"")  
    end
  end
end

def test_open_file(name)
#  $kaapi_log.print name,"\n"
  if $kaapi_run_test == 1
    f = File.open($kaapi_result_dir + "/" + name, "w")
    return f
  else
    return STDOUT
  end
end

def test_close_file(fs)
  if $kaapi_run_test == 1
    fs.close
  end
end

def run_weak_scalability(program)
  ninputs = [ 8192 ]
#  ninputs = [ 8192, 16384 ]
  nblocks = program['nb']
  nib = program['ib']
  ncpus = (4..11).to_a
  ngpu = 8
  win_size = 2
  cache = "lru_double"
#  naffinity= [ "default", "writer", "heft" ]
  naffinity= [ "writer", "heft", "gdual" ]

  if $kaapi_output_prefix != ""
    outprefix = "-" + $kaapi_output_prefix + "-"
  else
    outprefix = "-"
  end

  ncpus.each do |ncpu|
    ngpu = 12 - ncpu
    naffinity.each do |aff|
      output = "kaapi" + outprefix + program['name'] + "-" + aff + "-" + ncpu.to_s + "cpu" + ngpu.to_s + "gpu" +
	"-" + $kaapi_output_version + ".csv"
      foutput = test_open_file(output)
      ninputs.each do |n|
	nblocks.each do |b|
	  nib.each do |ib|
	    trace = "perfcounter" + n.to_s + "size" + b.to_s + "nb" + ib.to_s + "ib-" + output
	    ftrace = test_open_file(trace)
	    run_kaapi_program(program, n, b, ib, ncpu, ngpu, win_size, cache, aff, $niter, foutput, ftrace)
	    test_close_file(ftrace)
	  end
	end
      end
      test_close_file(foutput)
    end
  end
end

def run_strong_scalability(program)
  ninputs = (4096..20480).step(4096).to_a
  nblocks = program['nb']
  nib = program['ib']
  ncpus = [ 4 ]
  ngpu = 8
  win_size = 2
  cache = "lru_double"
  naffinity= [ "default", "writer", "heft" ]

  if $kaapi_output_prefix != ""
    outprefix = "-" + $kaapi_output_prefix + "-"
  else
    outprefix = "-"
  end

  ncpus.each do |ncpu|
    ngpu = 12 - ncpu
    naffinity.each do |aff|
      output = "kaapi" + outprefix + program['name'] + "-" + aff + "-" + ncpu.to_s + "cpu" + ngpu.to_s + "gpu" +
	"-" + $kaapi_output_version + ".csv"
      foutput = test_open_file(output)
      ninputs.each do |n|
	nblocks.each do |b|
	  nib.each do |ib|
	    trace = "perfcounter" + n.to_s + "size" + b.to_s + "nb" + ib.to_s + "ib-" + output
	    ftrace = test_open_file(trace)
	    run_kaapi_program(program, n, b, ib, ncpu, ngpu, win_size, cache, aff, $niter, foutput, ftrace)
	    test_close_file(ftrace)
	  end
	end
      end
      test_close_file(foutput)
    end
  end
end

def run_experiments()
  ib = [128]
  nb = [512]
  nprograms = [
    { 'file' => 'time_dpotrf_tile.xkaapi',	    'name' => 'PLASMA_dpotrf_Tile',	    'nb' => [800], 'ib' => ib },
    { 'file' => 'time_dgemm_tile.xkaapi',	    'name' => 'PLASMA_dgemm_Tile',	    'nb' => [1024], 'ib' => ib },
    { 'file' => 'time_dgetrf_incpiv_tile.xkaapi',  'name' => 'PLASMA_dgetrf_incpiv_Tile',  'nb' => [800], 'ib' => ib },
    { 'file' => 'time_dgeqrf_tile.xkaapi',	    'name' => 'PLASMA_dgeqrf_Tile',	    'nb' => nb, 'ib' => ib }
  #    { 'file' => 'time_dgeqrfrh_tile',	    'name' => 'PLASMA_dgeqrfrh_Tile',	    'nb' => nb, 'ib' => ib }
  ]

#  p nprograms

  nprograms.each do |prog|
    run_weak_scalability(prog)
#    run_strong_scalability(prog)
  end
end

options = {}
#options[:method] = []

OptionParser.new do |opts|
    opts.banner = "Usage: script.rb <options>"
    opts.separator ""

    opts.on("-r", "--run", "Run PLASMA benchmarks.") do 
      $kaapi_run_test = 1
    end

    opts.on("-w", "--warmup", "Warmup at each execution set.") do 
      $kaapi_warmup= 1
    end

    opts.on("-i", "--niter [n]", Integer, "Number of iterations.") do |n|
     $niter = n
    end

    opts.on("-p", "--prefix [name]", String, "Output variation.") do |s|
     $kaapi_output_prefix = s
    end

end.parse!

run_experiments
