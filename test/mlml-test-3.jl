require(Pkg.dir("MethylUtils", "test","mlml-test-2.jl"))
include(Pkg.dir("MethylUtils", "src","Features.jl"))

#d = Dzeros()
#mlml_parallel_performance()

## TIMING
TIMING_PARALLEL = false
TIMING = false
PROFILING = false

PARALLEL_CYTOSINE_COUNTS = true


println("Get random counts")
@everywhere meth_counts = rand_bs_ox_counts(10,100,5)
println("Finished generating random counts")
#println(meth_counts)

if TIMING_PARALLEL
  d = dzeros(10)
  println("Start timing (parallel) EM")
  cgs = @time mlml_parallel_performance(d,meth_counts)
  type_cgs = typeof(cgs)
  println("Got $type_cgs")
  println(cgs)
end

if TIMING
  println("Start timing EM")
  cgs = @time mlml_performance(meth_counts)
  println("EM Done in time: $elapsed")
end

if PROFILING
    println("Start profiling EM")
    @profile mlml_performance(meth_counts)
    println ("Done profiling")
end



println("Get random CytosineCounts")
@everywhere bs_counts = rand_CytosineCounts(10,100,5)
@everywhere ox_counts = rand_CytosineCounts(10,100,5)

println("Finished generating random counts")
println(bs_counts)
println(ox_counts)

if PARALLEL_CYTOSINE_COUNTS
  d = dzeros(10)
  println("Start timing (parallel) EM")
  cgs = @time mlml_parallel_performance(d,bs_counts, ox_counts)
  type_cgs = typeof(cgs)
  println("Got $type_cgs")
  println(cgs)
end
