#require(Pkg.dir("MethylUtils", "test","mlml-test-2.jl"))
#include(Pkg.dir("MethylUtils", "src","Features.jl"))

include(Pkg.dir("MethylUtils", "test","mlml-test-2.jl"))
include(Pkg.dir("MethylUtils", "src","mlml-bs-ox-typed.jl"))

require(Pkg.dir("MethylUtils", "src","Features.jl"))
require(Pkg.dir("MethylUtils", "src","mlml-bs-ox-typed.jl"))

println("Get random counts")
meth_counts = rand_bs_ox_counts(10,100,5)

bs_counts = rand_CytosineCounts(1000000,100,5)
ox_counts = rand_CytosineCounts(1000000,100,5)
println("Finished generating random counts")

#println(meth_counts)
#meth_counts
#println(bs_counts)
#println(ox_counts)
#@time a = pmap( EM,  bs_counts,ox_counts, 1:length(bs_counts) )
#println(a)
#println(length(a))

println("run mlml_run" )
@time run_mlml(bs_counts, ox_counts)
