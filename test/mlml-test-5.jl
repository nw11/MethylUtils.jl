include(Pkg.dir("MethylUtils", "src","Features.jl"))
include(Pkg.dir("MethylUtils", "src","Functions.jl"))

using Distributions

function rand_CytosineCounts(num_cgs, max_cov, min_cov=0)
   coverage_distribution = [ x for x= 1:max_cov]
   cgs = Array(CytosineCount,num_cgs)
   k=1
   while k <= num_cgs
       (Cs, Ts) = rand_c_counts(0.8, coverage_distribution)
       if ( Cs + Ts ) >= min_cov
           cgs[k] = CytosineCount(SequencePos(1,2),true,int16(Cs) ,int16(Ts))
           k +=1
       end
   end
   return cgs
end

function rand_c_counts(prob, coverage_distribution)
    num_reads = sample(coverage_distribution)
    binom = Binomial(num_reads,prob)
    Cs = rand(binom)
    Ts = num_reads - Cs
    return (Cs,Ts)
end

SF_B = SequenceFeatures(Dict(),Dict())
for k=1:10
    println( "Making CytsonineCounts" )
    rc = rand_CytosineCounts(100,100,5)
    SF_B.features["chr$k"] = rc
end

SF_O = SequenceFeatures(Dict(),Dict())
for k=1:10
    println( "Making CytsonineCounts" )
    rc = rand_CytosineCounts(100,100,5)
    SF_O.features["chr$k"] = rc
end

meth_c = run_mlml(SF_B,SF_O)
println(meth_c)
