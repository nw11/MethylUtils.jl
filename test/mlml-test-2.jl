# This test file is mostly about performance
# and probably will not be FactTested

using ProfileView
using Distributions
include(Pkg.dir("MethylUtils", "src","mlml-bs-ox-typed.jl"))
include(Pkg.dir("MethylUtils", "src","Features.jl"))

function rand_CytosineCounts(num_cgs, max_cov, min_cov=0)
   coverage_distribution = [ x for x= 1:max_cov]
   cgs = Array(CytosineCount,num_cgs)
   k=1
   while k <= num_cgs
       (Cs, Ts) = rand_c_counts(0.8, coverage_distribution)
       if ( Cs + Ts ) >= min_cov
           cgs[k] = CytosineCount(GenomePos("chr1",1,2),true,int16(Cs) ,int16(Ts))
           k +=1
       end
   end
   return cgs
end


function rand_bs_ox_counts(num_cgs, max_cov, min_cov=0)
   coverage_distribution = [ x for x= 1:max_cov]
   cgs = Array(Array{Int16},num_cgs)
   k=1
   while k <= num_cgs
       (M_B, u_B) = rand_c_counts(0.8, coverage_distribution)
       (m_O, U_B) = rand_c_counts(0.77,coverage_distribution)

       if( ( M_B + u_B ) >= min_cov && (m_O + U_B ) >= min_cov )
           meth_counts = map( int16, [M_B, u_B, m_O, U_B ] )
           cgs[k] = meth_counts
           k += 1
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

function start_point(M_B,u_B,m_O,U_O)
        starting_estimate_p_m = M_B / (M_B + u_B)
        starting_estimate_p_h = 1 - starting_estimate_p_m
        return ( starting_estimate_p_m, starting_estimate_p_h )
end

function mlml_performance(meth_counts)
    num_cgs = length(meth_counts)
    new_array = Array(Float64, num_cgs)
    for k = 1:num_cgs
        (M_B,u_B,m_O,U_O) = meth_counts[k]
        (starting_estimate_p_m, starting_estimate_p_h) = start_point(M_B,u_B,m_O,U_O)
        # if we have overshoot then calculate
        #p_m_O = m_O / ( m_O + U_O )
        #if (starting_estimate_p_m - p_m_O) < 0
        (em_p_m, em_p_h, iter) = EM(M_B, u_B, m_O, U_O, starting_estimate_p_m, starting_estimate_p_h)
        new_array[k] = em_p_m
        #end
    end
    return new_array
end

# THIS WILL HAVE TO BE IN ANOTHER FILE TO REQUIRE
function mlml_parallel_performance(d::DArray, meth_counts)
   DArray(size(d),procs(d)) do I
         yrange = I[1]
         ymin=yrange[1]
         #xrange = I[2]
         #xmin=xrange[1]

         # HERE
         #new_array = Array(Float64, size(yrange,1) )
         new_array = Array(Array{Float64}, size(yrange,1))
         # here should be 2-d for both em_p_m and em_p_h and column for conflict or iteration
         # or rather than matrix perhaps just push on array of arrays?

         num_cgs = length(meth_counts)
         for y in yrange
            (M_B,u_B,m_O,U_O) = meth_counts[y]
            (starting_estimate_p_m, starting_estimate_p_h) = start_point(M_B,u_B,m_O,U_O)
            (em_p_m, em_p_h, iter) = EM(M_B, u_B, m_O, U_O, starting_estimate_p_m, starting_estimate_p_h)
            #put onto new array
            new_array[ y - ymin + 1 ] = [ em_p_m, em_p_h]
         end
          return new_array
    end
end


function mlml_parallel_performance(d::DArray, bs_counts , ox_counts )
   DArray(size(d),procs(d)) do I
         yrange = I[1]
         ymin=yrange[1]
         #xrange = I[2]
         #xmin=xrange[1]

         # HERE
         #new_array = Array(Float64, size(yrange,1) )
         len = size(yrange,1)
         println("size of array is $len")
         new_array = Array(CytosineMethylation, len)
         # here should be 2-d for both em_p_m and em_p_h and column for conflict or iteration
         # or rather than matrix perhaps just push on array of arrays?

         #num_cgs = length(bs_counts)
         for y in yrange
            (M_B,u_B) = ( bs_counts[y].C_count, bs_counts[y].T_count)
            (m_O,U_O) = ( ox_counts[y].C_count, ox_counts[y].T_count)
            (starting_estimate_p_m, starting_estimate_p_h) = start_point(M_B,u_B,m_O,U_O)
            (em_p_m, em_p_h, iter) = EM(M_B, u_B, m_O, U_O, starting_estimate_p_m, starting_estimate_p_h)
            #put onto new array
            new_array[ y - ymin + 1 ] = CytosineMethylation( bs_counts[y].pos, bs_counts[y].strand, em_p_m, em_p_h )
         end
         return new_array
    end
end
