module FeaturesTest

    using FactCheck
    using DataFrames
    include(Pkg.dir("MethylUtils", "src","Features.jl"))


   # GenomePosition

    g1 = GenomePos("chr1","3","15")
    g2 = GenomePos("chr1",3,15)
    g3 = GenomePos(1,3,15)

    facts("GenomicPosition") do
        @fact g1.chr   => int16(1)
        @fact g1.start => int16(3)
        @fact g1.stop  => int16(15)
        @fact g2.chr   => int16(1)
        @fact g2.start => int16(3)
        @fact g2.stop  => int16(15)
        @fact g3.chr   => int16(1)
    end

   # CytosineCount ?

   # DO SOME TESTS

   # FeatureTable

   df = DataFrame(name = ["A","B","C","D"], chr = 1:4, start = 1:4, stop = 1:4 )
   ft = FeatureTable( ["organism"=>"mm10"], df )

   # FeatureRegions
    fr = FeatureRegions(ft)

    println(fr)

    regions = getoverlaps(fr, 2, 2, 2)
    println(regions)
end
