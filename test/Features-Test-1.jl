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

   df = DataFrame(name = ["A","B","C","D"], chr = ["chr1","chr2","chr1", "chr3"], start = 0:3, stop = 0:3 )
   ft = FeatureTable( ["organism"=>"mm10"], df )

   # FeatureRegions
   fr = FeatureRegions(ft)
   println(fr)
   regions = getoverlaps(fr, "chr2", 2, 2)
   println(regions)

   # load ensGene rows into DF
   # names_str = [ "bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart",
   #            "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2",
   ##           "cdsStartStat", "cdsEndStat", "exonFrames" ]
   # names = map(symbol, names_str)
   #DF = readtable( Pkg.dir("MethylUtils","testdata","ensGene1000lines.txt.gz"),
   #                separator='\t',
   #                names=names
   #                )
   # make a feature table
   #FT = FeatureTable({"organism" => "superman" }, DF)
   #FR = FeatureRegions(FT)

   # SequenceFeatures
   SF = SequenceFeatures( Dict(), Dict())
   feature = CytosineCount( SequencePos(100,200),true,15,90)
   add_sequence_feature( SF, "chr1", feature  )
   feature = CytosineCount( SequencePos(105,205),true,45,90)
   add_sequence_feature( SF, "chr1", feature  )
   println(SF)
end
