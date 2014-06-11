module FeatureTest

    using FactCheck

    include(Pkg.dir("MethylUtils", "src","FeatureServices.jl"))
    path = Pkg.dir("MethylUtils", "testdata","moabs-format-1000lines.txt")

    sf = load_moabs_format_as_sequence_features(path)
    println(sf)

    rf = read_moabs_format_as_sequence_features(path)

    #local_store_path =  Pkg.dir("MethylUtils", "testdata","store", "annotation")
    #cgs = load_moabs(path)

    #SET_LOCAL_STORE_PATH(local_store_path)

    #facts("load_moabs") do
    #    @fact length(cgs) => 1000
    #end

    #transcripts = ensgene_dataframe("mm10","ensGene1000lines.txt.gz" )
    #println(transcripts)

    #facts("ensgene_dataframe") do
    #    @fact size(transcripts,1) => 1000
    #end

    #genes = ensgene_regions("mm10","ensGene1000lines.txt.gz" )



    # println( genes )
    # println("length: ", size(genes.dataframe, 1))
    # We want to test for ensgene_regions that it is making the kind of ensgene_regions
    # we want, i.e. the region covering all transcripts under a gene id

    # set config path back to nothing
    #SET_LOCAL_STORE_PATH("nothing")
end
