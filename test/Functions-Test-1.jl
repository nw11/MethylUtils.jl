module FunctionTest1
   using FactCheck
   include(Pkg.dir("MethylUtils", "src","Functions.jl"))
   include(Pkg.dir("MethylUtils", "src","FeatureServices.jl"))


   println("start")
   local_store_path =  Pkg.dir("MethylUtils", "testdata","store", "annotation")
   SET_LOCAL_STORE_PATH(local_store_path)
   moabs_bs_path = Pkg.dir("MethylUtils", "testdata","moabs-format-1000lines.txt")
   organism_id = "mm10"
   ensgene_bodies =  ensgene_regions( organism_id, "ensGene1000lines.txt.gz")
   println(ensgene_bodies)
   println(moabs_bs_path)
   intersect_bed_path = intersect_bed_with_gene_regions( moabs_bs_path, ensgene_bodies)
   println(intersect_bed_path)
   # rm the file for testing
   SET_LOCAL_STORE_PATH("nothing")
end
