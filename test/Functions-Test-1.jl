module FunctionTest1
   using FactCheck
   include(Pkg.dir("MethylUtils", "src","Functions.jl"))
   include(Pkg.dir("MethylUtils", "src","FeatureServices.jl"))
   file1 = Pkg.dir("MethylUtils", "testdata", "moabs-format-1000lines.txt")
   file2 = Pkg.dir("MethylUtils", "testdata", "ensGene1000lines.txt.gz")

   ##### READ TEST DATA #####

   ## cytosine counts ##
   cytosine_counts         =   load_moabs_format_as_sequence_features(file1)
   println("Cytosine Counts")
   println(cytosine_counts)

   ## annotation file ##
   ensgene_gene_regions_ft =   make_ensgene_gene_feature_table(file2, "mm10")
   ensgene_gene_regions    =   FeatureRegions(ensgene_gene_regions_ft)

   ####====####





   ###### METHYLATION ESTIMATION FUNCTONS #################

   cytosine_methylation    =   run_mc_hmc_calc( cytosine_counts )

   println("Cytosines Methylation")
   println(cytosine_methylation)
   println(ensgene_gene_regions_ft)

   ####====####





   ###### TILING RELATED  ##############
   by_gene = cytosine_methylation_by_gene( cytosine_methylation, ensgene_gene_regions )
   println(by_gene)

   ####====####




   ###### UTILS ############

   # TEST INTERSECT BED
   println("Intersect $file1 and $file2")
   mm10_ensgene_path = Pkg.dir("MethylUtils", "testdata", "ensGene1000lines.bed")
   feature_table_to_bed(ensgene_gene_regions_ft, mm10_ensgene_path )
   path = intersect_files(file1, mm10_ensgene_path, overwrite=true)

   println("intersected files path: $path")

   ####====####

end



  # TEST MLML
  # test that we have order cytosines

  # include(Pkg.dir("MethylUtils", "src","FeatureServices.jl"))
  # println("start")
  # local_store_path =  Pkg.dir("MethylUtils", "testdata","store", "annotation")
  # SET_LOCAL_STORE_PATH(local_store_path)
  # moabs_bs_path = Pkg.dir("MethylUtils", "testdata","moabs-format-1000lines.txt")
  # organism_id = "mm10"
  # ensgene_bodies =  ensgene_regions( organism_id, "ensGene1000lines.txt.gz")
  # println(ensgene_bodies)
  # println(moabs_bs_path)
  # intersect_bed_path = intersect_bed_with_gene_regions( moabs_bs_path, ensgene_bodies)
  # println(intersect_bed_path)
  # # rm the file for testing
  # SET_LOCAL_STORE_PATH("nothing")

  #rm files for testing next time
  # Pkg.dir("MethylUtils", "testdata", "ensGene1000lines.bed")
  # Pkg.dir("MethylUtils", "testdata", "moabs-format-1000lines.txt.intersect.bed")
