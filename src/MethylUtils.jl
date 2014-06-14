module MethylUtils

# package code goes here
#export EM

export load_bis_seq_file

#Functions
export intersect_files
export run_mlml
export cytosine_methylation_by_gene

#Features
export feature_table_to_bed
export SequenceFeatures
export FeatureRegions
export FeatureTable
export calculate_methylation_level_and_serialize

#FeatureServices
export make_ensgene_gene_feature_table
export load_moabs_format_as_sequence_features
export read_moabs_format_as_sequence_features

include(Pkg.dir("MethylUtils", "src","Functions.jl"))

end
# module
