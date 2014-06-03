# FeatureServices.jl
#
# Gets features from sources
# whether this be files on disk in certain formats
# or web accessible repositories

using DataFrames
using GZip

include(Pkg.dir("MethylUtils", "src","Features.jl"))
include(Pkg.dir("MethylUtils", "src","Config.jl"))

function load_moabs(path; exclude_unmapped = true )
    #uncompressed_stream = ismatch(r"gz", path) ? GZip.open(path) : open(path)
    uncompressed_stream = open(path)
    cgs = Array(CytosineCount,0)
    idx = 0
    for line in eachline(uncompressed_stream)
       idx +=1
       # Parse line
       # Get the pos and point, i.e. this will result in two Annotation Points, for start and stop.
       (chr, start, stop, ratio,totalC, methC, strand, next, plus, plustotalC, plusmethC, minus, minustotalC,minusmethC ) = [ "" for i in 1:14 ]
       try
           (chr, start, stop, ratio,totalC, methC, strand, next, plus, plustotalC, plusmethC, minus, minustotalC,minusmethC ) = split(line, "\t")
       catch y
             error( "Error: Offending line ($idx) : $line" )
       end
       methC_int  = parseint(methC)
       totalC_int = parseint(totalC)
       cg = CytosineCount(GenomePos(chr, start, stop) , true, methC_int, totalC_int - methC_int )
       push!(cgs,cg)
       if idx % 1000000 == 0
          println("loaded: ",idx)
       end
   end
   return cgs
end

# this creates a file for the start
# and ends of gene regions from the transcripts
# ensgene_dataframe

function ensgene_regions(organism, ensgene_filename="ensGene.txt.gz")
     ensgene_df = ensgene_dataframe(organism,ensgene_filename)
     gene_dict = Dict{Any,Any}()
     nrows = size(ensgene_df,1)
     for row = 1:nrows
         gene = ensgene_df[row, :]
         #println(gene)
         gene_name = gene[:name2][1]
         chr       = gene[:chrom][1]
         start     = gene[:txStart][1]
         stop      = gene[:txEnd][1]
         strand    = gene[:strand][1]
         #println("$gene_name $start $stop")
         if haskey(gene_dict, gene_name) == true
             if gene_dict[gene_name][1] > start
                gene_dict[gene_name][1] = start
             end
             if gene_dict[gene_name][2] < stop
                gene_dict[gene_name][2] = stop
             end
         else
             gene_dict[gene_name] = [start, stop, chr, strand]
         end
     end


     CURRENT_ANNOTATION_CHR_INDEX = current_annotation_chr_index() # from config

     println("The length of the gene_dict is: ", length(gene_dict))
     names  = Array(ASCIIString,0)
     starts = Array(Int32,0)
     stops  = Array(Int32,0)
     chrs   = Array(Int16,0)
     strands = Array(Bool,0)
     for key in keys(gene_dict)
        push!(names,key)
        push!(starts,  gene_dict[key][1])
        push!(stops,   gene_dict[key][2])
        push!(chrs,    CURRENT_ANNOTATION_CHR_INDEX[organism][ gene_dict[key][3] ])
        if(gene_dict[key][4] == "+")
            push!(strands, true)
        else
            push!(strands, false)
        end
     end
     println( "lengths of names: ", length(names) )
     D = DataFrame(name=names,chr=chrs,start=starts,stop=stops,strand=strands)
     A = FeatureTable(["organism" => organism],D)
     return A
end

function ensgene_dataframe(organism, ensgene_filename)
    path = joinpath(LOCAL_STORE_PATH(), organism,ensgene_filename)
    names_str = [ "bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart",
              "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2",
              "cdsStartStat", "cdsEndStat", "exonFrames" ]
    names = map(symbol, names_str)

    #check that ensgene_filename exists
    if !isfile(path)
        error("$path does not exist on disk")
    end
    # try catch this and throw out path?
    ensgene_df = readtable(path,separator='\t', header=false, names=names)
    return ensgene_df
end

#### TO IMPLEMENT ####
# function load_moabs(path,include_intervals::IntervalTree)
#    error(" load_moabs(path,include_intervals) Not implemented")
# end
######################
