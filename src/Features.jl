## Features.jl
## -----------
## Feature Types and basic functionality
## related directly to IO of Feature Types
using DataFrames
using IntervalTrees

type GenomePos
     chr::Int16
     start::Int32
     stop::Int32

     function GenomePos( chr::String,start::String,stop::String, organism="none")
        start_int = parseint(start)
        stop_int  = parseint(stop)
        start     = int32(start_int)
        stop      = int32(stop_int)
        if organism == "none"
            chr       = chr_convert_int16(chr)
        else
            chr  = chr_convert_int16(chr, organism)
        end
        new(chr, start, stop)
    end
    function GenomePos(chr::Int,start::Int,stop::Int)
       new(int16(chr),int32(start),int32(stop))
    end
    function GenomePos(chr::String, start::Int,stop::Int)
       new( chr_convert_int16(chr),int32(start),int32(stop))
    end
end

type CytosineCount
    pos::GenomePos
    strand::Bool
    C_count::Int16
    T_count::Int16
end

type CytosineMethylation
    pos::GenomePos
    strand::Bool
    mc::Float32
    hmc::Float32
    #constructor to ensure float32s
    function CytosineMethylation(pos::GenomePos,strand::Bool, mc,hmc)
       new(pos,strand,float32(mc),float32(hmc))
    end
end

type FeatureTable
  description::Dict{ASCIIString,Any}
  dataframe::DataFrame
end

type FeatureRegions
  description::Dict{ASCIIString,Any}
  regions::Dict{Any,Any}

  function FeatureRegions(FT::FeatureTable)
     sequence_regions = Dict()
     DF = FT.dataframe
     nrows = size(DF,1)
     for row = 1:nrows
         gene = DF[row, :]
         #println(gene)
         gene_name = gene[:name][1]
         start     = gene[:start][1]
         stop      = gene[:stop][1]
         chr       = gene[:chr][1]
         insertintervals!(sequence_regions, chr, start, stop, gene )
     end
     new( ["organism"=>""], sequence_regions )
  end
end

# I think the reason we can't just stick anything into value
# and we have multiple functions for with different value types is because
# IntervalTree didn't let us do that. Check that.
function insertintervals!( sequence_regions_dict, chr, start, stop, value::DataFrame)
    if haskey(sequence_regions_dict, chr)
        interval_tree = sequence_regions_dict[chr]
        interval_tree[(start,stop)] = value
    else
        interval_tree = IntervalTree{Int64,DataFrame}()
        interval_tree[(start,stop)] = value
        sequence_regions_dict[chr] = interval_tree
    end
end

function insertintervals!( sequence_regions_dict, chr, start, stop, value::ASCIIString)
    if haskey(sequence_regions_dict, chr)
        interval_tree = sequence_regions_dict[chr]
        interval_tree[(start,stop)] = value
    else
        interval_tree = IntervalTree{Int64,ASCIIString}()
        interval_tree[(start,stop)] = value
        sequence_regions_dict[chr] = interval_tree
    end
end

function getoverlaps(feature_regions,chr,start,stop)
     # if (haskey(genome_intervals, chr) == true )
          interval_tree = feature_regions.regions[chr]
          iter = intersect(interval_tree, (start,stop))
          return iter
      #else
      #   println("WARN: No key found for $key")
      #   return []
      #end
end


function chr_convert_int16( chr)
    rgx1 = match(r"chr(.)",chr)
    if rgx1 != nothing
      # if it's just numbers
        chr_name = rgx1.captures[1]
        if chr_name == "X"
            return int16(200)
         end
        if chr_name == "Y"
            return int16(201)
        end
        if chr_name == "M"
            return int16(202)
        end
        rgx2 = match(r"^[0-9]+$", string( chr_name ))
        if rgx2 != nothing
           #chr_num = rgx2.captures[1]
            return convert(Int16, parseint(chr_name))
        else
           error("no valid chr number for $chr")
        end
    else
         error("Chromosome field not in valid format e.g. chr1")
    end
end

function organism (FT::FeatureTable)
   return FT.description["organism"]
end

function feature_table_to_bed( features::FeatureTable, path, delim="\t")
         # we have to translate the chromosome
         org = features.organism
         df = features.dataframe
         nrows = size(df,1)
         iostream = open(path,"w")
         count = 0
         for row =  1:nrows
         	  row    = df[row,:]
            name   = row[:name][1]
            start  = row[:start][1]
            stop   = row[:stop][1]
            strand = row[:strand][1]
            strand = strand ? "+" : "-"
         	  chr    = CURRENT_ANNOTATION_CHR_NAMES["mm10"][ row[:chr][1] ]
            line   = join([chr,start,stop,name,0,strand,"\n"],delim)
            print(iostream, line)
            count += 1
         end
         close(iostream)
         return count
end
