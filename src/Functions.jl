using Compose, Gadfly

include(Pkg.dir("MethylUtils", "src","Config.jl"))
include(Pkg.dir("MethylUtils", "src","FeatureServices.jl"))
include(Pkg.dir("MethylUtils", "src","Features.jl"))

function intersect_bed_with_gene_regions( moabs_bs_path, ensgene_bodies, method="bedtools" ; overwrite=false )
      genebodies_bed_path = joinpath(LOCAL_STORE_PATH(),ensgene_bodies.description["organism"], "genebodies.bed" )
      moabs_basename = basename(moabs_bs_path)
      intersected_bed_path = joinpath(LOCAL_STORE_PATH(),ensgene_bodies.description["organism"], "$moabs_basename.bed")

      if isfile(intersected_bed_path) && ( overwrite == false )
          println("Intersect file exists, overwrite is false, skipping intersection ...")
          return intersected_bed_path
      end

      count = feature_table_to_bed(ensgene_bodies, genebodies_bed_path )
      cmd = `sort -k1,1 -k2,2n $genebodies_bed_path` |> `bedtools intersect -sorted -wb -a - -b $moabs_bs_path` |> `cut --complement -f1-6` |> "$intersected_bed_path"
      output = readall(cmd)
      return intersected_bed_path
end

function serialize_cytosine_methylation(cg_methylation,org,id)
     path = joinpath(LOCAL_STORE_PATH(),org,"$id.jld")
     if(isfile(path))
         error("file already exists at this path: $path")
     end
     io = open(path,"w")
     serialize(io,cg_methylation)
end

function deserialize_cytosine_methylation( org,id )
    path = joinpath(LOCAL_STORE_PATH(),org,"$id.jld")
    cgs_meth = open( path ) do file
            deserialize(file)
    end
    return cgs_meth
end

function isserialized(org,id)
  path = joinpath(LOCAL_STORE_PATH(),org,"$id.jld")
  if isfile(path)
     return true
  else
     return false
  end
end

function cytosine_methylation_by_gene( cytosine_methylation, gene_regions::FeatureRegions )
	gene_meth = Array(CytosineMethylation,0)
	gene_dict = Dict()
  count = 0
	for cytosine in cytosine_methylation
        if cytosine.pos.chr > 100
            continue # temporary hack sort out the high chr codes
        end
        gene_iter = getoverlaps(gene_regions, cytosine.pos.chr, cytosine.pos.start, cytosine.pos.stop )
        #if (length( collect(gene_iter) ) > 0 )
        #	push!(gene_meth, cytosine)
        #end
        count +=1
        for (idx,gene) in gene_iter
        	   if ( haskey(gene_dict, gene[:name][1] ) )
                  push!( gene_dict[ gene[:name][1] ], cytosine)
             else
            	    A = Array(CytosineMethylation,0)
            	    push!(A,cytosine)
                  gene_dict[ gene[:name][1] ] = A
            end
        end
        if( count % 10000 == 0)
           println("cytosine_methylation_by_gene process $count ")
        end
    end
    return gene_dict
end

function tile_genes(cytosine_methylation_by_gene, gene_feature_table,  num_tiles =100, min_num_gene_cgs =50)

    gene_dataframe = gene_feature_table.dataframe
    nrows = size(gene_dataframe,1)
    tiled_genes = Array( Array{Float64,1}, num_tiles)
    gene_names = Array(String,0)
    skipped_genes = 0
    for i = 1:nrows
    	  gene_region = gene_dataframe[ i ,:]
        gene_name = gene_region[:name][1]

        if haskey(cytosine_methylation_by_gene, gene_name)
           gene_cgs = cytosine_methylation_by_gene[gene_name]
        else
        	continue
        end
        if length(gene_cgs) < min_num_gene_cgs
        	continue
        end
        gene_start = gene_region[:start][1]
        gene_stop  = gene_region[:stop][1]
        gene_orientation  = gene_region[:strand][1]

        ## MAKE TILES
        println("Make tiles for $gene_name")
        gene_tiles = gene_meth_tiles(gene_cgs,gene_start,gene_stop,gene_orientation, num_tiles)

        ## FILTER NANs
        num_nans = count( x -> isequal( x, NaN32 ), gene_tiles)
        #println("TILE: ", join( gene_tiles, "-") )
        #println("Number of nans $num_nans")
        if num_nans > 10
          skipped_genes += 1
          println("skipping $gene_name because it has $num_nans nans")
          continue
        end

        ## CALC MISSING VALUES
        gene_tiles = float64( interpolateNAN( gene_tiles ) )

        for  tile_num = 1:num_tiles
              if isdefined(tiled_genes, tile_num)
                  push!(tiled_genes[tile_num],  gene_tiles[tile_num] )
              else
                  tiled_genes[tile_num] = [ gene_tiles[tile_num] ]
              end
        end
        push!(gene_names, gene_name)
    end
    println("Length of df: ", length(tiled_genes[1]) )
    D = hcat( [ tiled_genes[i] for i =1:size(tiled_genes,1) ]... )
    println("skipped genes $skipped_genes")
    return (gene_names,D)
end

function gene_meth_tiles(gene_cgs::Vector{CytosineMethylation}, gene_start::Int32, gene_stop::Int32, gene_orientation::Bool, num_tiles)
    gene_length = gene_stop - gene_start
    bins = zeros(Float32, num_tiles)
    bins_count = zeros(Int32,num_tiles)
    #println("num of cgs in gene: ", length(gene_cgs))
    #println("gene region: $gene_stop, $gene_start, $gene_length")
    cg_skipped = 0
    for k=1:length(gene_cgs)

        bp_from_gene_start = gene_cgs[k].pos.start - gene_start
        bin = 0
        if bp_from_gene_start == -1 || bp_from_gene_start == 0
            # if on border of the start pos of a gene then you are in bin 1
            println("INFO: CG on border of gene start: gene_length, $gene_length, gene_start $gene_start, gene_stop $gene_stop, num_tiles $num_tiles, cg start ",   gene_cgs[k].pos.start )
            bin = 1
        elseif bp_from_gene_start < -1 || bp_from_gene_start > gene_length
            println("skipped cg: gene_length, $gene_length, gene_start $gene_start, gene_stop $gene_stop, num_tiles $num_tiles, cg start ",   gene_cgs[k].pos.start)
            cg_skipped +=1
            continue
        else
             gene_percentile = bp_from_gene_start / gene_length
             bin = iceil(  gene_percentile * num_tiles )
        end

       # if bin < 1
       #     error( "Gene tile failed (bin < 1 ($bin) ): bin assignment invalid: gene_length: $gene_length, gene_start $gene_start, gene_stop $gene_stop, num_tiles $num_tiles, cg start ",   gene_cgs[k].pos.start )
       # end

       # if bin > num_tiles
       #     error( "Gene tile failed (bin > num_tiles ($bin) ): bin assignment invalid: gene_length: $gene_length, gene_start $gene_start, gene_stop $gene_stop, num_tiles $num_tiles, cg start ",   gene_cgs[k].pos.start )
       # end

        try
           bins[ bin ] += gene_cgs[k].mc
        catch y
           if isa(y, BoundsError)
           	  error("Bounds Error: bin: $bin, k $k gene_start $gene_start cg_start ", gene_cgs[k].pos.start )
           end
        end
        bins_count[ bin ] +=1
    end
      bins_ave = bins ./ bins_count
      if( cg_skipped > 0 )
          println("CGs SKIPPED: $cg_skipped cg's which were outside gene boundaries")
      end
      if gene_orientation
          return bins_ave
      else
          return reverse(bins_ave)
      end
end

function interpolateNAN( a::Vector{Float32} )
    nans = 0
    last_value = NaN32
    new_a = Array(Float32, length(a))
    for k=1:length(a)
    	#println("value in a: ", a[k])
        if isequal(a[k], NaN32)
           #println("found nan")
           if nans == 0 # first nan we have seen in this run of possible consecutive nans
                if k != 1  # if not the start then assign last value - if is start remains nan
                   last_value = a[k-1]
                   #println("last_value: $last_value")
                end
           end

           nans += 1

           if k == length(a)
                #println("found nans at the end")
                for j=1:nans
                    idx = (k) - (nans - j)
                    new_a[idx] = last_value
                end
                nans=0
          end

        else
            if nans > 0 ## so fill up values of NANs
                if isequal(last_value,NaN32)
                	#CASE WITH NANS at start
                    #println("Found NaNs at beginning")
                    # set previous NaNs to the value we are at now
                    # ()
                    for j=1:nans
                    	new_a[j] = a[k]
                    end
                    nans=0
                    last_value = NaN32
                    new_a[k] = a[k]
                else
                    m = (last_value - a[k] )/(nans+1)
                    #println("m: $m")
                    m_sum = 0
                   for j=1:nans
           	           m_sum += m
           	           v = last_value - m_sum
           	           #println("add value $v")
                       idx = (k-1) - (nans - j) # fill in from first nan to k-1
           	           #println("At $idx")
                       new_a[idx] = v  ## INTERPOLATE ##
                   end
                   new_a[k] = a[k]
                   last_value = NaN32
                   nans=0
                end
                     # CASE WITH NANS AT END, i.e. k = length(a) and it's nan
            else
                    new_a[k] = a[k]
            end
        end
    end
    return new_a
end


function plot_matrix( tiled_genes, path )
     light_blue = color("light blue")
     dark_blue  = color("blue")
     sc = Scale.ContinuousColorScale(Scale.lab_gradient( light_blue,
                                                         dark_blue ) )
     num_rows = size(tiled_genes,1)
     draw(PNG( path, 10inch, 10inch), spy( tiled_genes,  sc, Scale.y_continuous(maxvalue=num_rows) ) )
     # or call the R code, because the time this takes can be excessive
end


# order_matrix
function order_matrix( tiled_genes, cluster_assignments)

    # we want each row to be together with its other cluster members
    # AND we want those groups to be sorted by the number that are in each cluster
    # This means changing the number of the cluster (e.g. 1,2,3 and 4) to the amount
    # that are in the cluster (e.g. 13,445,34 and 67) and sorting on that .

    clusters = unique(cluster_assignments)
    println(clusters)
    num_clusters = length(clusters)

    # obtain counts in each cluster label
    h = hist(cluster_assignments, [i for i=0:num_clusters])
    println(h)
    h = h[2]
    println("Got these clusters $h")
    # now make dictionary mapping the current labels (1 to n )
    # to a label that is the size of how many are in the cluster
    cluster_size = {i => h[i] for i = 1:num_clusters}
    println("cluster size: $cluster_size")
    cluster_assignments_size_labels = Int64[]
    println("do for loop btw: cluster_assignments_size_labels")
    for c in cluster_assignments
         c_size = cluster_size[c]
         push!( cluster_assignments_size_labels, c_size )
    end
    # println("pushed on to $cluster_assignments_size_labels")
    sort_order = sortperm(cluster_assignments_size_labels)
    return tiled_genes[ sort_order,:]
end
