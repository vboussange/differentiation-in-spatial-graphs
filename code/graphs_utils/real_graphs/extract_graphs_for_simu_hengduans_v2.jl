#=
This script extracts graph from habitat raster
in order to study diversificaiton in realistic graphs

The habitat chosen corresponds to shrubland

Difference with v1: here we directly extract raster from original tiff at the world level
This is to obtain projection in PlateCarree, which was working in the originnal version of Figure 6
=#

# /!\ this script uses PyCall, that relies on python
# Make sure to install Python, as well as the packages
# "cartopy", "rasterio", "networkx", "matplotlib", "shapely".
# Then uncomment the two following lines, making sure to direct
# ENV["PYTHON"] your your python installation (see https://github.com/JuliaPy/PyCall.jl)
# for more explanations.
# Those two lines can be uncommented after the first usepackage

# ENV["PYTHON"] = "/usr/local/anaconda3/envs/land2graph/bin/python"
# using Pkg; Pkg.build("PyCall")

using PyCall
using Printf
using LightGraphs
using FileIO
using DataFrames
using PyPlot
verbose = false
rasterio = pyimport("rasterio")
warp =  pyimport("rasterio.warp")
cd(@__DIR__)
include("../../graphs_utils/Land2Graph/src/graph_metrics.jl") # graph metrics
include("../../graphs_utils/src/graphs_utils.jl") #utility function
include("../../graphs_utils/Land2Graph/src/extract_graph_from_raster.jl") # function required to extract graphs
include("extract_metrics_1hab.jl") # do the graph extraction for each habitat raster file
function standardise_unif(x; shift=0.) 
    x = (x.- minimum(x))
    x = x / maximum(x)
    return x .- shift
end

filename_temp = "../../graphs_utils/Land2Graph/data/CHELSA_bio1_reprojected_nearest.tif"
filename_hab = "../../graphs_utils/Land2Graph/data/iucn_habitatclassification_fraction_lvl1__100_Forest__ver004.tif"
# loading rasters for (A) and (B)
data_src_temp = rasterio.open(filename_temp, "r")
data_src_hab = rasterio.open(filename_hab, "r")

xmin = 90.78596; xmax = 105.8718; ymin = 22.91373; ymax = 34.47213
# pixel coordinate for longitude / latitudes for top left corner
row, col = data_src_hab.index(xmin, ymax)
row2, col2 = data_src_hab.index(xmax, ymin)
# window_size = 10
area_threshold = 500

mytemp = data_src_temp.read(1, window = rasterio.windows.Window( col, row, col2-col, row2-row)) /10. .- 273.15 .|> Float64
myhab =  data_src_hab.read(1, window = rasterio.windows.Window( col, row, col2-col, row2-row)) .|> Float64

df_rasters = DataFrame("raster_g"=>[], 
                    "raster_B" =>[], 
                    "raster_sqrtk" => [], 
                    "raster_assortativity_temp" => [], 
                    "raster_temp" => [], 
                    "radius" => [], 
                    "window_size" => [])
for radius in [1., 1.5, 2., 2.9, 3., 4., 10], window_size in [7,8,9,10,11,12]
    cropx = 1:size(myhab,1) - (size(myhab,1) % window_size)
    cropy = 1:size(myhab,2) - (size(myhab,2) % window_size)
    mytemp = mytemp[cropx, cropy]
    myhab = myhab[cropx, cropy]
    raster_g, raster_B, raster_sqrtk, raster_assortativity_temp, raster_temp = extract_metrics_1hab(
                                                                                                    myhab, 
                                                                                                    mytemp, 
                                                                                                    window_size, 
                                                                                                    area_threshold, 
                                                                                                    radius = radius)
    push!(df_rasters, (raster_g, raster_B, raster_sqrtk, raster_assortativity_temp, raster_temp, radius, window_size))
end


######################################################################
######### extracting connected graphs with fixe # of vertices ########
######################################################################
M = 49
# in this case we look at all subgraphs within the window
df_interesting_graphs = DataFrame("graphs" => SimpleGraph{Int16}[], 
                                "pos_interesting_graph_unorganised" => Any[],
                                "geographic_coord" => Any[], #storing geographical coordinate of the window, returning indices top left corner of the large resolution grid
                                "temps" => Any[],
                                "radius" => [],
                                "window_size" => []
                                )
@showprogress for r in eachrow(df_rasters)
    rad = r.radius
    for (i,g) in enumerate(r.raster_g)
        if nv(g) > 0
            _con = connected_components(g)
            _idx = length.(_con) .== M
            if any(_idx)
                # indices of biggest component with M vertices
                idx_big = first(_con[_idx])
                B = r.raster_B[i]
                coord_graph = CartesianIndices(B)[B][idx_big]
                temp = r.raster_assortativity_temp[i][B][idx_big]
                # assiging coordinates to the graph nodes
                pos = Dict{Int,Array}()
                for j in 1:length(coord_graph)
                    pos[j-1] = [coord_graph[j][2]-1, coord_graph[j][1]-1] # coordinates of the graph in julia array reference frame
                end
                XY = CartesianIndices(r.raster_g)[i]
                push!(df_interesting_graphs,
                    (g[idx_big],pos,XY,temp,rad, r.window_size))
            end
        end
    end
end
df_interesting_graphs[!,"temp_bin"] = [round.(st) for st in standardise_unif.(df_interesting_graphs.temps)]
println("we have collected ",size(df_interesting_graphs,1), " graphs with M = $M vertices")
println("Average range of temperatures for graphs selected is ", mean(maximum.(df_interesting_graphs.temps) - minimum.(df_interesting_graphs.temps)))

# statistics of graphs
cls = characteristic_length.(df_interesting_graphs.graphs)
sqrtks = sqrtk.(df_interesting_graphs.graphs)
rθ_temp = assortativity.(df_interesting_graphs.graphs, df_interesting_graphs.temps)
rθ_temp_bin = assortativity.(df_interesting_graphs.graphs, df_interesting_graphs.temp_bin)

clf()
figure()
sc = scatter(sqrtks,
            cls, 
            # label = df.class[1],  
            # color = plt.cm.rainbow((i-1)/9),
            s=10.
            )
gcf()

###############################################################
######### Selecting graph with similar range of temp and prec #
###############################################################
# we actually want graphs for which range of temperature or precipitaiton is similar
# otherwise values of β_s will not be comparable (because difference cannot be solely explained by differences in geomtery)
# this is why we select subsets
ranges_temp =  maximum.(df_interesting_graphs.temps) - minimum.(df_interesting_graphs.temps)
subset_similar_temp = 3.5 .< ranges_temp .< 100

println("we have found ", count(subset_similar_temp), " graphs with similar range of temperatures")

###############################################################
######### plotting graphs with similar range of temp and prec #
###############################################################

figure()
hist(ranges_temp)
gcf()

fig, ax = subplots(2,2)
ax[1].hist(cls[subset_similar_temp])
ax[1].set_xlabel("Characteristic length")
ax[2].hist(sqrtks[subset_similar_temp])
ax[2].set_xlabel("sqrtk")
ax[3].hist(rθ_temp_bin[subset_similar_temp])
ax[3].set_xlabel("rθ_temp_bin")
ax[4].hist(rθ_temp[subset_similar_temp])
ax[4].set_xlabel("rθ_temp")
fig.tight_layout()
gcf()
close(fig)
# from above analysis, it loooks like choosing temperature as a proxy for environmental condition
# is more suitable as it covers a larger range of r_\theta values

#############################################
######### stratified sampling ###############
#############################################
# by number of edges
df_interesting_graphs[!,"ne"] = ne.(df_interesting_graphs.graphs)
dfg = groupby(df_interesting_graphs, :ne)
sub_df_interesting_graphs = vcat([dfg[i][1:1,:] for i in 1:2:length(dfg)]...)
# by characteristic length and sqrtk
df_interesting_graphs[!,"cls"] = round.(cls)
df_interesting_graphs[!,"sqrtk"] = round.(sqrtks * 100)
## /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\##
## problem with function `groupby`: ambiguous because defined in EvoId -> close julia and reopen
# make sure to delete temporary file df_interesting_graphs_temp as it is large!
## /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\##
using JLD2
@save "df_interesting_graphs_temp.jld2" df_interesting_graphs
## close and reopen
using JLD2
using PyCall
using Printf
using LightGraphs
using FileIO
using DataFrames
using PyPlot
verbose = false
rasterio = pyimport("rasterio")
warp =  pyimport("rasterio.warp")
cd(@__DIR__)
function standardise_unif(x; shift=0.) 
    x = (x.- minimum(x))
    x = x / maximum(x)
    return x .- shift
end
@load "df_interesting_graphs_temp.jld2" df_interesting_graphs
dfg = groupby(df_interesting_graphs, [:cls,:sqrtk])
using StatsBase
sub_df_interesting_graphs = dfg[1][sample(1:size(dfg[1],1), 2, replace=false),:]
for i in 2:length(dfg)
    sz = size(dfg[i],1)
    if sz > 2
        append!(sub_df_interesting_graphs, dfg[i][sample(1:sz, 2, replace=false),:])
    else
        append!(sub_df_interesting_graphs, dfg[i])
    end
end
# removing columns [:cls,:sqrtk,:ne] as they should not further be used (only for subsampling)
select!(sub_df_interesting_graphs, Not([:cls,:sqrtk,:ne]))

##################################################
#### standardising environmental conditions ######
##################################################
# standardising betwen -0.5 and 0.5
sub_df_interesting_graphs.temps = standardise_unif.(sub_df_interesting_graphs.temps, shift=0.5)
using JLD2
JLD2.save("realistic_graph_hengduans.jld2", Dict("sub_df_interesting_graphs" => sub_df_interesting_graphs))