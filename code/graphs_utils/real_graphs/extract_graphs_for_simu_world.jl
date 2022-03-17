#=
This script extracts graph from habitat raster
in order to study diversificaiton in realistic graphs

The habitat chosen corresponds to shrubland. 

WE SPAN THE ENTIRE WORLD HERE.
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

filename_hab = "/Users/victorboussange/ETHZ/projects/diversification_in_graphs/Land2Graph/Land2Graph_experiments/data/lvl1_frac_1km_ver004/iucn_habitatclassification_fraction_lvl1__300_Shrubland__ver004.tif"
filename_temp = "/Users/victorboussange/ETHZ/projects/diversification_in_graphs/Land2Graph/Land2Graph_experiments/data/chelsa/CHELSA_bio1_reprojected_nearest.tif"


# loading rasters for (A) and (B)
data_src_temp = rasterio.open(filename_temp, "r")

df = DataFrame("residualmodel" => Any[], "R2" => Float64[], "habitat" => String[], "window_size" => Int32[], "area_threshold" => Int32[])

habitat = "shrubland"
window_size = 10
area_threshold = 100

## local window_size, area_threshold
data_src_hab = rasterio.open(filename_hab, "r")

mytemp = data_src_temp.read(1) /10. .- 273.15 .|> Float64
myhab =  data_src_hab.read(1) .|> Float64

## cropping to adapt yaquan data
# for c in 1:length(myrich)
#     if isnan(myrich[c])
#         myhab[c] = 0 # we do not consider habitats for which we do not have richness values
#     end
# end
cropx = 1:size(myhab,1) - (size(myhab,1) % window_size)
cropy = 1:size(myhab,2) - (size(myhab,2) % window_size)
mytemp = mytemp[cropx, cropy]
myhab = myhab[cropx, cropy]
raster_g, raster_B, raster_sqrtk, raster_assortativity_temp, raster_temp = extract_metrics_1hab(
                                                                                                myhab, 
                                                                                                mytemp, 
                                                                                                habitat, 
                                                                                                window_size, 
                                                                                                area_threshold, 
                                                                                                radius = 2.)
# plotting distribution of # of vertices for connected graphs
figure()
distrib_nodes = length.(vcat(connected_components.(raster_g)...))
plt.hist(distrib_nodes, bins = 100)
plt.yscale("log")
gcf()
# we select connected graphs with 49 vertices
M = 49
count(distrib_nodes .== M)

######################################################################
######### extracting connected graphs with fixe # of vertices ########
######################################################################
# in this case we look at all subgraphs within the window
interesting_graph_unorganised = SimpleGraph{Int16}[]
pos_interesting_graph_unorganised = Any[]
geographic_coord = Any[] #storing geographical coordinate of the window, returning indices top left corner of the large resolution grid
temps = Any[]
for (i,g) in enumerate(raster_g)
    if nv(g) > 0
        _con = connected_components(g)
        _idx = length.(_con) .== M
        if any(_idx)
            # indices of biggest component with M vertices
            idx_big = first(_con[_idx])
            push!(interesting_graph_unorganised, g[idx_big])
            B = raster_B[i]
            coord_graph = CartesianIndices(B)[B][idx_big]
            temp = raster_assortativity_temp[i][B][idx_big]
            # assiging coordinates to the graph nodes
            pos = Dict{Int,Array}()
            for j in 1:length(coord_graph)
                pos[j-1] = [coord_graph[j][2]-1, coord_graph[j][1]-1] # coordinates of the graph in julia array reference frame
            end
            push!(pos_interesting_graph_unorganised, pos)
            push!(temps, temp)
            XY = CartesianIndices(raster_g)[i]
            push!(geographic_coord,XY)
        end
    end
end
length(interesting_graph_unorganised)

# plotting
_idx = 2
println("Graph selected has ", ne(interesting_graph_unorganised[_idx]), " edges.")
gx = to_nx(interesting_graph_unorganised[_idx])
pos = pos_interesting_graph_unorganised[_idx]
geo = geographic_coord[_idx]

fig = figure()
plt.imshow(myhab[(geo[1]*window_size - (window_size-1)):(geo[1]*window_size), (geo[2]*window_size - (window_size-1)):(geo[2]*window_size)].>area_threshold)
nx.draw_networkx(gx,pos,
                    node_color = temps[_idx],
                    with_labels = false,
                    # ax = ax,
                    )
gcf()
close(fig)

# statistics of graphs
cls = characteristic_length.(interesting_graph_unorganised)
sqrtks = sqrtk.(interesting_graph_unorganised)
rθ_temp = assortativity.(interesting_graph_unorganised, temps)

println("Average range of temperatures for graphs selected is ", mean(maximum.(temps) - minimum.(temps)))

# we actually want graphs for which range of temperature or precipitaiton is similar
# otherwise values of β_s will not be comparable (because difference cannot be solely explained by differences in geomtery)
# this is why we select subsets
ranges_temp =  maximum.(temps) - minimum.(temps)
subset_similar_temp = 3.5 .< ranges_temp .< 100

count(subset_similar_temp)

figure()
hist(ranges_temp)
gcf()

fig, ax = subplots(2,2)
ax[1].hist(cls[subset_similar_temp])
ax[1].set_xlabel("Characteristic length")
ax[2].hist(sqrtks[subset_similar_temp])
ax[2].set_xlabel("sqrtk")
ax[4].hist(rθ_temp[subset_similar_temp])
ax[4].set_xlabel("rθ_temp")
fig.tight_layout()
gcf()
close(fig)
# from above analysis, it loooks like choosing temperature as a proxy for environmental condition
# is more suitable as it covers a larger range of r_\theta values

# plotting graph coordinates
fig, ax = subplots(1)
ax.imshow(myhab)
[ax.scatter(geographic_coord[i][2], geographic_coord[i][1], c = "r", marker="+") for i in (1:length(geographic_coord))[subset_similar_temp]]
gcf()

##################################################
#### standardising environmental conditions ######
##################################################
# standardising betwen -0.5 and 0.5
function standardise_unif(x) 
    x = (x.- minimum(x))
    x = x / maximum(x)
    return x .- 0.5
end
temps = standardise_unif.(temps)
using JLD2
JLD2.save("realistic_graph_hengduans.jld2", Dict("geographic_coord" => geographic_coord[subset_similar_temp], 
                                            "graphs" => interesting_graph_unorganised[subset_similar_temp],
                                            "temps" => temps[subset_similar_temp],
                                            "vertices_poss" => pos_interesting_graph_unorganised[subset_similar_temp]))