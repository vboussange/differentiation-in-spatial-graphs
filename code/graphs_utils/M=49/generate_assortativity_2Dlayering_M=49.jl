#= 
In this script we generate graphs 
with different assortativity values
assortativity is created by disposing the vertices
on a 2D layer.
=#

cd(@__DIR__)
using LightGraphs, GraphPlot, Cairo, Fontconfig
using Colors, ProgressMeter
include("../src/graphs_utils.jl")
@load "graph_prop_M=$M.jld2" graphs_df
dfg = groupby(graphs_df,"class",sort=true)
graphs_df = vcat(dfg[[1,2,3,5,7,8,9,14]]...)

# scaling between (0 and 1)
habfun(x,y,ω) = (sin(2*π*ω*x) + sin(2*π*ω*y) + 2)/4

graphs_df[!,"θ_distrib"] = [Dict[] for _ in 1:size(graphs_df,1)]
graphs_df[!,"rθ"] = [Float64[] for _ in 1:size(graphs_df,1)]

@showprogress for r in eachrow(graphs_df)
    g = r.graph
    soptim_df = DataFrame(soptim = Dict[], rθ = Float64[])
    for i in 1:100, ω in [1e-1,5e-1,1e0, 5e0]
        locs_x, locs_y = spring_layout(g)
        soptim = Dict(1:nv(g) .=> round.(habfun.(locs_x, locs_y, ω)) .- 0.5) # scaling between -1/2 and 1/2
        rθ = soptim_correlation(g,soptim)
        if !(rθ in soptim_df.rθ)
            push!(soptim_df,(soptim,rθ))
        end
    end
    _s = size(soptim_df,1)
    sort!(soptim_df,:rθ)
    ## sampling 2 from the 0, 4th, 2nd and 3rd quartile
    # if _s > 2
    #     soptim_df = soptim_df[sample([1,_s,round(Int,(_s+1)/4),round(Int,3*(_s+1)/4)], 2, replace=false),:]
    # end
    ## taking all 0th, 4th, 2nd and 3rd quartile
    if _s > 4
        soptim_df = soptim_df[[1,_s,round(Int,(_s+1)/4),round(Int,3*(_s+1)/4)],:]
    end
    r.rθ = soptim_df.rθ; r.θ_distrib = soptim_df.soptim
end
    
figure()
r_theta_distrib = vcat(graphs_df.rθ...)
hist(r_theta_distrib)
xlabel(L"r_\theta")
gcf()












###########################################
###### playing around #####################
###########################################

# scaling between (0 and 1) for plotting
habfun(x,y,ω) = (sin(2*π*ω*x) + sin(2*π*ω*y) + 2)/4





M = 49
g = LightGraphs.grid(Int[sqrt(M), sqrt(M)])
locs_x, locs_y = spring_layout(g)

ω = 0.5
θ = round.(habfun.(locs_x, locs_y, ω))
gplot(g, locs_x, locs_y, nodefillc = ColorSchemes.balance[θ])

###############################
## visualising the 2d layout ##
###############################
using PyPlot
ω = 0.1
x = y = range(-1,1,length=100)
im = habfun.(x, y', ω)
imshow(im)
gcf()