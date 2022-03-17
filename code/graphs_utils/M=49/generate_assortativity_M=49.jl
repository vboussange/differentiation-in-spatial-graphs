#=
This script aims at generating
habitat distribution with wide 
difference in assortativity values

Compared to 2d_graphs, here we rely on sheer random trials to generate the distribution
=#
cd(@__DIR__)
using LightGraphs
include("../src/graphs_utils.jl")
using DataFrames
using JLD2
using Random
using PyPlot, LaTeXStrings
using ProgressMeter
calculate = true

M = 49
hetero = [-0.5,0.5]
@load "graph_prop_M=$M.jld2" graphs_df
dfg = groupby(graphs_df,"class",sort=true)
graphs_df = vcat(dfg[[1,2,3,5,7,8,9,14]]...)

graphs_df[!,"θ_distrib"] = [Dict[] for _ in 1:size(graphs_df,1)]
graphs_df[!,"rθ"] = [Float64[] for _ in 1:size(graphs_df,1)]
@showprogress for r in eachrow(graphs_df)
    g = r.graph
    soptim_df = DataFrame(soptim = [],rθ = [])
    for i in 1:100000
        _q = M ÷ length(hetero); _r = M % length(hetero)
        soptim = Dict(1:(M-_r) .=> shuffle(repeat(hetero,outer= _q)))
        merge!(soptim,Dict((M-_r+1):M .=> shuffle(hetero[1:_r])))
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

idx = 51
fig = figure()
nx.draw_networkx(to_nx(graphs_df.graph[idx]),
                    node_color = [graphs_df.θ_distrib[idx][2][v] for v in 1:M],
                    with_labels = false,
                    # ax = ax,
                    node_size = 50.
                    )
gcf()