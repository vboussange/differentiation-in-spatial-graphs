#=
This script loads database of graphs with M vertices
and compute topology metrics

=#

cd(@__DIR__)
using LightGraphs
include("../src/graphs_utils.jl")
using DataFrames
using JLD2
using Random
calculate = true

M = 7

if calculate
    file = readlines("allgraphs_M$(M)_final.txt")
    graphs_df = DataFrame("graph" => [])
    for l in file
        g = SimpleGraph(M)
        arr = eval(Meta.parse(l))
        for e in arr[1:end]
            b = add_edge!(g,e[1],e[2])
            if !b
                error("problem with adding edge")
            end
        end
        push!(graphs_df.graph,g)
    end
    funs = Symbol.([
                    local_clustering_coefficient,
                    # degree_centrality,
                    betweenness_centrality,
                    edge_betweeness_centrality_nx,
                    eigenvector_centrality,
                    rich_club_nx,
                    closeness_centrality,
                    ])
    sts = [:mean,:var]
    for fun in funs
        for st in sts
            println("Calculating ", String(st)," ", String(fun))
            prefix = st == :var ? "var." : "mean"
            graphs_df[!, prefix*"_"*String(fun)] = Float64[ eval(st)(eval(fun)(g)) for g in graphs_df.graph ]
        end
    end
    println("Calculating non standard other metrics")
    graphs_df[!, "density"] = density.(graphs_df.graph)
    graphs_df[!, "degree_var."] = var.(LightGraphs.degree.(graphs_df.graph))
    graphs_df[!, "degree_correlation"] = degree_correlation.(graphs_df.graph) # not sure this one works
    graphs_df[!, "cl"] = characteristic_length.(graphs_df.graph)
    graphs_df[!, "mean_neighb._deg."] = mean_neighb_deg.(graphs_df.graph)
    graphs_df[!, "var._neighb._deg."] = var_neighb_deg.(graphs_df.graph)
    graphs_df[!, "smax"] = s_metric.(graphs_df.graph,norm=true) # not sure this one works
    graphs_df[!, "s_elasticity"] = s_elasticity.(graphs_df.graph) # not sure this one works
    # degree correlation elasticity
    graphs_df[!, "modules"] = modules_nx.(graphs_df.graph) # number of modules
    graphs_df[!, "mean_deg._distrib."] = (mean ∘ collect ∘ values ∘ distance_distribution).(graphs_df.graph) # number of modules
    graphs_df[!, "var._deg._distrib."] = (var ∘ collect ∘ values ∘ distance_distribution).(graphs_df.graph) # number of modules
    graphs_df[!, "diameter"] = diameter_nx.(graphs_df.graph)
    # ncliques take a substantial time to be computed
    # graphs_df[!, "loop3"] = ncliques_nx.(graphs_df.graph,3)
    # graphs_df[!, "clique4"] = ncliques_nx.(graphs_df.graph,4)
    # 4stars
    graphs_df[!, "graph_energy"] = graph_energy.(graphs_df.graph)
    graphs_df[!, "algebraic_connectivity"] = algebraic_connectivity.(graphs_df.graph)

    ## extra
    graphs_df[!, "kk"] = [mean(degree(g))^2 / mean(degree(g).^2) for g in graphs_df.graph ]
    graphs_df[!, "sqrtk"] = [ mean(degree(g).^0.5)^2 / mean(degree(g)) for g in graphs_df.graph ]
    graphs_df.sqrtk_inv = 1 ./ graphs_df.sqrtk
    graphs_df[!, "heat_hetero"] = [ heat_hetero(g) for g in graphs_df.graph ]

    ################################
    ## Computing statistics of rθ ##
    ################################
    # graphs_df[!, [, "mean_rtheta",
    #             "var._rtheta"]] = zeros(size(graphs_df,1))
    # for (i,g) in enumerate(graphs_df.graph)
    #     asss = get_all_assortativity(g)
    #     graphs_df[i,"var._rtheta"] = var(asss)
    #     graphs_df[i,:mean_rtheta] = mean(asss)
    # end
    @save "graph_prop_M=$M.jld2" graphs_df
else
    @load "graph_prop_M=$M.jld2" graphs_df
    global graphs_df = graphs_df
end

metrics = names(graphs_df)[3:end]
for m in metrics
    println(m, "counts ", count(isnan.(graphs_df[:,m])), " nan ")
    println(m, "counts ", count(isinf.(graphs_df[:,m])), " inf ")
    println(m, "counts ", count(ismissing.(graphs_df[:,m])), " missing ")
end

## playing around
using PyPlot
clf()
sc = scatter(graphs_df.sqrtk, graphs_df.cl, 
            s=10.
            )
gcf()