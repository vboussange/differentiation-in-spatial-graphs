cd(@__DIR__)
using LightGraphs
include("../src/graphs_utils.jl")
using DataFrames
using JLD2
using Random
calculate = true

if calculate
    graphs_df = load("realistic_graph_hengduans.jld2", "sub_df_interesting_graphs")
    graphs = graphs_df.graphs    

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
            graphs_df[!, prefix*"_"*String(fun)] = Float64[ eval(st)(eval(fun)(g)) for g in graphs_df.graphs ]
        end
    end
    println("Calculating non standard other metrics")
    graphs_df[!, "density"] = density.(graphs_df.graphs)
    graphs_df[!, "degree_var."] = var.(degree.(graphs_df.graphs))
    graphs_df[!, "degree_correlation"] = degree_correlation.(graphs_df.graphs) # not sure this one works
    graphs_df[!, "cl"] = characteristic_length.(graphs_df.graphs)
    graphs_df[!, "mean_neighb._deg."] = mean_neighb_deg.(graphs_df.graphs)
    graphs_df[!, "var._neighb._deg."] = var_neighb_deg.(graphs_df.graphs)
    graphs_df[!, "smax"] = s_metric.(graphs_df.graphs,norm=true) # not sure this one works
    graphs_df[!, "s_elasticity"] = s_elasticity.(graphs_df.graphs) # not sure this one works
    # degree correlation elasticity
    graphs_df[!, "modules"] = modules_nx.(graphs_df.graphs) # number of modules
    graphs_df[!, "mean_deg._distrib."] = (mean ∘ collect ∘ values ∘ distance_distribution).(graphs_df.graphs) # number of modules
    graphs_df[!, "var._deg._distrib."] = (var ∘ collect ∘ values ∘ distance_distribution).(graphs_df.graphs) # number of modules
    graphs_df[!, "diameter"] = diameter_nx.(graphs_df.graphs)
    # ncliques take a substantial time to be computed
    # graphs_df[!, "loop3"] = ncliques_nx.(graphs_df.graph,3)
    # graphs_df[!, "clique4"] = ncliques_nx.(graphs_df.graph,4)
    # 4stars
    graphs_df[!, "graph_energy"] = graph_energy.(graphs_df.graphs)
    graphs_df[!, "algebraic_connectivity"] = algebraic_connectivity.(graphs_df.graphs)

    ## extra
    graphs_df[!, "kk"] = [mean(degree(g))^2 / mean(degree(g).^2) for g in graphs_df.graphs ]
    graphs_df[!, "sqrtk"] = [ mean(degree(g).^0.5)^2 / mean(degree(g)) for g in graphs_df.graphs ]
    graphs_df.sqrtk_inv = 1 ./ graphs_df.sqrtk
    graphs_df[!, "heat_hetero"] = [ heat_hetero(g) for g in graphs_df.graphs ]

    ################################
    ## Computing statistics of rθ ##
    ################################
    # graphs_df[!, "rθ_temp"] = zeros(size(graphs_df,1))
    # for r in eachrow(graphs_df)
    #     soptim = Dict(1:nv(r.graphs) .=> r.temps)
    #     r.rθ_temp = soptim_correlation(r.graphs, soptim)
    # end
    @save "graph_prop_realistic_graphs.jld2" graphs_df
else
    @load "graph_prop_realistic_graphs.jld2" graphs_df
    global graphs_df = graphs_df
end

metrics = ["var._local_clustering_coefficient", "mean_betweenness_centrality", "var._betweenness_centrality", "mean_edge_betweeness_centrality_nx", "var._edge_betweeness_centrality_nx", "mean_eigenvector_centrality", "var._eigenvector_centrality", "mean_rich_club_nx", "var._rich_club_nx", "mean_closeness_centrality", "var._closeness_centrality", "density", "degree_var.", "degree_correlation", "cl", "mean_neighb._deg.", "var._neighb._deg.", "smax", "s_elasticity", "modules", "mean_deg._distrib.", "var._deg._distrib.", "diameter", "graph_energy", "algebraic_connectivity", "kk", "sqrtk", "sqrtk_inv", "heat_hetero"]
for m in metrics
    println(m, "counts ", count(isnan.(graphs_df[:,m])), " nan ")
    println(m, "counts ", count(isinf.(graphs_df[:,m])), " inf ")
    println(m, "counts ", count(ismissing.(graphs_df[:,m])), " missing ")
end
## playing around
if false
    using Plots,GraphPlot
    histogram2d(
                graphs_df.var._betweenness_centrality,
                graphs_df.var._rtheta,
                color=:vibrant_grad,
                xaxis = "var._betweenness_centrality",
                yaxis = "var._rtheta")
    ## exploring betweenness_centrality
    gplot(graphs_df.graphs[argmax(graphs_df.var._betweenness_centrality)])
    gplot(graphs_df.graphs[argmin(graphs_df.var._betweenness_centrality)])
end
