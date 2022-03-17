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

M = 49

if calculate
    graphs_dict = load("graphs_M49.jld2", "graphs")
    graphs_df = DataFrame("graph" => [], "class" => [])
    for k in keys(graphs_dict)
        for g in graphs_dict[k]
            push!(graphs_df,(g,k))
        end
    end

    ## checking self loops ##
    for k in keys(graphs_dict)
        for g in graphs_dict[k]
            for n in vertices(g)
                if has_edge(g,n,n) 
                    println(k)
                    break
                end
            end
        end
    end

    if calculate
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
        graphs_df[!, "degree_var."] = var.(degree.(graphs_df.graph))
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
end

@load "graph_prop_M=$M.jld2" graphs_df
## playing around
using PyPlot
clf()
dfg = groupby(graphs_df,"class",sort=true)
for (i,df) in enumerate(dfg[[1,2,3,5,7,8,9,14]])
    sc = scatter(df.sqrtk, df.cl, label = df.class[1],  
            color = plt.cm.rainbow((i-1)/9),
            s=10.
            )
end
legend()
gcf()

# using GraphPlot, Fontconfig, Cairo
# dfg[5].class[1]
# gplot(dfg[13].graph[40])

for (i,df) in enumerate(dfg)
    println(i, df.class[1])
end