#= 
Generation of theoretical graphs with large number of vertices

We try to generate a minimum of 500 graphs
=#
cd(@__DIR__)
include("../../graphs_utils/src/graphs_utils.jl")
using LightGraphs
using RCall
R"""
library(OCNet)
"""
"""
    create_river_graph(x,y)
Wrapper around the package OCNet to create river like networks
From optimal channel networks
Note that we convert the directed graph to an undirected graph
We use the "Flow Direction" aggregation (see https://cran.r-project.org/web/packages/OCNet/vignettes/OCNet.html#application-metapopulation-model)
"""
function create_river_graph(x,y)
    W = rcopy(R"""
        OCN <- create_OCN($x,$y)
        as.matrix(OCN$FD$W)
        """)
    return SimpleDiGraph(W) |> SimpleGraph
end


M = 49

# we create a dictionary which maps graph class to an array of graphs
graphs = Dict{String,Vector{SimpleGraph}}()

graphs["star_graph"] = [star_graph(M)]
graphs["complete_graph"] = [complete_graph(M)]
graphs["lattice"] = [LightGraphs.grid(Int64[sqrt(M),sqrt(M)])]
graphs["line"] = [LightGraphs.grid(Int64[M,1])]
# graphs["binary_tree"] = [binary_tree(M)]

# small world
temp_graphs = SimpleGraph[]
for k in 2:10, p in range(0.1,0.9,length = 10)
    g  = watts_strogatz(M, k, p)
    while !is_connected(g)
        g  = watts_strogatz(M, k, p)
    end
    push!(temp_graphs,g)
end
graphs["watts_strogatz"] = temp_graphs
graphs["cycle"] = [cycle_graph(M)]
# Random geometric graphs 
# looks like r = 0.3 is the smallest radius possible to have a connected graph
temp_graphs = SimpleGraph[]
for r in range(0.3,0.9,length = 200), d in [2]
    g = to_julia(nx.random_geometric_graph(M,r,dim=d)) #r is the radius, from in the unit cube
    while !is_connected(g)
        g = to_julia(nx.random_geometric_graph(M,r,dim=d))#r is the radius, from in the unit cube
    end
    push!(temp_graphs,g)
end
graphs["random_geometric_graph"] = temp_graphs

# geometric waxman graph
temp_graphs = SimpleGraph[]
for alpha in range(0.3,0.9,length = 10), beta in range(0.3,0.9,length = 10)
    g = to_julia(nx.waxman_graph(M, alpha = alpha, beta = beta))
    while !is_connected(g)
        g = to_julia(nx.waxman_graph(M, alpha = alpha, beta = beta))
    end
    push!(temp_graphs,g)
end
graphs["waxman_graph"] = temp_graphs

# navigable_small_world_graph
temp_graphs = SimpleGraph[]
for p in 1:4, q in 2:9
    g = to_julia(nx.to_undirected(nx.navigable_small_world_graph(sqrt(M) |> Int, p=p, q=q)))
    while !is_connected(g)
        g = to_julia(nx.to_undirected(nx.navigable_small_world_graph(sqrt(M) |> Int, p=p, q=q)))
    end
    push!(temp_graphs,g)
end
graphs["navigable_small_world_graph"] = temp_graphs

temp_graphs = SimpleGraph[]
for ne in range(M,(M*(M-1)/2 - 1),length=50)
    g = SimpleGraph(M, floor(ne) |> Int)
    while !is_connected(g)
        g = SimpleGraph(M, floor(ne) |> Int)
    end
    push!(temp_graphs,g)
end
graphs["erdos_renyi"] = temp_graphs

temp_graphs = SimpleGraph[]
for _ in 1:50
    g = create_river_graph(sqrt(M) |> Int, sqrt(M) |> Int)
    while !is_connected(g)
        g = create_river_graph(sqrt(M) |> Int, sqrt(M) |> Int)
    end
    push!(temp_graphs,g)
end
graphs["river_graph"] = temp_graphs

temp_graphs = SimpleGraph[]
for n1 in 2:ceil(M/2)
    n2 = M-n1
    g = complete_bipartite_graph(n1|>Int,n2|>Int)
    push!(temp_graphs,g)
end
graphs["bipartite"] = temp_graphs

graphs["random_trees"] = [to_julia(nx.random_tree(M)) for _ in 1:50]

graphs["full_rary_tree"] = [to_julia(nx.full_rary_tree(r,M)) for r in 1:M-2] # M-1 corresponds to the star graph

# adding real graphs calculated from hengduans
graphs["graphs_hengduans"] = load("../real_graphs/realistic_graph_hengduans.jld2", "graphs")

graphs["lollipop_graph"] = [lollipop_graph(M-i,i) for i in 10:(M-10)]

println("We have collected ",length(vcat(values(graphs)...)), " graphs.")

using JLD2
@save "graphs_M49.jld2" graphs