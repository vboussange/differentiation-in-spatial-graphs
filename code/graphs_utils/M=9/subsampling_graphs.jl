using Random, StatsBase, LightGraphs
cd(@__DIR__)

M = 9

file = readlines("allgraphs_M$(M)_final.txt")
# subsampling the total list of graph
graphs = []
for l in file
    g = SimpleGraph(M)
    arr = eval(Meta.parse(l))
    for e in arr[1:end]
        er = add_edge!(g,e[1],e[2])
        if !er
            error("problem with adding edge")
        end
    end
    push!(graphs,g)
end

nb_ne_M_1 = count(ne.(graphs) .== M-1)
subsamp_graph = SimpleGraph[]
for _nv in (M-1):Int(M * (M-1)/2)
    subsamp_temp = graphs[ne.(graphs) .== _nv]
    nb = min(nb_ne_M_1, length(subsamp_temp))
    append!(subsamp_graph, sample(subsamp_temp, nb, replace = false))
end

@save "subsamp_graphs_M=$M.jld2" subsamp_graph
# dict_specs_graphs = Dict{String, SimpleGraph}()
# dict_specs_graphs["line"] = LightGraphs.grid([1,M])
# dict_specs_graphs["flake"] = star_graph(M)
# dict_specs_graphs["star"] = star_graph(M)
# dict_specs_graphs["star"] = star_graph(M)
