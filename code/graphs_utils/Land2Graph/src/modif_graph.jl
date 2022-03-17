"""
modif_graph(g::SimpleGraph, shortest_dist::Int)

Adds edges to the minimum spanning tree `g` according to the threshold `shortest_dist`.
"""
function modif_graph(g::SimpleGraph, shortest_dist::Int)
    for v in vertices(g)
        dij = dijkstra_shortest_paths(g,v)
        for (w,d) in dij.dists[3:end]
            if d < 2
                add_edge!(g,v,)
        end
    end
end