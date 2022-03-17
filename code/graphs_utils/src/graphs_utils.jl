#= 
    This script contains utils functions to calculate
    graph properties.
    We rely on the Python library networkx for some metrics.
    For a first use, do the following:
    ```julia
    using Pkg; Pkg.add("Conda")
    using Conda; Conda.add("networkx")
    # optionally
    Pkg.rm("Conda")
    ```
=#

using LightGraphs
using Statistics,StatsBase
using PyCall; nx = pyimport("networkx")
using LinearAlgebra

## LightGraphs.jl extension
function to_nx(g::AbstractGraph)
    return nx.Graph(Matrix(adjacency_matrix(g)))
end

function to_nx(g::AbstractGraph,soptim::Dict)
    gnx = to_nx(g)
    for a in soptim
        gnx.add_node(a[1]-1,soptim=a[2]+2) #for some reason values need to be positive
    end
    return gnx
end
"""
    to_julia(gnx::PyObject, no_self_loops = true)
Transform a netoworkx graph to a Graphs.SimpleGraph.
If `no_self_loops = true`, we reomve self loops
"""
function to_julia(gnx::PyObject, no_self_loops = true)
    A = nx.to_numpy_matrix(gnx)
    if no_self_loops
        A[diagind(A)] .= 0.
    end
    SimpleGraph(A)
end

function rw_laplacian_matrix(g::AbstractGraph)
    L = laplacian_matrix(g,dir=:out) ./ outdegree(g)
    replace!(L,NaN=>0)
end

## metric networks

function characteristic_length(graph::AbstractGraph)
    !is_connected(graph) && return Inf
    sum(floyd_warshall_shortest_paths(graph).dists) / (nv(graph) * (nv(graph) - 1))
end

import EvoId.groupby
"""
    function distance_distribution(graph::AbstractGraph)
Counts the number of shortest path with length i.
Returns a dictionary with keys corresponding to length of shortest path,
and value corresponding to the frequency (count).
Skipping length 0.
"""
function distance_distribution(graph::AbstractGraph)
    !is_connected(graph) && return Inf
    _l = floyd_warshall_shortest_paths(graph).dists
    _l = _l[.!(iszero.(_l[:]))]
    _g = groupby(x -> x,_l )
    return Dict(keys(_g) .=> length.(values(_g)))
end

function temperature(g::AbstractGraph)
    nei = neighbors.(Ref(g),vertices(g))
    temp = zeros(nv(g))
    for i in 1:nv(g)
        temp[i] = sum(1 ./ nei[i])
    end
    return temp
end

function heat_hetero(g::AbstractGraph)
    return var(temperature(g))
end

function density(g::AbstractGraph)
    m = ne(g); n = nv(g)
    return 2*m / (n*(n-1))
end

function mean_neighb_deg(g)
    sum = 0
    for v in vertices(g)
        sum += mean(degree.(Ref(g),neighbors(g,v)))
    end
    return sum / nv(g)
end

function var_neighb_deg(g::AbstractGraph)
    _m = []
    for v in vertices(g)
        push!(_m,mean(degree.(Ref(g),neighbors(g,v))))
    end
    return var(_m)
end

function degree_correlation(g::AbstractGraph)
    r1 = r2 = r3 = 0
    M = ne(g)
    for e in edges(g)
        ji = degree(g,src(e))
        ki = degree(g,dst(e))
        r1 += ji * ki
        r2 += ji + ki
        r3 += ji^2 + ki^2
    end
    return (r1/M - (0.5 * r2 / M)^2) / (0.5 * r3/M - (0.5 * r2 / M)^2)
end

function degree_correlation_nx(g::AbstractGraph)
    gnx = to_nx(g)
    return nx.degree_assortativity_coefficient(gnx)
end

"""
    rich_club(g,k)

Return the non-normalised [rich-club coefficient](https://en.wikipedia.org/wiki/Rich-club_coefficient) of graph `g`,
with degree cut-off `k`.

```jldoctest
julia> using LightGraphs
julia> g = star_graph(5)
julia> rich_club(g,1)
0.4
```
"""
function rich_club(g::AbstractGraph{T},k::Int) where T
    E = zero(T)
    for e in edges(g)
        if (outdegree(g,src(e)) >= k) && (indegree(g,dst(e)) >= k )
            E +=1
        end
    end
    N = count(degree(g) .>= k)
    if is_directed(g)
        return E / (N*(N-1))
    else
        return 2*E / (N*(N-1))
    end
end

function rich_club_nx(g::AbstractGraph)
    gnx = to_nx(g)
    return values(nx.rich_club_coefficient(gnx,normalized = false))
end

function mean_rich_club(g::AbstractGraph)
    return mean(rich_club.(Ref(g),
                            collect(1:maximum(degree(g))-1)
                            ))
end

function var_rich_club(g::AbstractGraph)
    return var(rich_club.(Ref(g),
                            collect(1:maximum(degree(g))-1)
                            ))
end

"""
    s_metric(g;norm=true)

Return the normalised s-metric of `g`.

The s-metric is defined as the sum of the product of degrees between pair of vertices
for every edge in `g`. [Ref](https://arxiv.org/abs/cond-mat/0501169)
In directed graphs, the paired values are the out-degree of source vertices
and the in-degree of destination vertices.
It is normalised by the maximum s_metric obtained from the family of graph
with similar degree distribution. s_max is computed from an approximation
formula as in https://journals.aps.org/pre/pdf/10.1103/PhysRevE.75.046102
If `norm=false`, no normalisation is performed.

# Examples
```jldoctest
julia> using LightGraphs

julia> s_metric(star_graph(4))
0.6
```
"""

function s_metric(g::AbstractGraph{T};norm=true) where T
    s = zero(T)
    for e in edges(g)
        s += outdegree(g,src(e)) * indegree(g,dst(e))
    end
    if norm
        sm = sum(degree(g).^3)/2
        return s/sm
    else
        return s
    end
end

function s_metric_nx(g::AbstractGraph)
    gnx = to_nx(g)
    nx.s_metric(gnx,normalized = false)
end

function smin(g::AbstractGraph)
    Z = vcat([fill(d,d) for d in sort(degree(g),rev=true)]...)
    return 0.5 * dot(Z,reverse(Z))
end

function s_elasticity(g::AbstractGraph)
    sm = sum(degree(g).^3)/2
    return abs(sm - smin(g))
end


function edge_betweeness_centrality_nx(g::AbstractGraph)
    gnx = to_nx(g)
    return values(nx.edge_betweenness_centrality(gnx))
end

function modules_nx(g::AbstractGraph)
    gnx = to_nx(g)
    return length(nx.algorithms.community.greedy_modularity_communities(gnx))
end

function average_shortest_path_length_nx(g::AbstractGraph)
    gnx = to_nx(g)
    return nx.average_shortest_path_length(gnx)
end

function diameter_nx(g::AbstractGraph)
    gnx = to_nx(g)
    return nx.diameter(gnx)
end

function ncliques_nx(g::AbstractGraph,k)
    gnx = to_nx(g)
    c = collect(nx.enumerate_all_cliques(gnx))
    return count(length.(c) .== k)
end

function graph_energy(g::AbstractGraph)
    a = Matrix(adjacency_matrix(g))
    return sum(abs.(eigen(a).values))
end

function algebraic_conn_nx(g::AbstractGraph)
    gnx = to_nx(g)
    return nx.algebraic_connectivity(gnx)
end


function algebraic_connectivity(g::AbstractGraph)
    a = Matrix(laplacian_matrix(g))
    return eigen(a).values[2]
end

# function autocov_graph(g::AbstractGraph,soptim::Dict)
#     _a = 0
#     for v in vertices(g)
#         for w in neighbors(g,v)
#             _a += soptim[v] * soptim[w]
#         end
#     end
#     return _a / (2 * ne(g))
# end


function autocov_graph(g::AbstractGraph,soptim::Dict)
    _a = 0
    for e in edges(g)
        _a += soptim[src(e)] * soptim[dst(e)]
    end
    return _a
end


function soptim_correlation(g::AbstractGraph,soptim::Dict)
    r1 = r2 = r3 = 0
    M = ne(g)
    for e in edges(g)
        ji = soptim[src(e)]
        ki = soptim[dst(e)]
        r1 += ji * ki
        r2 += ji + ki
        r3 += ji^2 + ki^2
    end
    return (r1/M - (0.5 * r2 / M)^2) / (0.5 * r3/M - (0.5 * r2 / M)^2)
end

function soptim_assortativity_nx(g::AbstractGraph,soptim::Dict)
    gnx = to_nx(g,soptim)
    return nx.numeric_assortativity_coefficient(gnx,"soptim")
end

function size_max_component(g::AbstractGraph,soptim::Dict)
    _group = groupby(x -> soptim[x],collect(1:M))
    _mc = []
    for _c in values(_group)
        push!(_mc,maximum(length.(connected_components(g[_c]))))
    end
    return _mc
end


## those last metrics have not been given enough thoughts
function clustering_coeff(g::AbstractGraph,soptim::Dict)
    _group = groupby(x -> soptim[x],collect(1:M))
    _mc = []
    for _c in values(_group)
        push!(_mc,global_clustering_coefficient(g[_c]))
    end
    return _mc
end

function s_metric(g::AbstractGraph,soptim::Dict)
    _group = groupby(x -> soptim[x],collect(1:M))
    _mc = []
    for _c in values(_group)
        push!(_mc,smax_metric(g[_c]))
    end
    return _mc
end


function ratio_size_components(g::AbstractGraph,soptim::Dict)
    _group = groupby(x -> soptim[x],collect(1:M))
    _mc = []
    for _c in values(_group)
        push!(_mc,mean(length.(connected_components(g[_c]))))
    end
    return maximum(_mc) / minimum(_mc)
end
