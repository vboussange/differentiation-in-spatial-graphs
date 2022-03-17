# Contains graphs metrics

function sqrtk(g) 
    _deg = degree(g)
    return mean(_deg.^0.5)^2 / mean(_deg)
end

using SimpleTraits
"""
    assortativity(g)
    assortativity(g, attributes)

Return the [assortativity coefficient](https://en.wikipedia.org/wiki/Assortativity)
of graph `g`, defined as the Pearson correlation of excess degree between
the end vertices of all the edges of the graph.

If `attributes` is provided, Pearson correlation is calculated
from atrribute values associated to each vertex and stored in
vector `attributes`.

The excess degree is equal to the degree of linked vertices minus one,
i.e. discounting the edge that links the pair.
In directed graphs, the paired values are the out-degree of source vertices
and the in-degree of destination vertices.

# Arguments
- `attributes` (optional) is a vector that associates to each vertex `i` a scalar value
`attributes[i]`.

# Examples
```jldoctest
julia> using LightGraphs

julia> assortativity(star_graph(4))
-1.0

julia> attributes = [-1., -1., 1., 1.]

julia> assortativity(star_graph(4),attributes)
-0.5
```
"""
function assortativity(g::AbstractGraph{T}) where T
    P = promote_type(Int64, T) # at least Int64 to reduce risk of overflow
    nue  = ne(g)
    sjk = sj = sk = sjs = sks = zero(P)
    for d in edges(g)
        j = P(outdegree(g, src(d)) - 1)
        k = P(indegree(g, dst(d)) - 1)
        sjk += j*k
        sj  += j
        sk  += k
        sjs += j^2
        sks += k^2
    end
    return assortativity_coefficient(g, sjk, sj, sk, sjs, sks, nue)
end

function assortativity(g::AbstractGraph{T}, attributes::Vector{R}) where {T<:Integer,R<:Real}
    P = promote_type(Float64, R) # at least Float64 to reduce risk of overflow
    nue  = ne(g)
    sjk = sj = sk = sjs = sks = zero(P)
    for d in edges(g)
        j = attributes[src(d)]
        k = attributes[dst(d)]
        sjk += j*k
        sj  += j
        sk  += k
        sjs += j^2
        sks += k^2
    end
    return assortativity_coefficient(g, sjk, sj, sk, sjs, sks, nue)
end

#=
assortativity coefficient for directed graphs:
see equation (21) in M. E. J. Newman: Mixing patterns in networks, Phys. Rev. E 67, 026126 (2003),
http://arxiv.org/abs/cond-mat/0209450
=#
@traitfn function assortativity_coefficient(g::::IsDirected, sjk, sj, sk, sjs, sks, nue)
    return (sjk - sj*sk/nue) / sqrt((sjs - sj^2/nue)*(sks - sk^2/nue))
end

#=
assortativity coefficient for undirected graphs:
see equation (4) in M. E. J. Newman: Assortative mixing in networks, Phys. Rev. Lett. 89, 208701 (2002),
http://arxiv.org/abs/cond-mat/0205405/
=#
@traitfn function assortativity_coefficient(g::::(!IsDirected), sjk, sj, sk, sjs, sks, nue)
    return (sjk/nue - ((sj + sk)/(2*nue))^2) / ((sjs + sks)/(2*nue) - ((sj + sk)/(2*nue))^2)
end