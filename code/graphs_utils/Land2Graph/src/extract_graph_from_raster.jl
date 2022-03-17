using LightGraphs

"""
    function extract_adjacency(N::BitMatrix)
    function extract_adjacency(N, radius::AbstractFloat)
Given `N` a matrix of suitable habitat (`true` if suitable habitat at index `i,j`), returns a graph.

If `radius` not specified, habitats are potentially connected to eight neighbours 
    (left, diagonal left top, top, diagonal right top, etc...).
"""
# we keep this function although it is not as general as `extract_adjacency(N, radius::AbstractFloat)`
# because it is supposedly faster (more specialised)
# Indeed, in `extract_adjacency(N, radius::AbstractFloat)`, all potential vertices are checked within the radius, while only half
# of it are required
# this is a downside that could be fixed in the future
function extract_adjacency(N)
    s1 = size(N,1)
    s2 = size(N,2)
    _li = LinearIndices(N)
    A = zeros(Bool,length(N), length(N))
    for i in 1:s1, j in 1:s2
        if N[i,j] # cell is not empty
            if i == s1 && j == s2
                break
            elseif i == s1
                padding = [(0, 1)]
            elseif j == 1
                padding = [(1, 0), (0, 1), (1,1)]
            elseif j == s2
                padding = [(1, -1), (1, 0)]
            else
                padding = [(0,1), (1, -1), (1, 0), (1,1)]
            end
            for _o in padding
                # @show _o
                if N[i + _o[1], j + _o[2]] # neighbour not empty?
                    A[ _li[i, j], _li[ i + _o[1], j + _o[2] ] ] = A[_li[ i + _o[1], j + _o[2] ],  _li[i, j],] = true 
                end
            end
        end
    end
    return A
end

function extract_adjacency(N, radius::AbstractFloat)
    rr = floor(Int,radius)
    s1 = size(N,1)
    s2 = size(N,2)
    _li = LinearIndices(N)
    A = zeros(Bool,length(N), length(N))
    for i in 1:s1, j in 1:s2
        if N[i,j] # cell is not empty
            for dist1 in -rr:rr
                for dist2 in -rr:rr
                    if 0 < dist1^2 + dist2^2 <= radius^2 
                        k = i + dist1
                        l = j + dist2
                        # checking for out of bound
                        if 1 <= k <= s1 && 1 <= l <= s2
                            A[_li[i, j], _li[k, l]] = A[ _li[k, l], _li[i, j]] = true
                        end
                    end
                end
            end
        end
    end
    return A
end

"""
    extract_graph(N; radius = nothing, habitat = 1)

Extract graph from raster `N`, with habitat defined by value `habitat`.

If `radius` not specified, habitats are potentially connected to eight neighbours 
    (left, diagonal left top, top, diagonal right top, etc...).
"""

function extract_graph(N; radius = nothing, habitat = 1)
    B = N .== habitat
    if isnothing(radius)
        A = extract_adjacency(B) #specialised method
    else
        A = extract_adjacency(B, radius)
    end
    g = SimpleGraph{Int16}(A)
    return g[LinearIndices(B)[B]]
end

"""
    extract_graph_1km(N; radius = nothing, area = 0)
Extract graph from fractional raster `N`, given threshold `area`.
Returns graph and coordinate BitMatrix for habitats.

If `radius` not specified, habitats are potentially connected to eight neighbours 
    (left, diagonal left top, top, diagonal right top, etc...).
"""
function extract_graph_1km(N; radius = nothing, area = 0)
    B = N .> area
    if isnothing(radius)
        A = extract_adjacency(B) #specialised method
    else
        A = extract_adjacency(B, radius)
    end
    g = SimpleGraph{Int16}(A)
    return g[LinearIndices(B)[B]], B
end
