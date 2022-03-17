using Test
cd(@__DIR__)
include("../src/extract_graph_from_raster.jl")

N = [0 1 0; 1 0 1;0 1 0] #raster matrix
g = extract_graph(N)
@test nv(g) == 4
@test ne(g) == 4

# gplot(g, nodelabel = 1:nv(g))

N = [1 1 1; 1 0 1;1 1 1] #raster matrix
g = extract_graph(N)
@test nv(g) == 8
@test ne(g) == 12

################################
## general method with radius ##
################################
N = [1 1 1; 1 0 1;1 1 1] #raster matrix
g = extract_graph(N, radius = 1.5)
@test nv(g) == 8
@test ne(g) == 12

N = [1 1 1; 1 0 1;1 1 1] #raster matrix
g = extract_graph(N, radius = 2.)
@test nv(g) == 8
@test ne(g) == 18

# gplot(g, nodelabel = 1:nv(g))
@testset "general vs specialised method" for _ in 1:10
    N = rand([0,1],20,20)
    g1 = extract_graph(N)
    g2 = extract_graph(N, 1.5)
    @test g1 == g2
    @test typeof(g) <: SimpleGraph
end
# gplot(g)

## checking the visually the graphs with coordinates
# using PyPlot
# clf()
# include("../../../graphs_utils/src/graphs_utils.jl") #utility function
# coord_graph = CartesianIndices(N)[N .|> Bool]
# # assiging coordinates to the graph nodes
# pos = Dict{Int,Array}()
# for j in 1:length(coord_graph)
#     pos[j-1] = [coord_graph[j][1]-1, coord_graph[j][2]-1] # coordinates of the graph in julia array reference frame
# end
# nx.draw_networkx(to_nx(g),pos)
# gcf()