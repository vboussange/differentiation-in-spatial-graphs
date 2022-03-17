cd(@__DIR__)
using DataFrames
include("code/graphs_utils/src/graphs_utils.jl")
using Test

M = 7
file = readlines("allgraphs_M7_final.txt")
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

# defining soptim
hetero = [-1,1]
_l = length(hetero)
soptim0 = Dict()
for i in 1:_l
    merge!(soptim0,Dict(i:_l:M .=> hetero[i]))
end
soptim1 = Dict(vcat(collect(1:4) .=> hetero[1], collect(5:8) .=> hetero[2]))


@testset "Assortativity" for g in graphs_df.graph
        @test isapprox(degree_correlation_nx(g),degree_correlation(g),atol = 1e-10) || isnan(degree_correlation(g))
        @test isapprox(soptim_assortativity_nx(g,soptim1),soptim_correlation(g,soptim1),atol = 1e-10)
end

@testset "Rich club coefficient" begin
    @testset "Small graphs" for _n in 5:10
        @test @inferred rich_club(star_graph(_n),1) ≈ 2 / _n
        @test @inferred rich_club(DiGraph(star_graph(_n)),1) ≈ 2 / _n
    end
    @testset "Directed ($seed)" for seed in [1, 2, 3], (n, ne) in [(14, 18), (10, 22), (7, 16)]
        g = erdos_renyi(n, ne; is_directed=true, seed=seed)
        _r = rich_club(g,1)
        @test @inferred rich_club(g,1) > 0.
    end
    @testset "Undirected ($seed)" for seed in [1, 2, 3], (n, ne) in [(14, 18), (10, 22), (7, 16)]
        g = erdos_renyi(n, ne; is_directed=false, seed=seed)
        _r = rich_club(g,1)
        @test @inferred rich_club(g,1) > 0.
    end
    # @testset "mean rich club vs networkx" for g in graphs_df.graph
    #     @test mean_rich_club(g) ≈  mean(rich_club_nx(g))
    # end
end

using Random, Statistics

@testset "S-metric" begin
    @testset "Directed ($seed)" for seed in [1, 2, 3], (_n, _ne) in [(14, 18), (10, 22), (7, 16)]
        g = erdos_renyi(_n, _ne; is_directed=true, seed=seed)
        sm = s_metric(g,norm=false)
        sm2 = sum([outdegree(g,src(d)) * indegree(g,dst(d)) for d in edges(g)])
        @test @inferred sm ≈ sm2
        sm = s_metric(g,norm=true)
        sm2 /= sum(degree(g).^3)/2
        @test @inferred sm ≈ sm2
    end
    @testset "Undirected ($seed)" for seed in [1, 2, 3], (_n, _ne) in [(14, 18), (10, 22), (7, 16)]
        g = erdos_renyi(_n, _ne; is_directed=false, seed=seed)
        sm = s_metric(g,norm=false)
        sm2 = sum([degree(g,src(d)) * degree(g,dst(d)) for d in edges(g)])
        @test @inferred sm ≈ sm2
        sm = s_metric(g,norm=true)
        sm2 /= sum(degree(g).^3)/2
        @test @inferred sm ≈ sm2
    end
end

@testset "Characteristic length" begin
    for g in graphs_df.graph
        @test average_shortest_path_length_nx(g) ≈ characteristic_length(g)
    end
end

@testset "Algebraic connectivity" begin
    for g in graphs_df.graph
        @test algebraic_connectivity(g) ≈ algebraic_conn_nx(g)
    end
end
