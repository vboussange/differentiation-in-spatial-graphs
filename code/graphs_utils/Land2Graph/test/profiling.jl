# profiling code
cd(@__DIR__)
using Pkg; Pkg.activate("../.")
# using PyPlot
using ArchGDAL; const AG = ArchGDAL
using Statistics
using ProgressMeter
using Glob, JLD2
using Dates
using BenchmarkTools
include("../src/assortativity.jl")
include("../src/extract_graph_from_raster.jl")

window_size = 100
# preallocating
data = ones(Int16,  window_size, window_size) #x, y , window size
temp = randn(window_size,window_size)

function sqrtk(g) 
    _deg = degree(g)
    return mean(_deg.^0.5)^2 / mean(_deg)
end

@btime ncells = count(data .> 0)
@time if ncells > 0
    g, B = extract_graph_1km(data)
    sqrtk(g) 
    assortativity(g,temp[B])
end
