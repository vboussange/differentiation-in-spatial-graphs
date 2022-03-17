#=
This script simulates transients for setting (1),
used in the conceptual figure.
=#

using Random
cd(@__DIR__)
using Dates
using JLD2
using EvoId,LightGraphs,UnPack
using LaTeXStrings

######################
## global param def ##
######################

g = SimpleGraph(Edge.([(1,2),(2,4),(2,5),(2,3),(3,6),(3,7)]))
nodes = 7
geospace = GraphSpace(g)
Tf = Float32
nodes = 7
dim_neutr = 300
μ = 1f-1
adaptivespace = RealSpace{1,Tf}()
neutralspace = RealSpace{dim_neutr,Tf}()
hetero = [-Tf(0.5), Tf(0.5)]
myspace = (geospace,neutralspace)
const K1 = 150f0;
soptim = Dict(vcat(collect(1:2:7) .=> hetero[1],collect(2:2:7) .=> hetero[2] ))
soptim_arr = [soptim[v] for v in 1:nodes]
@inbounds d(X, Y, t) = (X[1][] ≈ Y[1][]) ? 1f0 / K1 : 0f0
b(X,t) = 1f0
NMax = 3000
# tend = Tf(1.5)
tend = 1000.
hetero = [-Tf(0.5), Tf(0.5)]
D2 = Tf(5e-2)
D = [nothing, fill(D2,dim_neutr)]
mu = [Tf(0.05),fill(Tf(μ),dim_neutr)]

myagents = [Agent(myspace, [[xpos],
            D[2] .* randn(Tf, dim_neutr)]) for xpos in rand(1:nodes,nodes*Int(K1))]
w0 = World(myagents, myspace, D, mu, NMax, 0.)

# running simulation
@time s_setting1 = run!(w0, Gillepsie(), tend, b, d);

@save "s_setting1.jld2" s_setting1
