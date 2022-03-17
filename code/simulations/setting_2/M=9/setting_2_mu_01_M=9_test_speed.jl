# This script works with EvoId v5.0.0

using Random
cd(@__DIR__)
using Dates
using JLD2
using EvoId, LightGraphs, UnPack
using Revise
# name_sim =  "explo_G_2D_discrete_D2_1e-2"

######################
## global param def ##
######################

Tf = Float32
nodes = 9
dim_neutr = 300
adaptivespace = RealSpace{1,Tf}()
neutralspace = RealSpace{dim_neutr,Tf}()
K1 = Tf(150);
@inbounds d(X, Y, t) = (X[1][] ≈ Y[1][]) ? 1f0 / Tf(K1) : 0f0
NMax = 2000
tend = Tf(10)
# tend = 500.
hetero = [-Tf(1),Tf(1)]
D2 = Tf(5e-2)
D = [nothing, D2, fill(D2,dim_neutr)]
mu = [Tf(0.1),Tf(1e-1),fill(Tf(1e-1),dim_neutr)]
t_saving_cb = collect(range(0.,tend,length=300))

geospace = GraphSpace(star_graph(nodes))
myspace = (geospace,adaptivespace,neutralspace)
soptim = [rand([Tf(1), Tf(-1)]) for i in 1:nodes]
b(X,t) = max(one(Tf) - (X[2][] - soptim[Int(X[1][])])^2, zero(Tf))

myagents = [Agent(myspace, [[xpos], [Tf(soptim[xpos]) + D[2] * randn(Tf)],
                                D[3] .* randn(Tf, dim_neutr)]) for xpos in rand(1:nodes,Int(K1))]
w0 = World(myagents, myspace, D, mu, NMax, 0.)
function cb(w)
    Dict(
        # "alphau" => get_alpha_div(w,3),
        "betau" => get_beta_div(w,3),
        # "gammau" => mean(var(w,trait=3)),
        # "alphas" => get_alpha_div(w,2),
        "betas" => get_beta_div(w,2),
        # "gammas" => mean(var(w,trait=2)),
        "N" => length(w))
end
@time run!(w0, Gillepsie(), tend, b, d, cb=cb, t_saving_cb = copy(t_saving_cb));