#= 
This script works with EvoId v5.0.0

Here we consider the complete graph with M=9
We run simulation for 30 different values of 
m ranging from 0 to 1. 

Used in figure 4.
=# 

using Random, StatsBase
length(ARGS) > 0 ? seed = (parse(Int,ARGS[1])) : seed = 1
Random.seed!(seed)
cd(@__DIR__)
using Dates
using JLD2
using EvoId, LightGraphs, UnPack
using ProgressMeter
name_sim = split(splitpath(@__FILE__)[end],".")[1]
# name_sim =  "explo_G_2D_discrete_D2_1e-2"
isdir(name_sim) ?  nothing :  mkdir(name_sim)
include("../../../graphs_utils/src/graphs_utils.jl")

######################
## global param def ##
######################

Tf = Float32
nodes = 9
dim_neutr = 300
μ = 1f-1
ms = range(0f0, 1f0, length=30)
neutralspace = RealSpace{dim_neutr,Tf}()
const K1 = 150f0;
@inbounds d(X, Y, t) = (X[1][] ≈ Y[1][]) ? 1f0 / K1 : 0f0
b(X,t) = 1f0
NMax = 3000
tend = Tf(1.5)
# tend = 1000.
D2 = Tf(5e-2)
D = [nothing, fill(D2,dim_neutr)]
# We collect betas, betau that we average over time horizon [tend-50, tend]
t_saving_cb = collect(range(tend-99., tend, length=100))

function cb(w)
    Dict("alphau" => get_alpha_div(w,2),
        "betau" => get_beta_div(w,2),
        "gammau" => mean(var(w,trait=2)),
        "N" => length(w))
end

p_default = Dict{String,Any}();@pack! p_default = NMax,D

######################
## explo params def ##
######################
pars = []
g = star_graph(nodes)

for _m in ms
    p = copy(p_default)
    mu = [Tf(_m),fill(Tf(μ),dim_neutr)]        
    geospace = GraphSpace(g)
    myspace = (geospace,neutralspace)
    p = copy(p_default)
    @pack! p = myspace, mu
    push!(pars,p)
end
println("We only selected ", length(pars) / length(ms), " graphs")

# @show totgraph
# @show length(pars)
# distrib_rθ = [soptim_correlation(p["myspace"][1].g, p["soptim"])  for p in pars]
# using Plots; histogram(distrib_rθ)

println("############################")
println("Starting simulation of file")
println(@__FILE__)
println("Number of threads: $(Threads.nthreads())")
println("Number of simulations: $(length(pars))")
println(now())
println("############################")


function sim(p)
    @unpack myspace, mu, NMax = p
    myagents = [Agent(myspace, [[xpos], D[2] .* randn(Tf, dim_neutr)]) for xpos in rand(1:nodes, nodes * Int(K1))]
    w0 = World(myagents, myspace, D, mu, NMax, 0.)
    s = run!(w0,Gillepsie(), tend, b, d,cb=cb, t_saving_cb = copy(t_saving_cb));
    success=true
    return (mu[1][1], mu[2][1], D2, myspace[1].g, get_tend(s), mean(s["N"][2:end]), mean(s["alphau"][2:end]), mean(s["betau"][2:end]), mean(s["gammau"][2:end]), success)
end


df = DataFrame("m" => zeros(length(pars)), 
                "μ" => zeros(length(pars)), 
                "D" => zeros(length(pars)), 
                "graph" => Vector{SimpleGraph}(undef,length(pars)),
                "tend" => zeros(length(pars)), 
                "N" => zeros(length(pars)), 
                "alphau" => zeros(length(pars)), 
                "betau" => zeros(length(pars)),
                "gammau" => zeros(length(pars)), 
                "success" => fill(false,length(pars)))


# Threads.@sync 
progr = Progress(length(pars), showspeed = true, barlen = 10)
@Threads.threads for k in 1:length(pars)
    try
        df[k,:] = sim(pars[k]);
    catch e
        println("problem with p = $(pars[k])")
        println(e)
    end
    next!(progr)
end
@save joinpath(name_sim,name_sim*"_$(today())_seed_$(seed).jld2") df
chmod(name_sim,0o777,recursive=true)