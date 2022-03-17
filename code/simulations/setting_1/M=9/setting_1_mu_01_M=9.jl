#= 
This script works with EvoId v5.0.0

Here we explore graphs with 9 nodes.
We consider setting 1.
We only select 1/4 of the ensemble of graphs with 8 nodes.
We solve for 5 different values of m:  [1e-2, 5e-2, 1e-1, 5e-1, 1e0]
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
ms = Float32[1e-2, 5e-2, 1e-1, 5e-1, 1e0]
neutralspace = RealSpace{dim_neutr,Tf}()
const K1 = 150f0;
@inbounds d(X, Y, t) = (X[1][] ≈ Y[1][]) ? 1f0 / K1 : 0f0
b(X,t) = 1f0
NMax = 3000
# tend = Tf(1.5)
tend = 1000.
D2 = Tf(5e-2)
D = [nothing, fill(D2,dim_neutr)]
# We collect betas, betau that we average over time horizon [tend-99, tend]
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
sleep(rand()) # to make sure that the file is not loaded at the same time
graphs = load("../../../graphs_utils/M=$nodes/graph_prop_M=$nodes.jld2", "graphs_df")[:,"graph"]

pars = []
for g in graphs, _m in ms
    p = copy(p_default)
    # array is slightly faster than dictionary, so better using it for birth rate
    # we still use soptim to compute assortativity
    mu = [Tf(_m),fill(Tf(μ),dim_neutr)]        
    geospace = GraphSpace(g)
    myspace = (geospace,neutralspace)
    p = copy(p_default)
    @pack! p = myspace, mu
    push!(pars,p)
end
println("Total number of simulations ", length(pars) )

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
# |  ETA: 2 days, 4:41:44 (33.74  s/it)