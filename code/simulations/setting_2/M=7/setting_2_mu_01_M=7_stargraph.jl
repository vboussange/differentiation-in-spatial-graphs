#= 
This script works with EvoId v5.0.0

Here we run simulations on the star graphs for SI figure pde-vs-IBM-mresponse-setting2

=# 

using Random
length(ARGS) > 0 ? seed = (parse(Int,ARGS[1])) : seed = 1
Random.seed!(1) # we set seed to 1 here so that each run has graphs with exactly equivalent rθ 
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
nodes = 7
dim_neutr = 300
μ = 1f-1
ms = collect(range(0f0,1f0,length=30))
adaptivespace = RealSpace{1,Tf}()
neutralspace = RealSpace{dim_neutr,Tf}()
const K1 = 150f0;
@inbounds d(X, Y, t) = (X[1][] ≈ Y[1][]) ? 1f0 / K1 : 0f0
NMax = 3000
# tend = Tf(0.5)
tend = 1000.
hetero = [-Tf(0.5), Tf(0.5)]
D2 = Tf(5e-2)
D = [nothing, D2, fill(D2,dim_neutr)]
# We collect betas, betau that we average over time horizon [tend-99, tend]
t_saving_cb = collect(range(tend-99., tend, length=200))

function cb(w)
    Dict(
        "betas" => get_beta_div(w,2),
        "gammas" => mean(var(w,trait=2)),
        "betau" => get_beta_div(w,3),
        "gammau" => mean(var(w,trait=3)),
        "N" => length(w))
end

p_default = Dict{String,Any}();@pack! p_default = NMax,D

######################
## explo params def ##
######################
sleep(rand()) # to make sure that the file is not loaded at the same time
g = star_graph(nodes)


pars = []
totgraph = 0
for _m in ms
    soptim_arr = 0.5f0 * Float32[-1,1,-1,1,-1,1,-1]
    soptim = Dict(1:nodes .=> soptim_arr)
    rθ = soptim_correlation(g,soptim)
    mu = [Tf(_m),Tf(μ),fill(Tf(μ),dim_neutr)]        
    geospace = GraphSpace(g)
    myspace = (geospace,adaptivespace,neutralspace)
    b(X,t) = 1f0 - (X[2][] - soptim_arr[Int(X[1][])])^2
    p = copy(p_default)
    @pack! p = b, d, myspace, soptim, mu, rθ
    push!(pars,p)
end
# @show totgraph
# @show length(pars)


Random.seed!(seed)
println("############################")
println("Starting simulation of file")
println(@__FILE__)
println("Number of threads: $(Threads.nthreads())")
println("Number of simulations: $(length(pars))")
println(now())
println("############################")


function sim(p)
    @unpack b, d, myspace, soptim, rθ, D, mu, NMax = p
    myagents = [Agent(myspace, [[xpos], [D[2] * randn(Tf)],
                                D[3] .* randn(Tf, dim_neutr)]) for xpos in rand(1:nodes,nodes*Int(K1))]
    w0 = World(myagents, myspace, D, mu, NMax, 0.)
    s = run!(w0,Gillepsie(), tend, b, d,cb=cb, t_saving_cb = copy(t_saving_cb));
    success=true
    return (mu[1][1], mu[2][1], D2, myspace[1].g, soptim, rθ, get_tend(s), mean(s["N"][2:end]), mean(s["betas"][2:end]),mean(s["gammas"][2:end]), mean(s["betau"][2:end]), mean(s["gammau"][2:end]), success)
end


df = DataFrame("m" => zeros(length(pars)), 
                "μ" => zeros(length(pars)), 
                "D" => zeros(length(pars)), 
                "graph" => Vector{SimpleGraph}(undef,length(pars)),
                "soptim" => Vector{Any}(undef,length(pars)),
                "rθ" => zeros(length(pars)),
                "tend" => zeros(length(pars)), 
                "N" => zeros(length(pars)), 
                "betas" => zeros(length(pars)),
                "gammas" => zeros(length(pars)),
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
# |  ETA: 1 days, 22:04:17