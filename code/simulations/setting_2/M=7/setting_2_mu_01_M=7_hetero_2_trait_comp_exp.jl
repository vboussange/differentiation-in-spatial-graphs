#= 
This script works with EvoId v5.1.0

- graphs with 7 nodes
- trait competition, with d(x,y) = exp(-0.5(x-y)^2/σ_α^2)
- birth with b(x) = 1-p(x-θ)^2 where p = 1
- 2 different values for habitat: [-1/2, 1/2].
- 3 combinations of assortativity are selected when possible.

Conditions for diversification are given by Doebeli 2016: Adaptive diversification, 

σ_α < 1 / sqrt(2p)

This is why we explore 
σ_α = 0.5 < 1 / sqrt(2p) amd σ_α = 1 > 1 / sqrt(2p)
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
isdir(name_sim) ?  nothing :  mkdir(name_sim)
include("../../../graphs_utils/src/graphs_utils.jl")

######################
## global param def ##
######################

Tf = Float32
nodes = 7
dim_neutr = 300
μ = 1f-1
ms = Float32[1e-2, 5e-1]
σ_αs = Float32[0.5, 1]
adaptivespace = RealSpace{1,Tf}()
neutralspace = RealSpace{dim_neutr,Tf}()
const K1 = 150f0;
NMax = 2000
# tend = Tf(1.5)
tend = 2000.
hetero = [-Tf(0.5), Tf(0.5)]
D2 = Tf(5e-2)
D = [nothing, D2, fill(D2,dim_neutr)]
# We collect betas, betau that we average over time horizon [tend-50, tend]
t_saving_cb = collect(range(tend-99., tend, length=200))

function cb(w)
    Dict(
        "alphas" => get_alpha_div(w,2),
        "betas" => get_beta_div(w,2),
        # "gammas" => mean(var(w,trait=2)),
        "alphau" => get_alpha_div(w,3),
        "betau" => get_beta_div(w,3),
        # "gammau" => mean(var(w,trait=3)),
        "N" => length(w))
end

p_default = Dict{String,Any}();@pack! p_default = NMax,D

######################
## explo params def ##
######################
sleep(rand()) # to make sure that the file is not loaded at the same time
graphs = load("../../../graphs_utils/M=$nodes/graph_prop_M=$nodes.jld2", "graphs_df")[:,"graph"]

pars = []
totgraph = 0
for g in graphs
    soptim_df = DataFrame(soptim = [],rθ = [])
    for i in 1:1000
        _q = nodes ÷ length(hetero); _r = nodes % length(hetero)
        soptim = Dict(1:(nodes-_r) .=> shuffle(repeat(hetero,outer= _q)))
        merge!(soptim,Dict((nodes-_r+1):nodes .=> shuffle(hetero[1:_r])))
        rθ = soptim_correlation(g,soptim)
        if !(rθ in soptim_df.rθ)
            push!(soptim_df,(soptim,rθ))
        end
    end
    _s = size(soptim_df,1)
    global totgraph += _s
    # sorting will bias the sampling
    sort!(soptim_df,:rθ)

    if _s > 2
        soptim_df = soptim_df[[1 ,_s, round(Int,(_s+1)/2)],:]
    end
    for r in eachrow(soptim_df), _m in ms, σ_α in σ_αs
        p = copy(p_default)
        soptim  = r.soptim; rθ = r.rθ;
        # array is slightly faster than dictionary, so better using it for birth rate
        # we still use soptim to compute assortativity
        soptim_arr = [soptim[i] for i in 1:nodes]
        mu = [Tf(_m),Tf(μ),fill(Tf(μ),dim_neutr)]
        geospace = GraphSpace(g)
        myspace = (geospace,adaptivespace,neutralspace)
        b(X,t) = 1f0 - (X[2][] - soptim_arr[Int(X[1][])])^2
        d(X, Y, t) = (X[1][] ≈ Y[1][]) ? exp(- 0.5f0 * (X[2][] - Y[2][])^2 / σ_α^2) / K1 : 0f0
        p = copy(p_default)
        @pack! p = b, d, myspace, soptim, mu, rθ, σ_α
        push!(pars,p)
    end
end
println("Total number of simulations ", length(pars) )
println("We only selected ", length(pars) / length(ms) / length(σ_αs), " graphs")

# @show totgraph
# @show length(pars)
# distrib_rθ = [soptim_correlation(p["myspace"][1].g, p["soptim"])  for p in pars]
# using Plots; histogram(distrib_rθ)

Random.seed!(seed)
println("############################")
println("Starting simulation of file")
println(@__FILE__)
println("Number of threads: $(Threads.nthreads())")
println("Number of simulations: $(length(pars))")
println(now())
println("############################")


function sim(p)
    @unpack b, d, myspace, soptim, rθ, D, mu, σ_α, NMax = p
    myagents = [Agent(myspace, [[xpos], [D[2] * randn(Tf)],
    D[3] .* randn(Tf, dim_neutr)]) for xpos in rand(1:nodes,nodes*Int(K1))]
    w0 = World(myagents, myspace, D, mu, NMax, 0.)
    s = run!(w0,Gillepsie(), tend, b, d,cb=cb, t_saving_cb = copy(t_saving_cb));
    success=true
    return (mu[1][1], 
            mu[2][1], 
            D2, 
            σ_α, 
            myspace[1].g, 
            soptim, 
            rθ, 
            get_tend(s), 
            mean(s["N"][2:end]), 
            mean(s["alphas"][2:end]), 
            mean(s["betas"][2:end]), 
            # mean(s["gammas"][2:end]), 
            mean(s["alphau"][2:end]), 
            mean(s["betau"][2:end]), 
            # mean(s["gammau"][2:end]), 
            success)
end


df = DataFrame("m" => zeros(length(pars)), 
                "μ" => zeros(length(pars)), 
                "D" => zeros(length(pars)), 
                "σ_α" => zeros(length(pars)), 
                "graph" => Vector{SimpleGraph}(undef,length(pars)),
                "soptim" => Vector{Any}(undef,length(pars)),
                "rθ" => zeros(length(pars)),
                "tend" => zeros(length(pars)), 
                "N" => zeros(length(pars)),                
                "alphas" => zeros(length(pars)),
                "betas" => zeros(length(pars)),
                # "gammas" => zeros(length(pars)),
                "alphau" => zeros(length(pars)),
                "betau" => zeros(length(pars)),
                # "gammau" => zeros(length(pars)),
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
