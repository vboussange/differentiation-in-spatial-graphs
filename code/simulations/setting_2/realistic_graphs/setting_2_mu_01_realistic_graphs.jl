#= 
This script works with EvoId v5.0.0

Here we explore realistic graphs.
We consider habitats that we artificially scale from -1/2 to 1/2 
in order to respect similar setting as above
We run simulations for _m: [5e-2, 5e-1]

Graphs have been obtained from a region centered on the hengduan region
with varying radius and screening window size

2021-12-30: tend = 500, missed m = 1e0
2022-01-XX: tend = 1000, included m = 1e0
2022-02-XX: tend = 3000, included m = 1e0

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
dim_neutr = 300
μ = 1f-1
ms = Float32[1e-2, 5e-1,]
adaptivespace = RealSpace{1,Tf}()
neutralspace = RealSpace{dim_neutr,Tf}()
const K1 = 150f0;
@inbounds d(X, Y, t) = (X[1][] ≈ Y[1][]) ? 1f0 / K1 : 0f0
NMax = 100000
# tend = Tf(1.5)
tend = 2000.
D2 = Tf(5e-2)
D = [nothing, D2, fill(D2,dim_neutr)]
# We collect betas, betau that we average over time horizon [tend-50, tend]
t_saving_cb = collect(range(tend-99., tend, length=100))

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
graphs_df = load("../../../graphs_utils/real_graphs/realistic_graph_hengduans.jld2", "sub_df_interesting_graphs")
graphs = graphs_df.graphs
soptim_arr =  graphs_df.temps
pars = []
for (i,g) in enumerate(graphs)
    for _m in ms
        p = copy(p_default)
        soptim  = soptim_arr[i] .|> Tf;
        rθ = soptim_correlation(g,Dict(1:length(soptim) .=> soptim))
        mu = [Tf(_m),Tf(μ),fill(Tf(μ),dim_neutr)]        
        geospace = GraphSpace(g)
        myspace = (geospace,adaptivespace,neutralspace)
        b(X,t) = 1f0 - (X[2][] - soptim[Int(X[1][])])^2
        p = copy(p_default)
        @pack! p = b, d, myspace, soptim, mu, rθ
        push!(pars,p)
    end
end
println("total number of graphs to simulate is ", length(graphs))
# println("We only selected ", length(pars) / length(ms))

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
    @unpack b, d, myspace, soptim, rθ, D, mu, NMax = p
    nodes = nv(myspace[1].g)
    myagents = [Agent(myspace, [[xpos], [D[2] * randn(Tf)],
                                D[3] .* randn(Tf, dim_neutr)]) for xpos in rand(1:nodes,nodes*Int(K1))]
    w0 = World(myagents, myspace, D, mu, NMax, 0.)
    s = run!(w0,Gillepsie(), tend, b, d,cb=cb, t_saving_cb = copy(t_saving_cb));
    success=true
    return (mu[1][1], mu[2][1], D2, myspace[1].g, soptim, rθ, get_tend(s), mean(s["N"][2:end]), mean(s["alphas"][2:end]), mean(s["betas"][2:end]), mean(s["alphau"][2:end]), mean(s["betau"][2:end]), success)
end


df = DataFrame("m" => zeros(length(pars)), 
                "μ" => zeros(length(pars)), 
                "D" => zeros(length(pars)), 
                "graph" => Vector{SimpleGraph}(undef,length(pars)),
                "soptim" => Vector{Any}(undef,length(pars)),
                "rθ" => zeros(length(pars)),
                "tend" => zeros(length(pars)), 
                "N" => zeros(length(pars)), 
                "alphas" => zeros(length(pars)),
                "betas" => zeros(length(pars)),
                "alphau" => zeros(length(pars)),
                "betau" => zeros(length(pars)),
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
