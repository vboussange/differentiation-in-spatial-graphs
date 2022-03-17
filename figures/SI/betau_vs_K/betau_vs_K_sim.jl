#=
This script works with EvoId v5.0.0

Here we simulate transient dynamics to empirically 
determine the time scale of further simulations

Parameters:
- graphs : star graph and line graph
- m : 0.1 and 0.3 
- tend = 1000.
=# 
using Random
cd(@__DIR__)
using Dates
using JLD2
using EvoId,LightGraphs,UnPack
using ProgressMeter

######################
## global param def ##
######################

Tf = Float32
nodes = 7
dim_neutr = 300
μ = 1f-1
m = 1f-1
adaptivespace = RealSpace{1,Tf}()
neutralspace = RealSpace{dim_neutr,Tf}()
NMax = 10000
# tend = Tf(1.5)
tend = 1000.
D2 = Tf(5e-2)
D = [nothing, D2, fill(D2,dim_neutr)]
# We collect betas, betau that we average over time horizon [400, 500]
t_saving_cb = collect(range(tend-99., tend, length=200))

function cb(w)
        Dict(
                "alphas" => get_alpha_div(w,2),
                "betas" => get_beta_div(w,2),
                "gammas" => mean(var(w,trait=2)),
                "alphau" => get_alpha_div(w,3),
                "betau" => get_beta_div(w,3),
                "gammau" => mean(var(w,trait=3)),
                "N" => length(w))
end

p_default = Dict{String,Any}();@pack! p_default = NMax,D


pars = []
totgraph = 0
for K1 in Float32[50,150,300,500]
        for g in [complete_graph(nodes)] #
                soptim = 0.5f0 * Float32[-1,1,-1,1,-1,1,-1]
                mu = [Tf(m),Tf(μ),fill(Tf(μ),dim_neutr)]        
                geospace = GraphSpace(g)
                myspace = (geospace,adaptivespace,neutralspace)
                b(X,t) = 1f0 - (X[2][] - soptim[Int(X[1][])])^2
                @inbounds d(X, Y, t) = (X[1][] ≈ Y[1][]) ? 1f0 / K1 : 0f0
                p = copy(p_default)
                @pack! p = b, d, myspace, mu, soptim, K1
                push!(pars,p)
        end
end
nruns = 5
function sim(p)
        @unpack b, d, myspace, D, mu, NMax, soptim, K1 = p
        ss = []
        for i in 1:nruns
                myagents = [Agent(myspace, [[xpos], [Tf(soptim[xpos]) + D[2] * randn(Tf)],
                                                D[3] .* randn(Tf, dim_neutr)]) for xpos in rand(1:nodes,Int(K1))]
                w0 = World(myagents, myspace, D, mu, NMax, 0.)
                s = run!(w0,Gillepsie(), tend, b, d,cb=cb,
                         t_saving_cb = t_saving_cb
                                );
                push!(ss,s)
        end
        success=true
        return (mu[1][1], 
                mu[2][1], 
                D2, 
                K1,
                myspace[1].g, 
                soptim, 
                [get_tend(s) for s in ss], 
                [mean(s["N"][2:end]) for s in ss], 
                [mean(s["alphas"][2:end]) for s in ss], 
                [mean(s["betas"][2:end]) for s in ss], 
                [mean(s["gammas"][2:end]) for s in ss], 
                [mean(s["alphau"][2:end]) for s in ss], 
                [mean(s["betau"][2:end]) for s in ss], 
                [mean(s["gammau"][2:end]) for s in ss], 
                success)
end
    
#########################
## running experiments ##
#########################
progr = Progress(length(pars), showspeed = true, barlen = 10)
df = DataFrame("m" => zeros(length(pars)), 
                "μ" => zeros(length(pars)), 
                "D" => zeros(length(pars)), 
                "K1" => zeros(length(pars)),
                "graph" => Vector{SimpleGraph}(undef,length(pars)),
                "soptim" => Vector{Any}(undef,length(pars)),
                "tend" => [Float64[] for i in 1:length(pars)], 
                "N" => [Float64[] for i in 1:length(pars)], 
                "alphas" => [Float64[] for i in 1:length(pars)],
                "betas" => [Float64[] for i in 1:length(pars)],
                "gammas" => [Float64[] for i in 1:length(pars)],
                "alphau" => [Float64[] for i in 1:length(pars)],
                "betau" => [Float64[] for i in 1:length(pars)],
                "gammau" => [Float64[] for i in 1:length(pars)],
                "success" => fill(false,length(pars)))

# @Threads.threads 
@Threads.threads for k in 1:length(pars)
    try
        df[k,:] = sim(pars[k]);
    catch e
        println("problem with p = $(pars[k])")
        println(e)
    end
    next!(progr)
end

@save "betau_vs_K.jld2" df