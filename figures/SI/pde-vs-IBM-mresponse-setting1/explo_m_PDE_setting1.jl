cd(@__DIR__)
name_sim = split(splitpath(@__FILE__)[end],".")[1]
# name_sim =  "explo_G_2D_discrete_D2_1e-2"
isdir(name_sim) ?  nothing :  mkdir(name_sim)
using DiffEqOperators,LinearAlgebra,LightGraphs,Distributions
using DifferentialEquations,Random
using Plots,Printf
using UnPack,DataFrames,JLD2,Dates
import EvoId:gaussian
using IDEvol
# include(pwd()*"/utils.jl")
simu = true
df_aggreg = false

## Parameters used
mu= 0.1f0
K1 = 150f0 #scaled to 1 for decreasing numerical errors
M = 7;
dS = 0.02f0;
rS = 2f0;
σ_mu = 5f-2;
# c = 1e0 ;
c = 1f0 / K1;
tend = 1000f0
# tend = 1.5
tspan = (0f0,tend)

## rest of the simulation
g = star_graph(M)
S = collect(range(-rS,rS,step=dS)) #grid
N = length(S)
Δ_s = Float32.(Array(σ_mu^2/2 *mu * CenteredDifference{1}(2, 2, dS, N))[:,2:end-1])
X = 1:M
u0 = Float32.(repeat(K1 .* pdf.(Normal(0,σ_mu),S'),outer=(M,1)))

B(x,s) = 1f0

# this function is called to evaluate the convolution of local u and alpha evaluated at S[i]
function int_u(u::T,p::Dict) where T <: AbstractArray
        @unpack N,S,dS = p
        C = 0.5f0 * (u[1] + u[N])
        C += sum(u[2:N-1])
        return C*dS
end

lux = similar(u0);
lus = similar(u0);

# non linear term, i.e. the logistic summand of the IDE
function f(du,u,p,t)
    @unpack N,M,S,dS,X,m,Δ_s,Δ_x= p
    u[u[:,:] .< eps()] .= 0
    B_i = B.(X,S')
    for i in 1:M
        u_i = @view u[i,:]
        C_i = int_u((u_i), p)
        # here there is no b (1-m) because the -m is included in the laplacian
        du[i,:] .= u_i .* (B_i[i,:] .- C_i / K1  )
    end
    # @show size(u)
    # @show size(Δ_x')
    # @show size(B_i)
    # @show size(u .* B_i)
    mul!(lux, Δ_x', (u .* B_i))
    mul!(lus, u .* B_i, Δ_s )
    @. du = du - lux + lus
    return du
end

p_default = Dict{String,Any}()
@pack! p_default = M,S,dS,N,Δ_s,σ_mu,X

    #test
if false
    m = 1.0
    Δ_x = Float64.(m * rw_laplacian_matrix(g))
    @pack! p_default = Δ_x,m
    prob = ODEProblem(f,u0,tspan,p_default)
    @time sol = solve(prob)
    popsize(sol.u[end],p_default)
    Plots.plot(S,sol.u[end][1,:])
end


if true
    ms = collect(range(0f0,1f0,length=30))
    # ms = [0.5]
    allparams = []
    for m in ms
        p = copy(p_default)
        p["m"] = m
        p["Δ_x"] = Float32.(m * rw_laplacian_matrix(g))
        push!(allparams,p)
    end

    df_explo = DataFrame("sol" =>[],
                        "αdiv" => [],
                        "βdiv" => [],
                        "γdiv" => [],
                        "m" => [],
                        # "r_θ" => [],
                        "popsize" => [])


    ################SIMULATION

    for (i,p) in enumerate(allparams)
        prob = ODEProblem(f,u0,tspan,p)
        println("simulation $i / $(length(allparams))")
        try
            sol = solve(prob,alg=Tsit5())
            # we must multiply by K1, because it is scaled to 1f0 
            uend = copy(sol.u[end])
            # @show popsize(sol,p)
            push!(df_explo,(uend,αdiv(uend,p),βdiv(uend,p),γdiv(uend,p),p["m"],popsize(uend,p)))
        catch e
            println("problem with p = $p")
            println(e)
        end
    end
    @save joinpath(name_sim,name_sim*"_$(today()).jld2") df_explo
end
@load joinpath(name_sim,name_sim*"_$(today()).jld2") df_explo

# using PyPlot
# figure()
# plt.plot(S,df_explo.sol[1][1,:])
# gcf()