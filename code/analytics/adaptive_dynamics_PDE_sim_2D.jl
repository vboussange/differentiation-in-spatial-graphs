#= 
    This scripts calculates trait distribution
    to qualitatively illustrate claims about adaptive differentiation in figure
    2
=#

cd(@__DIR__)
name_sim = split(splitpath(@__FILE__)[end],".")[1]
using DiffEqOperators, LinearAlgebra, LightGraphs, Distributions
using DifferentialEquations,Random
using Printf
using UnPack,DataFrames,JLD2,Dates
import EvoId:gaussian
using IDEvol
include("pde_utils.jl")

## Parameters used
mu= 0.1
K1 = 1
M = 2;
dS = 0.02;
rS = 1.0;
σ_mu = 5e-2;
# c = 1e0 ;
c = 1e0 / K1;
tend = 500
# tend = 1.5
tspan = (0.0,tend)

## Overwritten later on
m = 0.5
r_θ = -0.5
## rest of the simulation

S = collect(range(-rS,rS,step=dS)) #grid
N = length(S)
Δ_s = Array( σ_mu^2/2 *mu * CenteredDifference{1}(2, 2, dS, N))[:,2:end-1]
X = 1:M
soptim = [-0.5,0.5]
u0 = vcat([K1 .* pdf.(Normal(so,σ_mu),S') for so in soptim]...)


p_default = Dict{String,Any}()
@pack! p_default = M,S,dS,N,Δ_s,σ_mu,X,soptim,m,r_θ

B(x,s,soptim) = 1e0 - (s - soptim[x])^2

# function to obtain arg of local maximum
function get_types(uend)
    du1 = uend[1,2:end]-uend[1,1:end-1]
    du10 = ((du1[2:end] .* du1[1:end-1]) .< 0)
    c = count(du10)
    println("derivative cancels at ", S[2:end-1][du10])
    if c == 3
        # the middle point where derivative is 0 corresponds to a local minimum
        return S[2:end-1][du10][[1,3]]
    elseif c == 1
        return S[2:end-1][du10]
    else 
        println("Cannot handle uend where derivative cancels $c times")
        return S[2:end-1][du10]
    end
end

# this function is called to evaluate the convolution of local u and alpha evaluated at S[i]
function int_u(u::T,p::Dict) where T <: AbstractArray
        @unpack N,S,dS = p
        C = 0.5 * (u[1] + u[N])
        C += sum(u[2:N-1])
        return C*dS
end

# non linear term, i.e. the logistic summand of the IDE
function f(du,u,p,t)
        @unpack N,M,S,dS,X,soptim,m,r_θ= p
        u[u[:,:] .< eps()] .= 0
        C_1 = int_u((@view u[1,:]), p)
        C_2 = int_u((@view u[2,:]), p)
        B1 = B.(X[1],S,Ref(soptim))
        B2 = B.(X[2],S,Ref(soptim))
        du[1,:] .= u[1,:] .* (B1 * (1-m) .- C_1 / K1  ) +   Δ_s * (u[1,:] .* B1)  + 0.5 * m * ((1 + r_θ) * u[1,:] .* B1 + (1 - r_θ) * u[2,:] .* B2)
        du[2,:] .= u[2,:] .* (B2 * (1-m) .-  C_2 / K1 ) +   Δ_s * (u[2,:] .* B2)  + 0.5 * m * ((1 + r_θ) * u[2,:] .* B2 + (1 - r_θ) * u[1,:] .* B1)
        # du .= du .* .!( u[:,:] .< eps() .* du .< 0.)
        return du
end



prob = ODEProblem(f,u0,tspan,p_default)
@time sol = solve(prob)

uend = copy(sol.u[end])
println("beta diversity is :", βdiv(uend,p_default))

df = DataFrame("uend" => [], "rθ" => Float64[], "m" => Float64[], "N" => Float64[], "β" => Float64[], "s1" => Float64[])
ms = range(0., 1., length=20)
rts = [-0.5, 0.5]
for m in ms, rt in rts
    p = copy(p_default)
    p["m"] = m
    p["r_θ"] = rt
    prob = ODEProblem(f,u0,tspan,p)
    println("m = ", m)
    @time sol = solve(prob,alg=Tsit5())
    uend = copy(sol.u[end])
    # s = get_types(uend)[1]
    # push!(df,(uend, rt, m, popsize(uend,p), βdiv(uend,p), s))

    push!(df,(uend, rt, m, popsize(uend,p), βdiv(uend,p), 0.))
end
dfg = groupby(df,"rθ")
@save "adaptive_dynamics_PDE_data.jld2" dfg S