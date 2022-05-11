cd(@__DIR__)
name_sim = split(splitpath(@__FILE__)[end],".")[1]
# name_sim =  "explo_G_2D_discrete_D2_1e-2"
using DiffEqOperators,LinearAlgebra,LightGraphs,Distributions
using DifferentialEquations,Random
using Plots,Printf
using UnPack,DataFrames,JLD2,Dates,CSV
import EvoId:gaussian
using IDEvol
# include(pwd()*"/utils.jl")
simu = true
df_aggreg = false

if simu
    ## Parameters used
    mu= 0.1
    K1 = 150
    M = 2;
    dS = 0.02;
    rS = 1.0;
    σ_mu = 5e-2;
    c = 1e0 / K1;
    tend = 500
    # tend = 1.5
    tspan = (0.0,tend)

    ## rest of the simulation

    S = collect(range(-rS,rS,step=dS)) #grid
    N = length(S)
    Δ_s = Array( σ_mu^2/2 *mu * CenteredDifference{1}(2, 2, dS, N))[:,2:end-1]
    X = 1:M
    soptim = [-0.5,0.5]
    u0 = vcat([K1 .* pdf.(Normal(so,σ_mu),S') for so in soptim]...)

    B(x,s,soptim) = 1. - (s - soptim[x])^2

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
            return du
    end

    m = 0
    r_θ = 0
    p_default = Dict{String,Any}()
    @pack! p_default = M,S,dS,N,Δ_s,σ_mu,X,soptim,m,r_θ

    #test
    if false
        prob = ODEProblem(f,u0,tspan,p_default)
        @time sol = solve(prob)
        Plots.plot(S,sol.u[end][1,:])
    end

    ms = [0.1,0.5]
    rtheta = collect(range(-1f0,1f0,length=30))
    allparams = []
    for m in ms, rt in rtheta
        p = copy(p_default)
        p["m"] = m
        p["r_θ"] = rt
        push!(allparams,p)
    end

    df_explo = DataFrame("uend" =>[],
                        "α" => [],
                        "β" => [],
                        "γ" =>[],
                        "m" => [],
                        "rθ" => [],
                        "npop" => [])


    ################SIMULATION

    for (i,p) in enumerate(allparams)
        prob = ODEProblem(f,u0,tspan,p)
        println("simulation $i / $(length(allparams))")
        try
            sol = solve(prob)
            uend = copy(sol.u[end])
            # @show popsize(sol,p)
            push!(df_explo,(uend,αdiv(uend,p),βdiv(uend,p),γdiv(uend,p),p["m"],p["r_θ"],popsize(uend,p)))
        catch e
            println("problem with p = $p")
            println(e)
        end
    end
    CSV.write("pde_meanfield_data_m_01-05.csv", df_explo[:,[:β, :m, :rθ, :npop]])
    @save "pde_meanfield_data_m_01-05.jld2" df_explo
end