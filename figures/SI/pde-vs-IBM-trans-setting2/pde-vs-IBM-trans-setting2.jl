#=
This script runs simulations of both pde and ibm models
in setting 2, to compare transients.
=#

cd(@__DIR__)
name_sim = split(splitpath(@__FILE__)[end],".")[1]
# name_sim =  "explo_G_2D_discrete_D2_1e-2"
# isdir(name_sim) ?  nothing :  mkdir(name_sim)
using DiffEqOperators,LinearAlgebra,LightGraphs,Distributions
using DifferentialEquations,Random
using Plots,Printf
using UnPack,DataFrames,JLD2,Dates
import EvoId:gaussian
using PyPlot
using KernelDensity
using EvoId
include("../../format.jl")
include("../../../code/analytics/pde_utils.jl")

cm_eth = ColorMap([c for c in eth_grad_std.colors])

function _scale(x)
    (x .- mean(x)) ./ std(x)
end

# l : for plotting local or global trait distribution
function plot_abm(ax, s, l = false)
    # @show i
    Y = []; X = []
    for i in 1:length(s[:])
        if l # only obtaining individuals on vertex 2 
            y = []
            for a in s[i]
                if a[1][] == 2
                    push!(y, a[2][])
                end
            end
        else
            y = [a[2][] for a in s[i]]
        end
        x = fill(s.tspan[i],length(y))
        Y = vcat(Y,y); X = vcat(X,x)
    end
    yy = _scale(Y)
    xy = Float64.(hcat(X,yy))
    k = InterpKDE(kde(xy))
    c = KernelDensity.pdf.(Ref(k),xy[:,1],xy[:,2])
    mesh = ax.scatter(X,#discarding some value
                    Y,
                    c = c,
                    s = 10,
                    cmap = cm_eth)
end

tend = 500f0
K1 = 150f0
M = 7;
g = star_graph(M)
D2 = 5f-2

###################################
########### param PDE #############
###################################
rS = 1f0;
dS = 0.02f0;
S = collect(range(-rS,rS,step=dS)) #grid
tspan = (0.0,tend)
N = length(S)
Δ_s = Float32.(Array( D2^2/2 * 1f-1 * CenteredDifference{1}(2, 2, dS, N))[:,2:end-1])
X = 1:M
soptim = 0.5f0 * Float32[-1,1,-1,1,-1,1,-1]
u0 = Float32.(repeat(K1 .* pdf.(Normal(0,D2),S'),outer=(M,1)))
B(x,s) = max(0f0,1f0 - (soptim[x] - s)^2)

###################################
########### param IBM #############
###################################
adaptivespace = RealSpace{1,Float32}()
myspace = (GraphSpace(g),adaptivespace)
@inbounds d(X, Y, t) = (X[1][] ≈ Y[1][]) ? 1f0 / K1 : 0f0
b(X,t) = 1f0 - (X[2][] - soptim[Int(X[1][])])^2
NMax = 3000
# tend = 1.5
D = [nothing, D2]
t_saving_cb = collect(range(0.f0,tend,length=300))

if false
    function sim(m)
        mu = [m,1f-1]
        Δ_x = Float32.(m * rw_laplacian_matrix(g))
        #############################
        ########### PDE #############
        #############################
        # this function is called to evaluate the convolution of local u and alpha evaluated at S[i]
        function int_u(u::T,p::Dict) where T <: AbstractArray
                @unpack N,S,dS = p
                C = 0.5 * (u[1] + u[N])
                C += sum(u[2:N-1])
                return C*dS
        end
        lux = similar(u0);
        lus = similar(u0);
        # non linear term, i.e. the logistic summand of the IDE
        function f(du,u,p,t)
            @unpack N,M,S,dS,X,m,Δ_s= p
            u[u[:,:] .< eps()] .= 0
            B_i = B.(X,S')
            for i in 1:M
                u_i = @view u[i,:]
                C_i = int_u((u_i), p)
                # here there is no b (1-m) because the -m is included in the laplacian
                du[i,:] .= u_i .* (B_i[i,:] .- C_i / K1  )
            end
            mul!(lux, Δ_x', (u .* B_i))
            mul!(lus, u .* B_i, Δ_s )
            @. du = du - lux + lus
            return du
        end
        #
        p_default_pde = Dict{String,Any}()
        @pack! p_default_pde = M,S,dS,N,Δ_s,X,soptim,m

        #test
        if false
            prob = ODEProblem(f,u0,tspan,p_default_pde)
            @time sol = solve(prob)
            Plots.plot(S,sol.u[end][1,:])
        end
        prob = ODEProblem(f,u0,tspan,p_default_pde)
        sol = solve(prob)

        #############################
        ########### IBM #############
        #############################

        function cb(w)
                Dict(
                    "alphau" => get_alpha_div(w,2),
                    "betau" => get_beta_div(w,2),
                    "gammau" => mean(var(w,trait=2)),
                    "N" => size(w))
        end
        myagents = [Agent(myspace, [[xpos], [D[2] * randn(Float32)]]) for xpos in rand(1:M,M*Int(K1))]
        w0 = World(myagents, myspace, D, mu, NMax, 0.)
        @time s = run!(w0,Gillepsie(),tend,b,d,dt_saving = 5f0);
        return s, sol
    end
    # @time s = run!(w0,Gillepsie(),tend,b,d,cb=cb,t_saving_cb = copy(t_saving_cb),dt_saving = 1.0);


    # migration regimes plotted
    m1 = 0.1f0
    m2 = 0.8f0

    s1, sol1 = sim(m1)
    s2, sol2 = sim(m2)
    @save "simulations.jld2" s1 sol1 s2 sol2

end
@load "simulations.jld2" s1 sol1 s2 sol2

sar = [s1, s2]
solar = [sol1, sol2]
mar = [m1, m2]
## Plotting for metapopulation
cla()

fig,ax= plt.subplots(2,2,figsize = (FIGSIZE_L[1],FIGSIZE_L[1]),sharey=true,sharex = true)
for (j,i) in enumerate([0,1])
    s = sar[j]
    sol = solar[j]
    m = mar[j]
    ax[1 + i ].set_ylabel("Local trait distrib.")
    tsspan_pde =  collect(0:5:tend)
    pdesolX = hcat([sum(sol(t),dims=1)' for t in 0:5:tend]...)
    pc = ax[1 + i].pcolormesh(tsspan_pde, S, pdesolX,cmap = cm_eth,
        # extent=[1,tend,S[1],S[end]]
        vmax = 4 * K1
        # norm = matplotlib.colors.LogNorm()
        )
    # fig.colorbar(pc,ax = ax[1])
    ax[1 + i].set_ylim(-1.0,1.0)
    ax[1 + i].set_xlim(0,tend)
    # ax.hist(get_x(w0,2),bins=50)
    plot_abm(ax[3 + i],s)
    ax[3 + i ].set_ylim(-1.0,1.0)
    ax[3 + i ].text(25,0.8,"m = $m",fontsize=12)
    fig.tight_layout()
    gcf()
end
ax[1].set_title("PDE simulation")
ax[3 ].set_title("IBM simulation")
ax[2 ].set_xlabel("Time");ax[4].set_xlabel("Time")

fig.savefig("pde-vs-IBM-trans-setting2_metapop_mu=$mar.pdf",
    dpi=500,
    bbox_inches = "tight",)
gcf()

## Plotting local only local trait distrib on vertex 2
cla()

fig,ax= plt.subplots(2,2,figsize = (FIGSIZE_L[1],FIGSIZE_L[1]),sharey=true,sharex = true)
for (j,i) in enumerate([0,1])
    s = sar[j]
    sol = solar[j]
    m = mar[j]
    ax[1 + i ].set_ylabel("Local trait density")
    tsspan_pde =  collect(0:5:tend)
    pdesolX = hcat([sol(t)[2,:] for t in 0:5:tend]...)
    pc = ax[1 + i].pcolormesh(tsspan_pde, S, pdesolX,cmap = cm_eth,
        # extent=[1,tend,S[1],S[end]]
        vmax = 4 * K1
        # norm = matplotlib.colors.LogNorm()
        )
    # fig.colorbar(pc,ax = ax[1])
    ax[1 + i].set_ylim(-1.0,1.0)
    ax[1 + i].set_xlim(0,tend)
    # ax.hist(get_x(w0,2),bins=50)
    plot_abm(ax[3 + i],s,true)
    ax[3 + i ].set_ylim(-1.0,1.0)
    ax[3 + i ].text(25,0.8,"m = $m",fontsize=12)
    fig.tight_layout()
    gcf()
end
ax[1].set_title("PDE simulation")
ax[3 ].set_title("IBM simulation")
ax[2 ].set_xlabel("Time");ax[4].set_xlabel("Time")

fig.savefig("pde-vs-IBM-trans-setting2_localpop.pdf",
            dpi=500,
            bbox_inches = "tight",)

gcf()
