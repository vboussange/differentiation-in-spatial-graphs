cd(@__DIR__)
name_sim = split(splitpath(@__FILE__)[end],".")[1]
# name_sim =  "explo_G_2D_discrete_D2_1e-2"
using DiffEqOperators,LinearAlgebra,LightGraphs,Distributions
using DifferentialEquations,Random
using Plots,Printf
using UnPack,DataFrames,JLD2,Dates
import EvoId:gaussian
using IDEvol
using PyPlot
include("../../format.jl")

@load "../../../code/simulations/setting_1/M=7/setting_1_mu_01_M=7_stargraph/setting_1_mu_01_M=7_stargraph_2022-01-18_aggreg.jld2" df_aggreg

@load "explo_m_PDE_setting1/explo_m_PDE_setting1_2022-01-19.jld2" df_explo

df_toplot_PDE = df_explo#[df_explo.r_θ .== -1, :]
df_toplot_ABM = df_aggreg
fig,ax = subplots(1,2,figsize = (FIGSIZE_S[2] * 2.1, FIGSIZE_S[1]), sharex = true)

# Q_TS_u diversity
ax[1].errorbar(df_toplot_ABM.m,
        df_toplot_ABM.Q_TS_u_mean,
        yerr = df_toplot_ABM.Q_TS_u_std,
        # c = mapper.to_rgba(log10(_m)),
        label = "IBM simulations",
        fmt = "o",
        ms = 4.)
ax[1].plot(df_toplot_PDE.m,
        df_toplot_PDE.βdiv ./ df_toplot_PDE.γdiv ,
        label = "PDE result",
        # yerr = df_toplot_var[:,:betas],
        # c = mapper.to_rgba(log10(_m)),
        # label = "m = "*(@sprintf "%2.2f" _m),
        # fmt = "o",
        ms = 4.)
ax[1].set_xlabel(L"m"); ax[1].set_ylabel(L"Q_{ST,u}");
gcf()

# beta diversity
ax[2].errorbar(df_toplot_ABM.m,
        df_toplot_ABM.N_mean,
        yerr = df_toplot_ABM.N_std,
        # c = mapper.to_rgba(log10(_m)),
        # label = "IBM simulations",
        fmt = "o",
        ms = 4.)
ax[2].plot(df_toplot_PDE.m,
        df_toplot_PDE.popsize,
        # label = "PDE result",
        # yerr = df_toplot_var[:,:betas],
        # c = mapper.to_rgba(log10(_m)),
        # label = "m = "*(@sprintf "%2.2f" _m),
        # fmt = "o",
        ms = 4.)
ax[2].set_xlabel(L"m"); ax[2].set_ylabel(L"N");
fig.legend(loc="upper center", bbox_to_anchor=(0.5, 1.1),
          ncol=3,
          fancybox=true,
          # shadow=true
          )
# [_ax.set_xscale("log") for _ax in ax]
fig.tight_layout()
fig.savefig("pde-vs-IBM-mresponse-setting1.pdf",
        dpi=1200,
        bbox_inches = "tight",)
gcf()
