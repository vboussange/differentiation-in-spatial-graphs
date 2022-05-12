#=
Plotting std effect on Q_ST_u and Q_ST_s in setting 2 with trait dependent competition
for M = 7 to display in SI

v3 is last version as of 20-01-2022
=#

cd(@__DIR__)
using EvoId,JLD2,FileIO
using DataFrames
using Printf;#pyplot()
using Glob
using LightGraphs
using LaTeXStrings
using KernelDensity
using GLM, Interpolations, Polynomials
using Plots:ColorGradient
using ColorSchemes
using PyPlot
cm_eth = ColorMap([c for c in eth_grad_std.colors]);
include("../../../code/graphs_utils/src/graphs_utils.jl")
# isdir("img") ? nothing : mkdir("img")
## graphs prop

include("../../format.jl")
M = 7
date_sim = "2022-01-18" #t=1000
# date_sim = "2022-02-26" #t=3000
# date_sim = "2022-03-17" #t=2000

@load "../../../code/simulations/setting_2/M=$M/setting_2_mu_01_M=$(M)_hetero_2_trait_comp_exp/setting_2_mu_01_M=$(M)_hetero_2_trait_comp_exp_$(date_sim)_aggreg.jld2" df_aggreg

df_aggreg_g_sigmaa = groupby(df_aggreg,:σ_α,sort=true)

#####################
######## β_u ########
#####################
for _df_aggreg in df_aggreg_g_sigmaa

    _df_aggreg_g = groupby(_df_aggreg,:m,sort=true)
    totsize = size(_df_aggreg_g[1],1)

    fig,axs = subplots(1,2,figsize = (FIGSIZE_S[2] * 2.1, FIGSIZE_S[1]))

    axa,axb = axs
    # calculating GLM for each m and extracting coefficients for Q_ST_u
    sqrtk_coeff = []; sqrtk_coeff_err = []; cl_coeff = []; cl_coeff_err = []; rθ_coeff = []; rθ_coeff_err = []; lms_u = []; r2s = []
    for i in 1:length(_df_aggreg_g)
        df_temp = DataFrame(_df_aggreg_g[i][:,["Q_ST_u_mean","sqrtk", "cl", "rθ"]]); [df_temp[!,n] = _scale(df_temp[:,n]) .|> Float64 for n in names(df_temp)]
        mylm = lm(@formula(Q_ST_u_mean ~ sqrtk + cl + rθ), df_temp)
        push!(lms_u, mylm)
        push!(r2s, r2(mylm))
        sqrtk = coeftable(mylm).cols[1][2]
        push!(sqrtk_coeff,sqrtk);
        push!(sqrtk_coeff_err,[sqrtk - coeftable(mylm).cols[5][2],coeftable(mylm).cols[6][2] - sqrtk]);
        cl = coeftable(mylm).cols[1][3]
        push!(cl_coeff,cl)
        push!(cl_coeff_err,[cl - coeftable(mylm).cols[5][3],coeftable(mylm).cols[6][3] - cl]);
        rθ = coeftable(mylm).cols[1][4]
        push!(rθ_coeff,rθ)
        push!(rθ_coeff_err,[rθ - coeftable(mylm).cols[5][4],coeftable(mylm).cols[6][4] - rθ]);
    end
    ms = unique(df_aggreg.m)
    cols = ["tab:blue", "tab:orange", "tab:red"]
    axb.errorbar(ms,sqrtk_coeff, yerr = hcat(sqrtk_coeff_err...), label = L"h_d", capsize = 2., capthick = 1., c = cols[1])
    axb.errorbar(ms,cl_coeff, yerr = hcat(cl_coeff_err...), label = L"cl", capsize = 2., capthick = 1., c = cols[2])
    axb.errorbar(ms,rθ_coeff, yerr = hcat(rθ_coeff_err...), label = L"r_\theta", capsize = 2., capthick = 1., c = cols[3])
    axb.hlines(0., 1e-2, 1e0, colors="grey", linestyles = "--", label = "")

    rθ = -1. # rθ
    p = 1; θ=0.5
    mstar = 1. ./ (1 .- rθ ) .* 4 * p * θ^2/(1 + 3*p * θ^2) #critical threshold
    @show mstar
    axb.axvline(mstar,label=L"m^*",linestyle="--" )

    # axb.legend()
    axb.set_xscale("log")
    axb.set_ylabel("Standardized effect on "*L"Q_{ST,u}")
    axb.set_xlabel(L"m")
    # showing r2 for m = 0.05
    # axb.text(0.1, 0.9, 
    #         L"m = 0.05: R^2 = "*@sprintf("%1.2f",r2s[2]),
    #         transform=axb.transAxes,
    #         fontsize=12)
    gcf()


    # calculating GLM for each m and extracting coefficients for Q_ST_s
    sqrtk_coeff = []; sqrtk_coeff_err = []; cl_coeff = []; cl_coeff_err = []; rθ_coeff = []; rθ_coeff_err = []; lms_s = []; r2s = []
    for i in 1:length(_df_aggreg_g)
        df_temp = DataFrame(_df_aggreg_g[i][:,["Q_ST_s_mean","sqrtk", "cl", "rθ"]]); [df_temp[!,n] = _scale(df_temp[:,n]) .|> Float64 for n in names(df_temp)]
        mylm = lm(@formula(Q_ST_s_mean ~ sqrtk + cl + rθ), df_temp)
        push!(lms_s, mylm)
        push!(r2s, r2(mylm))
        sqrtk = coeftable(mylm).cols[1][2]
        push!(sqrtk_coeff,sqrtk);
        push!(sqrtk_coeff_err,[sqrtk - coeftable(mylm).cols[5][2],coeftable(mylm).cols[6][2] - sqrtk]);
        cl = coeftable(mylm).cols[1][3]
        push!(cl_coeff,cl)
        push!(cl_coeff_err,[cl - coeftable(mylm).cols[5][3],coeftable(mylm).cols[6][3] - cl]);
        rθ = coeftable(mylm).cols[1][4]
        push!(rθ_coeff,rθ)
        push!(rθ_coeff_err,[rθ - coeftable(mylm).cols[5][4],coeftable(mylm).cols[6][4] - rθ]);
    end
    ms = unique(df_aggreg.m)
    cols = ["tab:blue", "tab:orange", "tab:red"]
    axa.errorbar(ms,sqrtk_coeff, yerr = hcat(sqrtk_coeff_err...), label = L"h_d", capsize = 2., capthick = 1., c = cols[1])
    axa.errorbar(ms,cl_coeff, yerr = hcat(cl_coeff_err...), label = L"cl", capsize = 2., capthick = 1., c = cols[2])
    axa.errorbar(ms,rθ_coeff, yerr = hcat(rθ_coeff_err...), label = L"r_\theta", capsize = 2., capthick = 1., c = cols[3])
    axa.hlines(0., 1e-2, 1e0, colors="grey", linestyles = "--", label = "")

    rθ = -1. # rθ
    p = 1; θ=0.5
    mstar = 1. ./ (1 .- rθ ) .* 4 * p * θ^2/(1 + 3*p * θ^2) #critical threshold
    @show mstar
    axa.axvline(mstar,label=L"m^*",linestyle="--" )

    axa.legend()
    axa.set_xscale("log")
    axa.set_ylabel("Standardized effect on "*L"Q_{ST,s}")
    axa.set_xlabel(L"m")
    # showing r2 for m = 0.05
    # axa.text(0.1, 0.9, 
    #         L"m = 0.05: R^2 = "*@sprintf("%1.2f",r2s[2]),
    #         transform=axa.transAxes,
    #         fontsize=12)
    gcf()


    _let = ["a","b"]
    for (i,ax) in enumerate([axa,axb])
        _x = -0.2
        ax.text(_x, 1.05, _let[i],
            fontsize=12,
            fontweight="bold",
            va="bottom",
            ha="left",
            transform=ax.transAxes ,
        )
    end
    gcf()

    axa.set_facecolor("None")
    axb.set_facecolor("None")
    fig.set_facecolor("None")
    ###################
    ### annotating ####
    ##################
    fig.tight_layout()
    fig.savefig("setting2_2plots_M=$(M)_sigma_a=$(_df_aggreg.σ_α[1]).pdf",
                dpi=1200,
                bbox_inches = "tight")
    display(fig)

    ###########################
    #### printing latex table #
    ###########################
    using RegressionTables
    regtable([lms_s;lms_u]...; renderSettings = latexOutput("setting2_trait_dep_comp_coefficients_M=$(M)_sigma_a=$(_df_aggreg.σ_α[1])_$date_sim.txt"),
                print_estimator_section=false,
                regression_statistics=[:nobs,:r2],
                labels = Dict("Q_ST_u_mean" => "\$Q_{ST,u}\$",
                            "Q_ST_s_mean" => "\$Q_{ST,s}\$",
                                "sqrtk"=> "\$h_d\$",
                                "rθ"=> "\$r_\\theta\$",
                                "cl" => "\$cl\$",
                                "betas_mean" => "\$\\beta_s\$",
                                "__LABEL_STATISTIC_N__" => "Number of sim."),
                groups = repeat([@sprintf("%1.2f",df.m[1]) for df in _df_aggreg_g],2),
                number_regressions = false)
    gcf()
end