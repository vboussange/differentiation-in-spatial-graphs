#=
Plotting std effect on Q_ST_u and Q_ST_s in setting 2
for M = 7/9 to display in SI

v4 is last version as of 14-03-2022
where we show only effect for low and high migration regimes
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
include("../../code/graphs_utils/src/graphs_utils.jl")
# isdir("img") ? nothing : mkdir("img")
## graphs prop

include("../format.jl")
M = 9
# date_sim = "2022-01-14" # tend = 1000
date_sim = "2022-02-27" # tend = 3000
idx_m_toplot = [1,4]

# M = 7
# date_sim = "2022-02-26" # tend = 3000
# date_sim = "2022-01-13" # tend = 1000
# date_sim = "2022-02-18" # tend = 2000
# idx_m_toplot = [1,5]

@load "../../code/simulations/setting_2/M=$M/setting_2_mu_01_M=$(M)_hetero_2_[-onehalf,onehalf]/setting_2_mu_01_M=$(M)_hetero_2_[-onehalf,onehalf]_$(date_sim)_aggreg.jld2" df_aggreg
df_aggreg_g = groupby(df_aggreg,:m,sort=true)
totsize = size(df_aggreg_g[1],1)

fig,axs = subplots(1,2,figsize = (FIGSIZE_S[2] * 2.1, FIGSIZE_S[1]))

axa,axb = axs

# calculating GLM for each m and extracting coefficients for Q_ST_u
sqrtk_coeff = []; sqrtk_coeff_err = []; cl_coeff = []; cl_coeff_err = []; rθ_coeff = []; rθ_coeff_err = []; lms_u = []; r2s = []
for i in 1:length(df_aggreg_g)
    df_temp = DataFrame(df_aggreg_g[i][:,["Q_ST_u_mean","sqrtk", "cl", "rθ"]]); [df_temp[!,n] = _scale(df_temp[:,n]) .|> Float64 for n in names(df_temp)]
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
axb.barh(1 .+ collect(0.2 .* (1:2) .- 0.4),
        sqrtk_coeff[idx_m_toplot], 
        xerr = hcat(sqrtk_coeff_err...)[idx_m_toplot], 
        capsize = 2., 
        color = cols[1:2], 
        align="center", 
        height=0.1,
        # edgecolor="black"
        )
axb.barh(2 .+ collect(0.2 .* (1:2) .- 0.4),
        cl_coeff[idx_m_toplot], 
        xerr = hcat(cl_coeff_err...)[:,idx_m_toplot], 
        # label = "m = $(@sprintf("%1.2f",ms[1])), "*legs[i], 
        capsize = 2., 
        color = cols[[1,2]], 
        align="center", 
        height=0.1,
        # edgecolor="black"
        )

axb.barh(3 .+ collect(0.2 .* (1:2) .- 0.4),
        rθ_coeff[idx_m_toplot], 
        xerr = hcat(rθ_coeff_err...)[:,idx_m_toplot], 
        capsize = 2., 
        color = cols[[1,2]], 
        align="center", 
        height=0.1,
        # edgecolor="black"
        )
axb.vlines(0., 0.4, 3.4, colors="grey", linestyles = "--", label = "")
axb.set_yticks(1:3)
axb.set_yticklabels([L"h_d",L"\langle l \rangle", L"r_\theta"],fontsize=12)
# axb.legend()
axb.set_xlabel("Standardized effect on "*L"Q_{ST,u}")
gcf()


# calculating GLM for each m and extracting coefficients for Q_ST_s
sqrtk_coeff = []; sqrtk_coeff_err = []; cl_coeff = []; cl_coeff_err = []; rθ_coeff = []; rθ_coeff_err = []; lms_s = []; r2s = []
for i in 1:length(df_aggreg_g)
    df_temp = DataFrame(df_aggreg_g[i][:,["Q_ST_s_mean","sqrtk", "cl", "rθ"]]); [df_temp[!,n] = _scale(df_temp[:,n]) .|> Float64 for n in names(df_temp)]
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
axa.barh(1.,
        sqrtk_coeff[idx_m_toplot[1]], 
        xerr = hcat(sqrtk_coeff_err...)[idx_m_toplot[1]], 
        label = "m = $(@sprintf("%1.2f",ms[idx_m_toplot[1]]))", 
        capsize = 2., 
        color = cols[1], 
        align="center", 
        height=0.1,
        # edgecolor="black"
        )
axa.barh(1.2,
        sqrtk_coeff[idx_m_toplot[2]], 
        xerr = hcat(sqrtk_coeff_err...)[idx_m_toplot[2]], 
        label = "m = $(@sprintf("%1.2f",ms[idx_m_toplot[2]]))", 
        capsize = 2., 
        color = cols[2], 
        align="center", 
        height=0.1,
        # edgecolor="black"
        )
axa.barh(2 .+ collect(0.2 .* (1:2) .- 0.4),
        cl_coeff[idx_m_toplot], 
        xerr = hcat(cl_coeff_err...)[:,idx_m_toplot], 
        # label = "m = $(@sprintf("%1.2f",ms[1]))", 
        capsize = 2., 
        color = cols[[1,2]], 
        align="center", 
        height=0.1,
        # edgecolor="black"
        )
axa.barh(3 .+ collect(0.2 .* (1:2) .- 0.4),
        rθ_coeff[idx_m_toplot], 
        xerr = hcat(rθ_coeff_err...)[:,idx_m_toplot], 
        capsize = 2., 
        color = cols[[1,2]], 
        align="center", 
        height=0.1,
        # edgecolor="black"
        )
axa.vlines(0., 0.4, 3.4, colors="grey", linestyles = "--", label = "")
axa.set_yticks(1:3)
axa.set_yticklabels([L"h_d",L"\langle l \rangle", L"r_\theta"],fontsize=12)
axa.legend(loc="upper left")
axa.set_xlabel("Standardized effect on "*L"Q_{ST,s}")

[ax.set_xticks([-0.5,0.,0.5]) for ax in [axa,axb]]
[ax.set_xlim([-1.,1.]) for ax in [axa,axb]]

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
fig.savefig("setting2_2plots_M=$(M).pdf",
            dpi=1200,
            bbox_inches = "tight")

###########################
#### printing latex table #
###########################
# we do not want to show too many ms
# so we drop m=mstar corresponding to index 4
# if M == 7
#     df_aggreg_g = df_aggreg_g[[1,2,3,5,6]]
#     lms_s = lms_s[[1,2,3,5,6]]
#     lms_u = lms_u[[1,2,3,5,6]]
# end
using RegressionTables
regtable([lms_s[idx_m_toplot];lms_u[idx_m_toplot]]...; renderSettings = latexOutput("setting2_coefficients_M=$(M)_$(date_sim).txt"),
            print_estimator_section=false,
            regression_statistics=[:nobs,:r2],
            labels = Dict("Q_ST_u_mean" => "\$Q_{ST,u}\$",
                        "Q_ST_s_mean" => "\$Q_{ST,s}\$",
                            "sqrtk"=> "\$h_d\$",
                            "rθ"=> "\$r_\\theta\$",
                            "cl" => "\$cl\$",
                            "Q_ST_u_residuals" => "\$Q_{ST,u} - b N\$",
                            "__LABEL_STATISTIC_N__" => "Number of sim."),
            groups = repeat([@sprintf("%1.2f",df.m[1]) for df in df_aggreg_g[idx_m_toplot]],2),
            number_regressions = false)
gcf()