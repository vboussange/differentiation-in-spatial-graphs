#=
Plotting std effect on Q_ST_u in setting 1 and 2
for M= 7 with different time horizon

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
df_aggreg_t1000_set1 = load("../../../code/simulations/setting_1/M=$M/setting_1_mu_01_M=$(M)/setting_1_mu_01_M=$(M)_2022-01-09_aggreg.jld2", "df_aggreg")
df_aggreg_t2000_set1 = load("../../../code/simulations/setting_1/M=$M/setting_1_mu_01_M=$(M)/setting_1_mu_01_M=$(M)_2022-02-16_aggreg.jld2", "df_aggreg")
df_aggreg_t1000_set2 = load("../../../code/simulations/setting_2/M=$M/setting_2_mu_01_M=7_hetero_2_[-onehalf,onehalf]/setting_2_mu_01_M=7_hetero_2_[-onehalf,onehalf]_2022-01-13_aggreg.jld2", "df_aggreg")
df_aggreg_t2000_set2 = load("../../../code/simulations/setting_2/M=$M/setting_2_mu_01_M=7_hetero_2_[-onehalf,onehalf]/setting_2_mu_01_M=7_hetero_2_[-onehalf,onehalf]_2022-02-18_aggreg.jld2", "df_aggreg")


fig,axs = subplots(2,2,figsize = (FIGSIZE_S[2] * 2.1, FIGSIZE_S[1].*2))

axa = axs[1]
axb = axs[3]
axc = axs[4]

############################
######## axa ###############
############################
legs = [ L"t = 1000", L"t = 2000"]
hatch = ["","///"]

for (i,df_aggreg) in enumerate([df_aggreg_t1000_set1, df_aggreg_t2000_set1])

    df_aggreg_g = groupby(df_aggreg,:m,sort=true)[[1,5]]
    totsize = size(df_aggreg_g[1],1)

    # calculating GLM for each m and extracting coefficients for Q_ST_u
    
    sqrtk_coeff = []; sqrtk_coeff_err = []; cl_coeff = []; cl_coeff_err = []; r2s = []; lms = []
    for i in 1:length(df_aggreg_g)
        df_temp = DataFrame(df_aggreg_g[i][:,["Q_ST_u_mean","sqrtk", "cl"]]); [df_temp[!,n] = _scale(df_temp[:,n]) .|> Float64 for n in names(df_temp)]
        mylm = lm(@formula(Q_ST_u_mean ~ sqrtk + cl), df_temp)
        push!(r2s,r2(mylm))
        push!(lms,mylm)
        sqrtk = coeftable(mylm).cols[1][2]
        push!(sqrtk_coeff,sqrtk);
        push!(sqrtk_coeff_err,[sqrtk - coeftable(mylm).cols[5][2],coeftable(mylm).cols[6][2] - sqrtk]);
        cl = coeftable(mylm).cols[1][3]
        push!(cl_coeff,cl)
        push!(cl_coeff_err,[cl - coeftable(mylm).cols[5][3],coeftable(mylm).cols[6][3] - cl]);
    end
    cols = ["tab:blue", "tab:orange", "tab:red"]
    axa.barh(1 .+ collect(0.2 .- 0.4.+0.1*(i-1)),
            sqrtk_coeff[1], 
            xerr = hcat(sqrtk_coeff_err...)[1], 
            label = "m = $(@sprintf("%1.2f",ms[1])), "*legs[i], 
            capsize = 2., 
            color = cols[1], 
            align="center", 
            height=0.1,
            hatch=hatch[i], 
            edgecolor="black")
    axa.barh(1 .+ collect(0.4 .- 0.4.+0.1*(i-1)),
            sqrtk_coeff[2], 
            xerr = hcat(sqrtk_coeff_err...)[2], 
            label = "m = $(@sprintf("%1.2f",ms[2])), "*legs[i], 
            capsize = 2., 
            color = cols[2], 
            align="center", 
            height=0.1,
            hatch=hatch[i], 
            edgecolor="black")
    axa.barh(2 .+ collect(0.2 .* (1:2) .- 0.4.+0.1*(i-1)),
            cl_coeff, 
            xerr = hcat(cl_coeff_err...), 
            # label = "m = $(@sprintf("%1.2f",ms[1])), "*legs[i], 
            capsize = 2., 
            color = cols[1:2], 
            align="center", 
            height=0.1,
            hatch=hatch[i], 
            edgecolor="black")
    # axa.errorbar(ms,sqrtk_coeff, yerr = hcat(sqrtk_coeff_err...), 
    #             label = i == 1 ?  L"h_d, "*legs[i] : nothing, 
    #             capsize = 2., capthick = 1., c = cols[1], linestyle = lst[i])
    # axa.errorbar(ms,cl_coeff, yerr = hcat(cl_coeff_err...), 
    #             label = i == 1 ? L"\langle l \rangle, "*legs[i] : nothing, 
    #             capsize = 2., capthick = 1., c = cols[2], linestyle = lst[i])
end
axa.vlines(0., 0.4, 2.4, colors="grey", linestyles = "--", label = "")
# axa.legend(loc="upper center", bbox_to_anchor=(0.5, 1.6),
#             ncol=2,
#             fancybox=true,)
# axa.set_xscale("log")
axa.set_xlabel("Standardized effect on "*L"Q_{ST,u}")
axa.set_yticks(1:2)
axa.set_yticklabels([L"h_d",L"\langle l \rangle"],fontsize=12)# showing r2 for m = 0.05
axa.set_title("Setting with no selection")
gcf()


############################
######## axb ###############
############################
for (i,df_aggreg) in enumerate([df_aggreg_t1000_set2, df_aggreg_t2000_set2])

    df_aggreg_g = groupby(df_aggreg,:m,sort=true)[[1,5]]
    totsize = size(df_aggreg_g[1],1)

    # calculating GLM for each m and extracting coefficients for Q_ST_s
    sqrtk_coeff = []; sqrtk_coeff_err = []; cl_coeff = []; cl_coeff_err = []; rθ_coeff = []; rθ_coeff_err = []; lms_s = []; r2s = []
    for i in 1:length(df_aggreg_g)
        df_temp = DataFrame(df_aggreg_g[i][:,["Q_ST_u_mean","N_mean","sqrtk", "cl", "rθ"]]); [df_temp[!,n] = _scale(df_temp[:,n]) .|> Float64 for n in names(df_temp)]
        mylm = lm(@formula(Q_ST_u_mean ~ sqrtk + cl + rθ), df_temp)
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
    ms = unique(df_aggreg.m)[[1,5]]
    cols = ["tab:blue", "tab:orange", "tab:red"]
    axb.barh(1 .+ collect(0.2 .* 1 .- 0.4.+0.1*(i-1)),
            sqrtk_coeff[1], 
            xerr = hcat(sqrtk_coeff_err...)[1], 
            # label = "m = $(@sprintf("%1.2f",ms[1])), "*legs[i], 
            capsize = 2., 
            color = cols[1], 
            align="center", 
            height=0.1,
            hatch=hatch[i], 
            edgecolor="black")
    axb.barh(1 .+ collect(0.2 .* 2 .- 0.4.+0.1*(i-1)),
            sqrtk_coeff[2], 
            xerr = hcat(sqrtk_coeff_err...)[2], 
            # label = "m = $(@sprintf("%1.2f",ms[2])), "*legs[i], 
            capsize = 2., color = cols[2], 
            align="center", height=0.1,
            hatch=hatch[i], 
            edgecolor="black")
    axb.barh(2 .+ collect(0.2 .* (1:2) .- 0.4.+0.1*(i-1)),
            cl_coeff, 
            xerr = hcat(cl_coeff_err...), 
            capsize = 2., 
            color = cols[1:2], 
            align="center", 
            height=0.1,
            hatch=hatch[i], 
            edgecolor="black")
    axb.barh(3 .+ collect(0.2 .* (1:2) .- 0.4.+0.1*(i-1)),
            rθ_coeff, 
            xerr = hcat(rθ_coeff_err...), 
            capsize = 2., 
            color = cols[1:2], 
            align="center", 
            height=0.1,
            hatch=hatch[i],
            edgecolor="black")
end
axb.vlines(0., 0.4, 3.4, colors="grey", linestyles = "--", label = "")
axb.set_yticks(1:3)
axb.set_yticklabels([L"h_d",L"\langle l \rangle", L"r_\theta"],fontsize=12)
# axb.set_xscale("log")
axb.set_xlabel("Standardized effect on "*L"Q_{ST,u}")
# axb.set_xlabel(L"m")
axb.set_title("Setting with heterogeneous selection")

############################
######## axc ###############
############################
for (i,df_aggreg) in enumerate([df_aggreg_t1000_set2, df_aggreg_t2000_set2])

    df_aggreg_g = groupby(df_aggreg,:m,sort=true)[[1,5]]
    totsize = size(df_aggreg_g[1],1)

    # calculating GLM for each m and extracting coefficients for Q_ST_s
    sqrtk_coeff = []; sqrtk_coeff_err = []; cl_coeff = []; cl_coeff_err = []; rθ_coeff = []; rθ_coeff_err = []; lms_s = []; r2s = []
    for i in 1:length(df_aggreg_g)
        df_temp = DataFrame(df_aggreg_g[i][:,["Q_ST_s_mean","N_mean","sqrtk", "cl", "rθ"]]); [df_temp[!,n] = _scale(df_temp[:,n]) .|> Float64 for n in names(df_temp)]
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
    ms = unique(df_aggreg.m)[[1,5]]
    cols = ["tab:blue", "tab:orange", "tab:red"]
    axc.barh(1 .+ collect(0.2 .* 1 .- 0.4.+0.1*(i-1)),
            sqrtk_coeff[1], 
            xerr = hcat(sqrtk_coeff_err...)[1], 
            # label = "m = $(@sprintf("%1.2f",ms[1])), "*legs[i], 
            capsize = 2., 
            color = cols[1], 
            align="center", 
            height=0.1,
            hatch=hatch[i], 
            edgecolor="black")
    axc.barh(1 .+ collect(0.2 .* 2 .- 0.4.+0.1*(i-1)),
            sqrtk_coeff[2], 
            xerr = hcat(sqrtk_coeff_err...)[2], 
            # label = "m = $(@sprintf("%1.2f",ms[2])), "*legs[i], 
            capsize = 2., color = cols[2], 
            align="center", height=0.1,
            hatch=hatch[i], 
            edgecolor="black")
    axc.barh(2 .+ collect(0.2 .* (1:2) .- 0.4.+0.1*(i-1)),
            cl_coeff, 
            xerr = hcat(cl_coeff_err...), 
            capsize = 2., 
            color = cols[1:2], 
            align="center", 
            height=0.1,
            hatch=hatch[i], 
            edgecolor="black")
    axc.barh(3 .+ collect(0.2 .* (1:2) .- 0.4.+0.1*(i-1)),
            rθ_coeff, 
            xerr = hcat(rθ_coeff_err...), 
            capsize = 2., 
            color = cols[1:2], 
            align="center", 
            height=0.1,
            hatch=hatch[i],
            edgecolor="black")
end
axs[2].axis("off")
axc.vlines(0., 0.4, 3.4, colors="grey", linestyles = "--", label = "")
axc.set_yticks(1:3)
axc.set_yticklabels([L"h_d",L"\langle l \rangle", L"r_\theta"],fontsize=12)
# axb.set_xscale("log")
axc.set_xlabel("Standardized effect on "*L"Q_{ST,s}")
# axb.set_xlabel(L"m")
axc.set_title("Setting with heterogeneous selection")


fig.legend(bbox_to_anchor=(0.4, 0.4),
            ncol=1,
            fancybox=true,)


# showing r2 for m = 0.05
# axb.text(0.1, 0.9, 
#         L"m = 0.05: R^2 = "*@sprintf("%1.2f",r2s[2]),
#         transform=axb.transAxes,
#         fontsize=12)
gcf()

[ax.set_xticks([-0.5,0.,0.5]) for ax in axs]
[ax.set_xlim([-1.,1.]) for ax in axs]

_let = [L"\textbf{a}",L"\textbf{b}", L"\textbf{c}", L"\textbf{d}"]
for (i,ax) in enumerate([axa,axb,axc])
    _x = -0.2
    ax.text(_x, 1.05, _let[i],
        fontsize=12,
        fontweight="bold",
        va="bottom",
        ha="left",
        transform=ax.transAxes ,
    )
end

axa.set_facecolor("None")
axb.set_facecolor("None")
fig.set_facecolor("None")
###################
### annotating ####
##################
fig.tight_layout()
fig.savefig("time_effect_Q_ST_u.pdf",
            dpi=1200,
            bbox_inches = "tight")

fig.savefig("time_effect_QST_u_QST_s.png", dpi=500)
gcf()
