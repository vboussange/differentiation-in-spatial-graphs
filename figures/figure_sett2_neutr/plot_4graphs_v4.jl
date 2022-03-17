#=
Plotting last figure in results
For setting 2

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

using CSV
# data for figure a
df_mresponse_sett2 = JLD2.load("../../code/simulations/setting_2/M=9/setting_2_mu_01_M=9_complete_hetero_2_[-onehalf,onehalf]/setting_2_mu_01_M=9_complete_hetero_2_[-onehalf,onehalf]_2022-02-21_aggreg.jld2", "df_aggreg")
df_mresponse_sett1 = JLD2.load("../../code/simulations/setting_1/M=9/setting_1_mu_01_M=9_complete/setting_1_mu_01_M=9_complete_2022-02-21_aggreg.jld2", "df_aggreg")

# data for figure b
# df_line_rtheta = JLD2.load("../../code/simulations/setting_2/M=9/setting_2_mu_01_M=9_ring_hetero_2_[-onehalf,onehalf]/setting_2_mu_01_M=9_ring_hetero_2_[-onehalf,onehalf]_2022-02-21_aggreg.jld2", "df_aggreg")
df_line_rtheta = JLD2.load("../../code/simulations/setting_2/M=9/setting_2_mu_01_M=9_line_hetero_2_[-onehalf,onehalf]/setting_2_mu_01_M=9_line_hetero_2_[-onehalf,onehalf]_2022-02-21_aggreg.jld2", "df_aggreg")
df_line_rtheta_g = groupby(df_line_rtheta,:m,sort = true)
# data figure c,d
M = 7
date_sim = "2022-02-26"
@load "../../code/simulations/setting_2/M=$M/setting_2_mu_01_M=$(M)_hetero_2_[-onehalf,onehalf]/setting_2_mu_01_M=$(M)_hetero_2_[-onehalf,onehalf]_$(date_sim)_aggreg.jld2" df_aggreg
df_aggreg_g = groupby(df_aggreg,:m,sort=true)
totsize = size(df_aggreg_g[1],1)

_norm = matplotlib.colors.Normalize(vmin=-1.5, vmax=0., clip=true)
mapper = matplotlib.cm.ScalarMappable(norm=_norm, cmap=cm_eth)

function plot_sett1_sett2(ax)

        # cmap = plt.get_cmap("tab10")

        ax.set_xlabel(L"m");
        ax.set_ylabel(L"Q_{ST}");
        ax.set_yscale("log")

        ax.errorbar(df_mresponse_sett1.m,
                df_mresponse_sett1[:,"Q_ST_u_mean"],
                yerr = df_mresponse_sett1[:,"Q_ST_u_std"],
                c = mapper.to_rgba(log10(0.04)),
                fmt = "o",
                label = L"Q_{ST,u}"*", sett.\n without selection",
                ms = 3.)
        ax.errorbar(df_mresponse_sett2.m,
                df_mresponse_sett2.Q_ST_u_mean,
                yerr = df_mresponse_sett2.Q_ST_u_std,
                c = mapper.to_rgba(log10(0.1)),
                fmt = "o",
                label = L"Q_{ST,u}"*", sett.\n with selection",
                ms = 3.)
        # axa = ax.twinx()
        # ax.set_xscale("log")

        ax.errorbar(df_mresponse_sett2.m,
                df_mresponse_sett2.Q_ST_s_mean,
                yerr = df_mresponse_sett2.Q_ST_s_std,
                c = mapper.to_rgba(log10(1.0)),
                fmt = "o",
                label = L"Q_{ST,s}"*", sett.\n with selection",
                ms = 3.)
end

function plot1(ysym,ylabel,axb;
                        ax = nothing,
                        legend = nothing,
                        g = nothing,
                        soptim_df = nothing,
                        seed = 1,
                        bboxx = -0.3)

    # skipping m=5e-2
    for df in df_line_rtheta_g[[1,3,4]]
        _m = df.m[1]
        x = Float64.(df.rθ); y = Float64.(df[:,ysym*"_mean"])
        p = Polynomials.fit(x,y,1)
        axb.errorbar(x,
                y,
                yerr = df[:,ysym*"_std"],
                c = mapper.to_rgba(log10(_m)),
                label = "m = "*(@sprintf "%2.2f" _m),
                fmt = "o",
                ms = 4.)
        axb.plot(x,
                    p.(x),
                    c = mapper.to_rgba(log10(_m)),
                    alpha = 0.5
                    )
    end
    axb.set_xlabel(L"r_\theta")
    axb.set_ylabel(ylabel)
    if !isnothing(legend)
        axb.legend(#bbox_to_anchor=(bboxx,0.6)
        )
    end
    if !isnothing(ax)
        (xmin,xmax) = axb.get_xlim()
        (ymin,ymax) = axb.get_ylim()
        yoffset = 0.1 .* (ymax - ymin)
        xoffset = 0.05 .* (xmax - xmin)

        y = ymin -yoffset
        g = to_nx(g)
        for i in 1:size(soptim_df,1)
            x = soptim_df.ass[i]
            soptim = soptim_df.soptim[i]
                # pos = nx.nx_agraph.graphviz_layout(g, prog="neato")
            pos = nx.spring_layout(g,seed=seed)
            dx_dy = [x,0] #.+ ([xoffset, yoffset])
            for k in pos
                pos[k[1]] = k[2] .* xoffset .+ dx_dy
            end
            for h in [-1,1]
                _c = h == -1 ? "tab:blue" : "tab:red"
                node_list = [k-1 for (k,v) in soptim if v == h]
                nx.draw_networkx(g,pos,
                                nodelist = node_list,
                                edgecolors = _c,
                                node_size = 25.,
                                node_color = _c,
                                with_labels = false,
                                ax = ax,
                                )
            end
        end
        ax.axis("off")
        ax.set_ylim([1.2*y for y in ax.get_ylim()])
        axb.tick_params(left=true, bottom=true, labelleft=true, labelbottom=true)
    end
    axb.set_xticks(unique(df_line_rtheta.rθ))
    axb.xaxis.set_major_formatter(plt.FormatStrFormatter("%.2f"))
    # ax.set_xlim(xmin,xmax )
end

figsize = FIGSIZE_L
clf()
fig = plt.figure(constrained_layout=true,
                # figsize = (FIGSIZE_S[2] * 2.1, 2*FIGSIZE_S[1])    
                figsize = figsize  
                )
width_ratios = [0.5, 0.5]
height_ratios = [0.5,0.5]
gs = fig.add_gridspec(2, 2, width_ratios = width_ratios, height_ratios = height_ratios)
axa = fig.add_subplot(py"$(gs)[0,0]")
axb = fig.add_subplot(py"$(gs)[0,1]")
axc = fig.add_subplot(py"$(gs)[1,0]")
axd = fig.add_subplot(py"$(gs)[1,1]")
gcf()

plot1("Q_ST_u",L"Q_{ST,u}", axb,
        # ax=axx1,
        # soptim_df=soptim_df,
        # g=g,
        seed=53,
        legend = true,
        bboxx = -0.15)

rθ = df_mresponse_sett2.rθ[1] # rθ
p = 1; θ=0.5
mstar = 1. ./ (1 .- rθ ) .* 4 * p * θ^2/(1 + 3*p * θ^2) #critical threshold
plot_sett1_sett2(axa)
axa.axvline(mstar,label=L"m^*",linestyle="--" )
axa.legend(framealpha = 1.)
gcf()

# calculating GLM for each m and extracting coefficients for Q_ST_u
sqrtk_coeff = []; sqrtk_coeff_err = []; cl_coeff = []; cl_coeff_err = []; rθ_coeff = []; rθ_coeff_err = []; lms = []; r2s = []
for i in 1:length(df_aggreg_g)
    df_temp = DataFrame(df_aggreg_g[i][:,["Q_ST_u_mean","sqrtk", "cl", "rθ"]]); [df_temp[!,n] = _scale(df_temp[:,n]) .|> Float64 for n in names(df_temp)]
    mylm = lm(@formula(Q_ST_u_mean ~ sqrtk + cl + rθ), df_temp)
    push!(lms, mylm)
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
axd.barh(1.,
        sqrtk_coeff[1], 
        xerr = hcat(sqrtk_coeff_err...)[1], 
        label = "m = $(@sprintf("%1.2f",ms[1]))", 
        capsize = 2., 
        color = cols[1], 
        align="center", 
        height=0.15,
        # edgecolor="black"
        )
axd.barh(1.2,
        sqrtk_coeff[5], 
        xerr = hcat(sqrtk_coeff_err...)[5], 
        label = "m = $(@sprintf("%1.2f",ms[5]))", 
        capsize = 2., 
        color = cols[2], 
        align="center", 
        height=0.15,
        # edgecolor="black"
        )
axd.barh(2 .+ collect(0.2 .* (1:2) .- 0.4),
        cl_coeff[[1,5]], 
        xerr = hcat(cl_coeff_err...)[:,[1,5]], 
        # label = "m = $(@sprintf("%1.2f",ms[1])), "*legs[i], 
        capsize = 2., 
        color = cols[[1,2]], 
        align="center", 
        height=0.15,
        # edgecolor="black"
        )

axd.barh(3 .+ collect(0.2 .* (1:2) .- 0.4),
        rθ_coeff[[1,5]], 
        xerr = hcat(rθ_coeff_err...)[:,[1,5]], 
        capsize = 2., 
        color = cols[[1,2]], 
        align="center", 
        height=0.15,
        # edgecolor="black"
        )
axd.vlines(0., 0.4, 3.4, colors="grey", linestyles = "--", label = "")
axd.set_yticks(1:3)
axd.set_yticklabels([L"h_d",L"\langle l \rangle", L"r_\theta"],fontsize=12)
# axd.legend()
axd.set_xlabel("Standardized effect on "*L"Q_{ST,u}")
gcf()


# calculating GLM for each m and extracting coefficients for Q_ST_s
sqrtk_coeff = []; sqrtk_coeff_err = []; cl_coeff = []; cl_coeff_err = []; rθ_coeff = []; rθ_coeff_err = []; lms = []; r2s = []
for i in 1:length(df_aggreg_g)
    df_temp = DataFrame(df_aggreg_g[i][:,["Q_ST_s_mean","sqrtk", "cl", "rθ"]]); [df_temp[!,n] = _scale(df_temp[:,n]) .|> Float64 for n in names(df_temp)]
    mylm = lm(@formula(Q_ST_s_mean ~ sqrtk + cl + rθ), df_temp)
    push!(lms, mylm)
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
axc.barh(1.,
        sqrtk_coeff[1], 
        xerr = hcat(sqrtk_coeff_err...)[1], 
        # label = "m = $(@sprintf("%1.2f",ms[1]))", 
        capsize = 2., 
        color = cols[1], 
        align="center", 
        height=0.15,
        # edgecolor="black"
        )
axc.barh(1.2,
        sqrtk_coeff[5], 
        xerr = hcat(sqrtk_coeff_err...)[5], 
        label = "m = $(@sprintf("%1.2f",ms[5]))", 
        capsize = 2., 
        color = cols[2], 
        align="center", 
        height=0.15,
        # edgecolor="black"
        )
axc.barh(2 .+ collect(0.2 .* (1:2) .- 0.4),
        cl_coeff[[1,5]], 
        xerr = hcat(cl_coeff_err...)[:,[1,5]], 
        label = "m = $(@sprintf("%1.2f",ms[1]))", 
        capsize = 2., 
        color = cols[[1,2]], 
        align="center", 
        height=0.15,
        # edgecolor="black",
        )

axc.barh(3 .+ collect(0.2 .* (1:2) .- 0.4),
        rθ_coeff[[1,5]], 
        xerr = hcat(rθ_coeff_err...)[:,[1,5]], 
        capsize = 2., 
        color = cols[[1,2]], 
        align="center", 
        height=0.15,
        # edgecolor="black",
        )
axc.vlines(0., 0.4, 3.4, colors="grey", linestyles = "--", label = "")
axc.set_yticks(1:3)
axc.set_yticklabels([L"h_d",L"\langle l \rangle", L"r_\theta"],fontsize=12)
axc.legend(loc="upper left")
axc.set_xlabel("Standardized effect on "*L"Q_{ST,s}")

[ax.set_xticks([-0.5,0.,0.5]) for ax in [axc,axd]]
[ax.set_xlim([-1.,1.]) for ax in [axc,axd]]

gcf()



_let = ["a","b","c","d"]
for (i,ax) in enumerate([axa,axb,axc,axd])
    _x = -0.2
    ax.text(_x, 1.05, _let[i],
        fontsize=12,
        fontweight="bold",
        va="bottom",
        ha="left",
        transform=ax.transAxes ,
    )
end

axc.set_facecolor("None")
axd.set_facecolor("None")
axa.set_facecolor("None")
axb.set_facecolor("None")
fig.set_facecolor("None")
###################
### annotating ####
##################
fig.tight_layout()
fig.savefig("setting2_4plots_M=$(M)_v4.pdf",
            dpi=1200,
            bbox_inches = "tight")
gcf()