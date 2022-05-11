#= 
Plots for figure 2.

In comparision to `plot_v2.jl`, 
here we plot between population variance scaled by the total variance

v3 : only m=0.01 and m= 0.5 displayed
=#

cd(@__DIR__)
using EvoId,JLD2
using DataFrames
using Printf;#pyplot()
using Glob
using LightGraphs
using LaTeXStrings
using KernelDensity
using GLM
using Polynomials # for plotfit function
using ColorSchemes
using PyPlot
cm_eth = ColorMap([c for c in eth_grad_std.colors]);
include("../../code/graphs_utils/src/graphs_utils.jl")

idx_m_toplot = [1,4]

# isdir("img") ? nothing : mkdir("img")
## graphs prop

M = 7
date_sim = "2022-01-09"

# M = 9
# date_sim = "2022-01-13"

# suffix = "t=1000_"
suffix = ""

# @load "../../code/simulations/setting_1/M=$M/setting_1_mu_01_M=$M/setting_1_mu_01_M=$(M)_2021-12-30_aggreg.jld2" df_aggreg
@load "../../code/simulations/setting_1/M=$M/setting_1_mu_01_M=$M/setting_1_mu_01_M=$(M)_$(date_sim)_aggreg.jld2" df_aggreg
# @load "../../code/simulations/setting_1/M=$M/setting_1_mu_01_M=$(M)_dim_neutr_500/setting_1_mu_01_M=$(M)_dim_neutr_500_$(date_sim)_aggreg.jld2" df_aggreg
# @load "../../code/simulations/setting_1/M=$M/setting_1_mu_01_M=$(M)_tend_600/setting_1_mu_01_M=$(M)_tend_600_$(date_sim)_aggreg.jld2" df_aggreg
@load "../../code/graphs_utils/M=$M/graph_prop_M=$M.jld2" graphs_df


function scatter_plot_with_graph(Xlab,
                        Ylab,
                        df,
                        xlab,
                        ylab,
                        m_toplot;
                        graphs_toplot=nothing,
                        ax = nothing,
                        seed = nothing,
                        ba = 0.7,
                        pos_inset = 4,
                        )
    X = df[:,Xlab]; Y = df[:,Ylab]
    x = _scale(X)
    y = _scale(Y)
    xy =hcat(x,y)
    k = InterpKDE(kde(xy))
    c = KernelDensity.pdf.(Ref(k),xy[:,1],xy[:,2])

    mesh = ax.scatter(X,#discarding some value
                    Y,
                    c = c,
                    s = 10,
                    cmap = cm_eth
                    )
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(L"m = "*"$(m_toplot)")

    if Ylab == "bar_N_mean"
        cbaxes = mplt.inset_axes(ax, width="3%", height="30%", loc=pos_inset)
        cb = plt.colorbar(mesh, cax = cbaxes,shrink = 0.1)
        cb.set_label("graph density",fontsize = 12)
        cbaxes.yaxis.set_ticks_position("none") 
        cbaxes.yaxis.set_label_position(:left)
    end

    (xmin,xmax) = ax.get_xlim()
    (ymin,ymax) = ax.get_ylim()
    xoffset = 0.1 .* (xmax - xmin)
    yoffset = 0.1 .* (ymax - ymin)
    Ylab == "bar_N_mean" ? xoo = 0.005 : xoo = 0. #adding little offset for specific fig a
    ax.set_xlim(xmin - xoffset +xoo,xmax + xoffset + xoo)
    ax.set_ylim(ymin - yoffset,ymax + yoffset)
    ax.set_box_aspect(ba)
    #
    # gcf()
    # fig.savefig("betau_m_$(m_toplot)_density.png")
    if !isnothing(graphs_toplot)
        idx_graphs = [] # indices in graphs_toplot
        text_graphs = []
        text_graphs = []
        offsets = []
        for k in graphs_toplot
             _arg = findfirst(g -> g == graphs_df.graph[k[2][1]],df.graph)
             if !isnothing(_arg)
                 push!(idx_graphs,_arg)
                 push!(text_graphs,k[1])
                 push!(offsets,k[2][2])
             else
                 println("Could not find" , k[1])
             end
        end
        println(text_graphs)
        for (o,i) in enumerate(idx_graphs)
            x = X[i]
            y = Y[i]
            ax.scatter(x,#discarding some value
                        y,
                        edgecolors = "black",
                        c = "None",
                        s = 100,
                        cmap = cm_eth,
                        zorder = 10)
            g = to_nx(df.graph[i])
            # pos = nx.nx_agraph.graphviz_layout(g, prog="neato")
            pos = nx.spring_layout(g,seed=seed)
            dx_dy = [x,y] .+ ([xoffset, yoffset] .*offsets[o])
            for k in pos
                pos[k[1]] = k[2].* 0.06 .* [xmax-xmin ,(xmax-xmin) * ax. get_data_ratio() / ba ] .+ dx_dy
            end
            nx.draw_networkx(g,pos,
                            # edgecolors="tab:green",
                            node_size = 25.,
                            node_color = "tab:gray",
                            # linewidths = 3.,
                            # horizontalalignment = "right",
                            # verticalalignment = "baseline",
                            # alpha = 0.
                            # options,
                            with_labels = false,
                            ax = ax,
                            )
            # gcf()
            xy_ar = (x,y) .+ (dx_dy[1]-x,dx_dy[2]-y) .* 0.5
            ax.annotate("", xy=xy_ar, xytext=(x, y),arrowprops=Dict("arrowstyle" =>"<-"))
        end
    end
    ax.tick_params(left=true, bottom=true, labelleft=true, labelbottom=true)
end

function plotfit(X,Y,ax,pol, filter_left = false)
    xfit = Float64.(X); yfit = Float64.(Y)
    # p = curve_fit((t,p) ->  p[1] * exp.(-p[2] * t), xfit, yfit, [0.,0.])
    p = Polynomials.fit(xfit,yfit,pol)
    filter_left ? xeval = sort!(xfit)[2:end] : xeval = sort!(xfit)[1:end]
    ax.plot(xeval,p.(xeval),c="tab:blue")
end

include("../format.jl")
figsize = FIGSIZE_L
clf()
fig = plt.figure(constrained_layout=true,
                # figsize = (FIGSIZE_S[2] * 2.1, 2*FIGSIZE_S[1])
                figsize = figsize   
                )
width_ratios = [0.5, 0.5]
height_ratios = [0.5,0.5]
gs = fig.add_gridspec(2, 2, width_ratios = width_ratios, height_ratios = height_ratios)
ax1 = fig.add_subplot(py"$(gs)[0,0]")
ax2 = fig.add_subplot(py"$(gs)[0,1]")
ax3 = fig.add_subplot(py"$(gs)[1,0]")
ax4 = fig.add_subplot(py"$(gs)[1,1]")
gcf()

df_aggreg[!,:bar_N_mean] = df_aggreg.N_mean ./ M
df_aggreg[!,:bar_N_std] = df_aggreg.N_std ./ M

df_aggreg_g = groupby(df_aggreg,:m,sort=true)
totsize = size(df_aggreg_g[1],1)
# caculating average population size

### first plots
# dict of graphs to plot, with index from graphs_df and offset
graphs_toplot = Dict("star" => (1, [-1.5,0.]),
                "lollipop" => (2, [-1.2,1.4]),
                "line" => (225, [0.0,-2.]),
                "flake" => (216,[1.,-1.5]),
                "river 1" => (4,[.5,1.5]),
                "river 2" => (9,[2.,0.]),
                "complete" => (totsize,[0.,1.5]), #line
                # "bipartite 25" => 19,
                # "bipartite 34" => 212,
                # "lattice 4" => 721,
                "ring" => (422,[1,-1.])
                )
# only plotting graph for M=7 on main manuscript
M > 7 ? graphs_toplot = nothing : nothing;
scatter_plot_with_graph("cl",
                "Q_ST_u_mean",
                df_aggreg_g[1],
                L"\langle l \rangle",
                L"Q_{ST,u}",
                0.01,
                graphs_toplot = graphs_toplot,
                seed = 1,
                ba = height_ratios[1] * figsize[2] / (width_ratios[2] * figsize[1]) ,
                ax = ax2)
plotfit(df_aggreg_g[1][:,"cl"],df_aggreg_g[1].Q_ST_u_mean,ax2,1)
df_temp = df_aggreg_g[1][:,["Q_ST_u_mean","cl"]] .|> Float64
ax2.text(0.05, 0.95, 
        L"R^2 = "*@sprintf("%1.2f",r2(lm(@formula(Q_ST_u_mean ~ cl), df_temp))),
        transform=ax2.transAxes,
        fontsize=12)
(ymin,ymax) = ax2.get_ylim()
ax2.set_ylim(ymin+0.025,ymax+0.025)
gcf()

### first plots
# dict of graphs to plot, with index from graphs_df and offset
graphs_toplot = Dict("star" => (1, [2.,0,]),
                "lollipop" => (2, [-2.5,0.5]),
                "line" => (225, [-0.5,1.]),
                # "flake" => (216,[-1.,1.5]),
                "river 1" => (4,[-1.5,0.5]),
                "river 2" => (9,[1.,-1.]),
                "complete" => (totsize,[.8,.8]), #line
                # "bipartite 25" => 19,
                # "bipartite 34" => 212,
                # "lattice 4" => 721,
                "ring" => (422,[0.,-2.])
                )
M > 7 ? graphs_toplot = nothing : nothing;
scatter_plot_with_graph("sqrtk",
                "bar_N_mean",
                df_aggreg_g[4],
                L"h_d",
                L"\bar{N}",
                0.5,
                graphs_toplot = graphs_toplot,
                seed = 1,
                ax = ax1,
                ba= height_ratios[1] * figsize[2] / (width_ratios[1]* figsize[1])  )
npop(sqrtk) = sqrtk * 150 
ax1.plot(df_aggreg_g[4].sqrtk,npop.(df_aggreg_g[4].sqrtk),c="tab:blue")

df_temp = df_aggreg_g[4][:,["N_mean","sqrtk"]] .|> Float64
ax1.text(0.05, 0.95,
        L"R^2 = "*@sprintf("%1.2f",r2(lm(@formula(N_mean ~ sqrtk), df_temp))),
        transform=ax1.transAxes,
        fontsize=12)
gcf()
graphs_toplot = Dict("star" => (1, [0.,-2.]),
                "lollipop" => (2, [1.,.6]),
                "line" => (225, [1.5,0.]),
                "flake" => (216,[1.5,1.5]),
                "river 1" => (4,[-1.5,-1.5]),
                "river 2" => (9,[1.,1.]),
                "complete" => (totsize,[-2.,-1.]), #line
                # "bipartite 25" => 19,
                # "bipartite 34" => 212,
                # "lattice 4" => 721,
                "ring" => (422,[0.5,2.])
                )
M > 7 ? graphs_toplot = nothing : nothing;
scatter_plot_with_graph("sqrtk",
                "Q_ST_u_mean",
                df_aggreg_g[1],
                L"h_d",
                L"Q_{ST,u}",
                0.01,
                graphs_toplot = graphs_toplot,
                seed = 1,
                ax = ax3,
                ba= height_ratios[1] * figsize[2] / (width_ratios[1]* figsize[1]),
                pos_inset=8)
plotfit(df_aggreg_g[1].sqrtk,df_aggreg_g[1].Q_ST_u_mean,ax3,1, true)

df_temp = df_aggreg_g[1][:,["Q_ST_u_mean","sqrtk"]] .|> Float64
ax3.text(0.05, 0.95,
        L"R^2 = "*@sprintf("%1.2f",r2(lm(@formula(Q_ST_u_mean ~ sqrtk), df_temp))),
        transform=ax3.transAxes,
        fontsize=12)
gcf()

# calculating GLM for each m and extracting coefficients
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
ms = unique(df_aggreg.m)
cols = ["tab:blue", "tab:orange", "tab:red"]
ax4.barh(0.8,
        sqrtk_coeff[idx_m_toplot[1]], 
        xerr = hcat(sqrtk_coeff_err...)[idx_m_toplot[1]], 
        label = "m = $(@sprintf("%1.2f",ms[idx_m_toplot[1]]))", 
        capsize = 2., 
        color = cols[1], 
        align="center", 
        height=0.1,
        # edgecolor="black"
        )
ax4.barh(1.0,
        sqrtk_coeff[idx_m_toplot[2]], 
        xerr = hcat(sqrtk_coeff_err...)[idx_m_toplot[2]],
        label = "m = $(@sprintf("%1.2f",ms[idx_m_toplot[2]]))", 
        capsize = 2., 
        color = cols[2], 
        align="center", 
        height=0.1,
        # edgecolor="black"
        )
ax4.barh(2 .+ collect(0.2 .* (1:2) .- 0.4),
        cl_coeff[idx_m_toplot], 
        xerr = hcat(cl_coeff_err...)[:,idx_m_toplot], 
        # label = "m = $(@sprintf("%1.2f",ms[1])), "*legs[i], 
        capsize = 2., 
        color = cols[[1,2]], 
        align="center", 
        height=0.1,
        # edgecolor="black"
        )
ax4.legend(loc = "upper left")
ax4.vlines(0., 0.4, 2.4, colors="grey", linestyles = "--", label = "")
ax4.set_yticks(1:2)
ax4.set_yticklabels([L"h_d",L"\langle l \rangle"],fontsize=12)
ax4.set_xticks([-0.5,0.,0.5])
ax4.set_xlim([-1.,1.])
# axb.legend()
ax4.set_xlabel("Standardized effect on "*L"Q_{ST,u}")
gcf()

################
##### 3d figure ##
################

_let = [L"\textbf{a}",L"\textbf{b}", L"\textbf{c}", L"\textbf{d}"]
for (i,ax) in enumerate([ax1,ax2,ax3,ax4])
    _x = -0.2
    ax.text(_x, 1.05, _let[i],
        fontsize=12,
        va="bottom",
        ha="left",
        transform=ax.transAxes ,
    )
    # ax.set_xticks(-1.:0.5:1.)
end
gcf()

ax3.set_facecolor("None")
ax2.set_facecolor("None")
ax1.set_facecolor("None")
fig.set_facecolor("None")
###################
### annotating ####
##################
fig.savefig("setting1_neutr_$(suffix)M=$M.pdf",
            dpi=1200,
            bbox_inches = "tight",
            )
display(fig)

###########################
#### printing latex table #
###########################
## 2 variate models
if true
    using RegressionTables
    regtable(lms[idx_m_toplot]...; renderSettings = latexOutput("setting1_neutr_$(suffix)M=$M.txt"),
                    print_estimator_section=false,
                    regression_statistics=[:nobs,:r2],
                    labels = Dict("Q_ST_u_mean" => "\$Q_{ST,u}\$",
                                "Q_ST_s_mean" => "\$Q_{ST,s}\$",
                                    "sqrtk"=> "\$h_d\$",
                                    "rθ"=> "\$r_\\theta\$",
                                    "cl" => "\$cl\$",
                                    "Q_ST_u_residuals" => "\$Q_{ST,u} - b \bar{N}\$",
                                    "__LABEL_STATISTIC_N__" => "Number of sim."),
                    groups = [@sprintf("%1.2f",df.m[1]) for df in df_aggreg_g[idx_m_toplot]],
                    number_regressions = false)

    gcf()

    ## 1 variate models
    ## characteristic length
    lms1var = []
    for i in (1:length(df_aggreg_g))[idx_m_toplot]
        df_temp = DataFrame(df_aggreg_g[i][:,["Q_ST_u_mean", "cl"]]); [df_temp[!,n] = _scale(df_temp[:,n]) .|> Float64 for n in names(df_temp)]
        mylm = lm(@formula(Q_ST_u_mean ~ cl), df_temp)
        push!(lms1var,mylm)
    end
    ## sqrtk
    for i in (1:length(df_aggreg_g))[idx_m_toplot]
        df_temp = DataFrame(df_aggreg_g[i][:,["Q_ST_u_mean", "sqrtk"]]); [df_temp[!,n] = _scale(df_temp[:,n]) .|> Float64 for n in names(df_temp)]
        mylm = lm(@formula(Q_ST_u_mean ~ sqrtk), df_temp)
        push!(lms1var,mylm)
    end
    ## sqrtk corrected
    for i in (1:length(df_aggreg_g))[idx_m_toplot]
        df_temp = DataFrame(df_aggreg_g[i][:,["Q_ST_u_mean","N_mean", "sqrtk"]]); [df_temp[!,n] = _scale(df_temp[:,n]) .|> Float64 for n in names(df_temp)]
        mylmres = lm(@formula(Q_ST_u_mean ~ N_mean), df_temp)
        df_temp[!,"Q_ST_u_residuals"] = residuals(mylmres)
        mylm = lm(@formula(Q_ST_u_residuals ~ sqrtk), df_temp)
        push!(lms1var,mylm)
    end
    regtable(lms1var...; renderSettings = latexOutput("setting1_neutr_$(suffix)M=$(M)_1var.txt"),
                print_estimator_section=false,
                regression_statistics=[:nobs,:r2],
                labels = Dict("Q_ST_u_mean" => "\$Q_{ST,u}\$",
                            "Q_ST_s_mean" => "\$Q_{ST,s}\$",
                                "sqrtk"=> "\$h_d\$",
                                "rθ"=> "\$r_\\theta\$",
                                "cl" => "\$cl\$",
                                "Q_ST_u_residuals" => "\$Q_{ST,u} - b \\bar{N}\$",
                                "__LABEL_STATISTIC_N__" => "Number of sim."),
                groups = repeat([@sprintf("%1.2f",df.m[1]) for df in df_aggreg_g[idx_m_toplot]],3),
                number_regressions = false)
end


###########################
#### extra analysis #######
###########################
#= 
although we find that 
there are some discrepancies looking at the 
individual level for graphs by taking different time horizons
we find that this can be explained by the large variance 
in graphs with low number of vertices
=#
# df_aggreg[df_aggreg.graph .== Ref(graphs_df.graph[graphs_toplot["line"][1]]),:]
# df_aggreg[df_aggreg.graph .== Ref(graphs_df.graph[graphs_toplot["lollipop"][1]]),:]
# df_aggreg[df_aggreg.graph .== Ref(graphs_df.graph[graphs_toplot["star"][1]]),:]