#=
Plotting setting 2 adapt against PDE meanfield approximation

Here results are taken from 
the IBM model simulation `df_ibm`
=#
cd(@__DIR__)
using LightGraphs
using DataFrames
using CSV
using GLM
using Statistics
using JLD2
using EvoId
using Printf
using Random
using Polynomials
using DataStructures
using KernelDensity
using PyPlot
cm_eth = ColorMap([c for c in eth_grad_std.colors]);
include("../../code/graphs_utils/src/graphs_utils.jl")
include("../format.jl")

M = 7
date_sim = "2022-02-26"
@load "../../code/simulations/setting_2/M=$M/setting_2_mu_01_M=$(M)_hetero_2_[-onehalf,onehalf]/setting_2_mu_01_M=$(M)_hetero_2_[-onehalf,onehalf]_$(date_sim)_aggreg.jld2" df_aggreg

@load "../../code/graphs_utils/M=7/graph_prop_M=7.jld2" graphs_df
@load "../../code/simulations/setting_2/explo_r&m-pde/pde_meanfield_data_m_01-05.jld2" df_explo

# selecting migration regime
df_aggreg_g = groupby(df_aggreg,"m",sort=true)
df_explo_g = groupby(df_explo,"m",sort=true)

df_ibm = df_aggreg_g[3] |> DataFrame # 3 for m = 0.1, 5 for m = 0.5
df_pde = df_explo_g[1] |> DataFrame  # 1 for m = 0.1, 2 for m = 0.5
m = df_ibm.m[1]

# adding bar{N}
df_ibm.bar_N = df_ibm.N_mean ./ 7
df_pde.bar_N = df_pde.npop ./ 7

# loading illustrative graphs
date_sim = "2022-01-19"
@load "../../code/simulations/setting_2/M=$M/setting_2_mu_01_M=$(M)_hetero_2_illustrative_graphs/setting_2_mu_01_M=$(M)_hetero_2_illustrative_graphs_$(date_sim)_aggreg.jld2" df_aggreg
df_aggreg_g = groupby(df_aggreg,"m",sort=true)

df_ibm_illustr = df_aggreg_g[1] |> DataFrame  # 1 for m = 0.1, 3 for m = 0.5
sort!(df_ibm_illustr,:rθ)
df_ibm_illustr.bar_N = df_ibm_illustr.N_mean ./ 7

# making sure that pde and ibm data are from similar migration regime
@assert df_ibm.m[1] ≈ df_pde.m[1]
@assert df_ibm.m[1] ≈ df_ibm_illustr.m[1]

# merging df_ibm_illustr and df_ibm, making sure that graphs are not counted twice
for r in eachrow(df_ibm_illustr)
    _arg = findfirst((r.rθ .== df_ibm.rθ) .& (Ref(r.graph) .== df_ibm.graph))
    if isnothing(_arg)
        println("adding graph", r.graph)
        append!(df_ibm, r |>DataFrame)
    else #replacing values, to make sure that df_ibm and df_ibm_illustr are coherent
        df_ibm[_arg,:] = r
    end
end

function plot_scatter(myx::Symbol,
                        myy::Symbol,
                        xlab,
                        ylab,
                        m_toplot;
                        graphs_toplot=nothing,
                        offsets=nothing,
                        ax = nothing,
                        seed = nothing,
                        ba = 0.7
                        )
    x = _scale(df_ibm[:,myx])
    y = _scale(df_ibm[:,myy])
    xy =hcat(x,y)
    k = InterpKDE(kde(xy))
    c = KernelDensity.pdf.(Ref(k),xy[:,1],xy[:,2])

    isnothing(ax) ? ((fig, ax) = plt.subplots()) : nothing
    mesh = ax.scatter(df_ibm[:,myx],#discarding some value
                    df_ibm[:,myy],
                    c = c,
                    s = 10,
                    cmap = cm_eth)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(L"m = "*"$(m_toplot)")
    ax.set_box_aspect(ba)
    #
    # gcf()
    # fig.savefig("betau_m_$(m_toplot)_density.png")
    return ax, mesh
end

function plotting_graph(df,
                        myx::Symbol,
                        myy::Symbol;
                        ax = nothing,
                        seed = nothing,
                        ba = 0.7
                        )
    (xmin,xmax) = ax.get_xlim()
    (ymin,ymax) = ax.get_ylim()
    xoffset = 0.1 .* (xmax - xmin)
    yoffset = 0.1 .* (ymax - ymin)
    scale_gr = 0.17

    for r in eachrow(df)
        x = r[myx]
        y = r[myy]

        ax.scatter(x,#discarding some value
                    y,
                    edgecolors = "black",
                    c = "None",
                    s = 50,
                    cmap = cm_eth,
                    zorder = 10)
        g = to_nx(r.graph)
        # pos = nx.nx_agraph.graphviz_layout(g, prog="neato")
        pos = nx.spring_layout(g,seed=seed)
        offsets = r.offsets
        dx_dy = [x,y] .+ ([xoffset, yoffset] .* offsets)
        for k in pos
            pos[k[1]] = k[2].* scale_gr .* [1 , ax.get_data_ratio() / ba] .+ dx_dy
        end
        # iterating over red and blues nodes
        _c = [r.soptim[i] < 0 ? "tab:blue" : "tab:red" for i in 1:length(r.soptim)]
        nx.draw_networkx(g,pos,
                        # nodelist = node_list,
                        edgecolors = _c,
                        node_size = 25.,
                        node_color = _c,
                        # linewidths = 3.,
                        # horizontalalignment = "right",
                        # verticalalignment = "baseline",
                        # alpha = 0.
                        # options,
                        with_labels = false,
                        ax = ax,
                        )
        xy_ar = (x,y) .+ (dx_dy[1]-x,dx_dy[2]-y) .* 0.55
        ax.annotate("", xy=xy_ar, xytext=(x, y),arrowprops=Dict("arrowstyle" =>"<-"))
    end
        # gcf()
    ax.tick_params(left=true, bottom=true, labelleft=true, labelbottom=true)
    ## if vertical
    # ax.set_xlim(xmin - xoffset,xmax + xoffset)
    # ax.set_ylim(ymin - yoffset,ymax + yoffset)
    ## if horizontal

    ax.set_xlim(xmin - xoffset*1.7,xmax + xoffset)
    myy == :N ? yyoff = 0. : yyoff = 1.
    ax.set_ylim(ymin - yoffset,ymax + yoffset.* yyoff)
    return ax, mesh
end



#####################################
####### plotting ############
####################################
# adding offsets for plotting
df_ibm_illustr[!,"offsets"] = [[-1.0,0.],[1.3,-1.3],[-1.5,1.5],[-1.5,1.5],[1.,-2.],[0.,-1.5],[0.5,-1.5],[-2.5,0.],] .* 1.3 #last is river

## if vertical
# fig,ax = subplots(2,1,figsize = (FIGSIZE_S[1], FIGSIZE_S[2] * 2),sharex = true)
## if horizontal
fig,ax = subplots(1,2,figsize = (FIGSIZE_S[2] * 2.1, FIGSIZE_S[1]))
ax1 = ax[1]
ax3 = ax[2]

cbar3 = mplt.inset_axes(ax1,width = "3%", height = "23%",loc=4,bbox_to_anchor = (0., 0.05, 1, 1),bbox_transform=ax1.transAxes )
# cbar4 = mplt.inset_axes(ax3,width = "3%", height = "23%",loc=4,bbox_to_anchor = (0., 0.05, 1, 1),bbox_transform=ax3.transAxes )

graphs_toplot = OrderedDict("star" => 1,
                # "lollipop" => 2,
                "line" => 225,
                # "flake" => 216,
                "river 1" => 4,
                # "river 2" => 9,
                "complete" => 853,
                # "bipartite 25" => 19,
                # "bipartite 34" => 212,
                # "lattice 4" => 721,
                "ring" => 422
                )

_,_mesh = plot_scatter(:rθ,
                        :Q_ST_s_mean,
                        L"r_\theta",
                        L"Q_{ST,s}",
                        @sprintf("%1.1f",m),
                        ax = ax1,
                        seed = 1,
                        ba = 1.)
# only plotting for m=0.1 to display in main text
if m < 0.3
    plotting_graph(df_ibm_illustr,
                            :rθ,
                            :Q_ST_s_mean,
                            ax = ax1,
                            seed = 1,
                            ba = 1.)
end
gcf()
#####################
### from pde data ###
#####################
ax1.plot(df_pde[:,"rθ"],df_pde[:,"β"] ./ (df_pde[:,"α"] + df_pde[:,"β"]))
#####################
### from IBM      ###
#####################
# plotfit(df_ibm.ass,df_ibm.Q_ST_s_mean,ax1,2)
_cb = fig.colorbar(_mesh,cax = cbar3)
_cb.set_label("graph density",fontsize =8)
# _cb.ax.set_xlabel("graph_density",labelpad = 5)
# _cb.ax.xaxis.set_label_position("top")
cbar3.yaxis.set_ticks_position("none") 
cbar3.yaxis.set_label_position("left")
ax1.set_xlabel(L"r_\theta")
gcf()

#####################################
####### Population size, rtheta ############
####################################
_,_mesh = plot_scatter(:rθ,
                        :bar_N,
                        L"r_\theta",
                        L"\bar{N}",
                        @sprintf("%1.1f",m),
                        ax = ax3,
                        seed = 1,
                        ba = 1.)
# only plotting for m=0.1 to display in main text
if m < 0.3
    plotting_graph(df_ibm_illustr,
                            :rθ,
                            :bar_N,
                            ax = ax3,
                            seed = 1,
                            ba = 1.)
end
#####################
### from pde data ###
#####################
ax3.plot(df_pde[:,"rθ"],df_pde[:,"bar_N"] /2. * 7.)
gcf()
#####################
### from IBM      ###
#####################

# _cb = fig.colorbar(_mesh,cax = cbar4)
# _cb.set_label("graph dens.",fontsize =8)
# # _cb.ax.set_xlabel("graph_density",labelpad = 5)
# # _cb.ax.xaxis.set_label_position("top")
# cbar4.yaxis.set_ticks_position("left")
# cbar4.yaxis.set_label_position("left")
# cbar4.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))


ax3.set_title("")
gcf()

fig.set_facecolor("None")

_let = [L"\textbf{a}",L"\textbf{b}", L"\textbf{c}", L"\textbf{d}"]
for (i,ax) in enumerate([ax1,ax3])
    _x = -0.2
    ax.text(_x, 1.0, _let[i],
        fontsize=12,
        fontweight="bold",
        va="bottom",
        ha="left",
        transform=ax.transAxes ,
    )
    ax.set_xticks(-1.:0.5:1.)
end

fig.tight_layout()
display(fig)

fig.savefig("sett2_adapt_IBM_horizontal_$(@sprintf("%1.1f",m)).pdf",
            dpi=1200,
            # bbox_inches = "tight",
            )
# fig.set_facecolor("None")
# ax1.set_facecolor("None")
# ax3.set_facecolor("None")
# ax3.set_facecolor("None")
# gcf()


# fig.savefig("slice_N_mu_$(m_toplot)_rtheta.png")
