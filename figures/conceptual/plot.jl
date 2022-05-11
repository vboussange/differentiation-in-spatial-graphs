# /!\ for some reason, only works in atom
using Random
cd(@__DIR__)
using Dates
using JLD2
using Revise
using EvoId,LightGraphs,UnPack
using LaTeXStrings,KernelDensity
# using IDEvol
using PyPlot,PyCall
cm_eth = ColorMap([c for c in eth_grad_std.colors]);
include("../../code/graphs_utils/src/graphs_utils.jl")
include("../format.jl")

@load "../../code/simulations/transients/s_setting1.jld2" s_setting1
@load "../../code/simulations/transients/s_setting2.jld2" s_setting2


g = SimpleGraph(Edge.([(1,2),(2,4),(2,5),(2,3),(3,6),(3,7)]))
g = to_nx(g)
pos = nx.spring_layout(g,seed=3)
fig, ax = subplots(2,1,figsize = FIGSIZE_M)
ax_neutr = ax[1]; ax_adapt = ax[2]

node_size = 400
indiv_s = 5
winsx = 0.5; winsy = 0.8
xins = [0.3, -0.3];
yins = [-0.9,0.2];
xinsp = collect(-1.:0.05:1.);yinsp = collect(-1.:0.05:1.2) # range of inset plots

t = collect(1:7)
indposx = @. 0.05 * cos(2 * π / 7 * t); push!(indposx,0.)
indposy = @. 0.05 / 0.5 * sin(2 * π / 7 * t); push!(indposy,0.)

## adapt
vertex = [2,3]
ys = []
ym = []
for v in vertex
    push!(ys,[])
    push!(ym,[])
    for a in s_setting2[end]
        if a[1][] == v
            push!(ys[end],a[2][])
            push!(ym[end],a[3][4])
        end
    end
end

labels = Dict(collect(0:6) .=> string.(collect(1:7)))
nx.draw_networkx_nodes(g,pos,nodelist = collect(0:2:6),
                edgecolors="tab:blue",
                node_size = node_size,
                node_color = "tab:blue",
                # horizontalalignment = "right",
                # verticalalignment = "baseline",
                # alpha = 0.
                # options,
                # with_labels = false,
                ax = ax_adapt,
                )
nx.draw_networkx_nodes(g,pos,nodelist = collect(1:2:6),
                edgecolors="tab:red",
                linewidths = 1.,
                node_size = node_size,
                node_color = "tab:red",
                # horizontalalignment = "right",
                # verticalalignment = "baseline",
                # alpha = 0.
                # options,
                # with_labels = false,
                ax = ax_adapt,
                )
nx.draw_networkx_edges(g, pos, alpha=0.5, width=4,ax=ax_adapt)

##
for (i,v) in enumerate(vertex)
    xys =hcat(ys[i],ym[i])
    ks = InterpKDE(kde(xys))
    YS = KernelDensity.pdf.(Ref(ks),xinsp,yinsp')
    YS[YS .<= 0.01] .= 0
    i == 2 ? (_left = 1) : (_left = 0)
    ax_adapt_ins = ax_adapt.inset_axes([xins[i], yins[i], winsx, winsy], transform = ax_adapt.transData)
    line = matplotlib.lines.Line2D([pos[v-1][1], xins[i] + _left * winsx ], [pos[v-1][2], yins[i]],zorder=-10,c = "black",lw = 1) # -1 because python indexing starts with 0
    ax_adapt.add_line(line)
    line = matplotlib.lines.Line2D([pos[v-1][1], xins[i] + _left * winsx ], [pos[v-1][2], yins[i] + winsy],zorder=-10,c = "black",lw = 1)
    ax_adapt.add_line(line)
    # xinsp = collect(-4:0.05:4)
    cbar = ax_adapt_ins.contour(yinsp,xinsp,YS,
                                cmap = cm_eth,
                                levels = 3,
                                alpha = .8,
                                linewidths = 1.,
                                )
    ax_adapt_ins.scatter(ym[i],ys[i],s = 0.2,zorder=100,c="black",alpha=1)
    # cbar = plt.colorbar(cbar, ax = ax_adapt_ins, aspect = 15, anchor = (1,1), shrink = 0.8)
    # cbar.set_ticks([])
    # cbar.set_label("density",fontsize = 8)
    # ax_neutr_ins.axis("off")
    ax_adapt_ins.tick_params(left=true, bottom=false, labelleft=true, labelbottom=false)
    ax_adapt_ins.set_xlabel(L"u_1",fontsize = 8)
    sl = ax_adapt_ins.set_ylabel(L"s",fontsize = 8)
    plt.setp(sl,backgroundcolor="white")
    ax_adapt_ins.xaxis.set_label_position("top")
    ax_adapt_ins.yaxis.set_label_position("right")
    # ax_adapt_ins.axhline(1.,c = "r", ds = "steps", alpha = 0.8)
    ax_adapt_ins.set_yticks([-0.5,0.5])
    ytl = ax_adapt_ins.set_yticklabels([L"\theta_\circ", L"\theta_\bullet"])
    plt.setp(ytl[1],backgroundcolor="white",color = "tab:blue")
    plt.setp(ytl[2],color = "tab:red")
end

# ax_adapt_ins.
ax_adapt.text(.05, 0.9,"Setting with \n heterogeneous \n selection" ,
    # fontsize=12,
    # fontweight="bold",
    va="bottom",
    ha="left",
    transform=ax_adapt.transAxes ,
    fontsize = 10
)
for v in values(pos)
    ax_adapt.scatter(v[1] .+ indposx,
                v[2] .+ indposy ,
                zorder = 10,
                s = indiv_s,
                c = "w",
                # cmap = ColorMap(eth_grad_std.colors)
                )
end

ax_adapt.axis("off")
ax_adapt.margins(0.10)

##neutr
vertex = [2,3]
ym1 = []
ym2 = []
for v in vertex
    push!(ym1,Float64[])
    push!(ym2,Float64[])
    for a in s_setting1[end]
        if a[1][] == v
            push!(ym1[end],a[2][3])
            push!(ym2[end],a[2][2])
        end
    end
end

nx.draw_networkx_nodes(g,pos,
                edgecolors="tab:grey",
                node_size = node_size,
                node_color = "w",
                linewidths = 3.,
                # horizontalalignment = "right",
                # verticalalignment = "baseline",
                # alpha = 0.
                # options,
                # with_labels = false,
                ax = ax_neutr,
                )
    nx.draw_networkx_edges(g, pos, alpha=0.5, width=4,ax=ax_neutr)

pos_higher = Dict(collect(keys(pos)) .=> collect(values(pos)) .+ Ref([-0.10,0.14]))
for v in values(pos)
    ax_neutr.scatter(v[1] .+ indposx,
                v[2] .+ indposy,
                zorder = 10,
                s = indiv_s,
                c = "tab:grey")
end

## adding distribution
xins = xins .- 0.015 # for some reason the spacing is not equivalent between adapt and neutr figure
for (i,v) in enumerate(vertex)

    xys =hcat(ym1[i],ym2[i])
    ks = InterpKDE(kde(xys))
    xinsp = collect(-1.7:0.05:1.5)
    YS = KernelDensity.pdf.(Ref(ks),xinsp,xinsp')
    YS[YS .<= 0.01] .= 0
    i == 2 ? (_left = 1) : (_left = 0)
    ax_neutr_ins = ax_neutr.inset_axes([xins[i], yins[i], winsx, winsy], transform = ax_neutr.transData)
    line = matplotlib.lines.Line2D([pos[v-1][1], xins[i] + _left * winsx ], [pos[v-1][2], yins[i]],zorder=-10,c = "black",lw = 1) # -1 because python indexing starts with 0
    ax_neutr.add_line(line)
    line = matplotlib.lines.Line2D([pos[v-1][1], xins[i] + _left * winsx ], [pos[v-1][2], yins[i] + winsy],zorder=-10,c = "black",lw = 1)
    ax_neutr.add_line(line)
    # YS = pdf.(Ref(ks),yinsp)
    # ax_neutr_ins.plot(yinsp,YS)
    # ax_neutr_ins.fill_between(yinsp, YS)
    ax_neutr_ins.contour(xinsp,xinsp,YS,
                                cmap = cm_eth,
                                levels = 10,
                                alpha = .8,
                                linewidths = 1.,
                                )
    ax_neutr_ins.scatter(ym2[i],ym1[i],s = 0.2,zorder=100,c="black",alpha=1)

    ax_neutr_ins.yaxis.set_label_position("right")
    ax_neutr_ins.xaxis.set_label_position("top")
    # ax_neutr_ins.axis("off")
    ax_neutr_ins.tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
    ax_neutr_ins.set_xlabel(L"u_1",fontsize = 8)
    ul = ax_neutr_ins.set_ylabel(L"u_2",fontsize = 8)
    plt.setp(ul,backgroundcolor="white")
end
# ax_neutr_ins.set_box_aspect(1)
ax_adapt.text(.05, 1.0,"Setting with \n no selection" ,
    # fontsize=12,
    # fontweight="bold",
    va="bottom",
    ha="left",
    transform=ax_neutr.transAxes ,
    fontsize = 10
)

ax_neutr.axis("off")
ax_neutr.margins(0.10)

_let = [L"\textbf{a}",L"\textbf{b}", L"\textbf{c}", L"\textbf{d}"]
for (i,ax) in enumerate([ax_neutr,ax_adapt])
    ax.text(0., 1.05, _let[i],
        fontsize=12,
        fontweight="bold",
        va="bottom",
        ha="left",
        transform=ax.transAxes ,
    )
    # ax.set_xticks(-1.:0.5:1.)
end

fig.set_facecolor("None")
# fig.tight_layout()

fig.savefig("conceptual_v3.pdf",
            dpi=1200,
            bbox_inches = "tight")
display(fig)
