using Random
cd(@__DIR__)
using Dates
using JLD2
using Revise
using EvoId,LightGraphs,UnPack
using LaTeXStrings,KernelDensity
using IDEvol
using PyPlot,PyCall
include("code/graphs_utils/src/graphs_utils.jl")
include("../format.jl")
@load "code/simulations/setting_1/transients/concept_figure_sim_neutr.jld2" s_neutr

@load "code/simulations/setting_2/transients/concept_figure_sim_adapt.jld2" s_adapt


g = SimpleGraph(Edge.([(1,2),(2,4),(2,5),(2,3),(3,6),(3,7)]))
g = to_nx(g)
pos = nx.spring_layout(g,seed=3)
fig, ax_adapt = subplots(figsize = FIGSIZE_M)

node_size = 800
indiv_s = 5
winsx = 0.8; winsy = 0.8
xins = [0.3, -0.6];
yins = [-0.9,0.4];
xinsp = collect(-2.:0.05:2.);yinsp = collect(-0.8:0.05:0.8) # range of inset plots

t = collect(1:10)
indposx = @. 0.1 * cos(2 * π / 10 * t)
indposy = @. 0.1 * sin(2 * π / 10 * t)

## adapt
vertex = [2,3]
ys = []
ym = []
for v in vertex
    push!(ys,[])
    push!(ym,[])
    for a in s_adapt[end]
        if a[1] == v
            push!(ys[end],a[2])
            push!(ym[end],a[3][5])
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
                                cmap = ColorMap(eth_grad_std.colors),
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
    ax_adapt_ins.set_yticks([-1,1])
    ytl = ax_adapt_ins.set_yticklabels([L"\theta_\bullet", L"\theta_\bullet"])
    plt.setp(ytl[1],backgroundcolor="white",color = "tab:blue")
    plt.setp(ytl[2],color = "tab:red")
end


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


fig.set_facecolor("None")

fig.savefig("conceptual_onlyadapt.png",
            dpi=1200,
            bbox_inches = "tight")
