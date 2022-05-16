using LightGraphs
using PyPlot
using Printf
cd(@__DIR__)
include("../../../code/graphs_utils/src/graphs_utils.jl")
include("../../format.jl")

fig, axx = plt.subplots(2,3,figsize=(FIGSIZE_L[1],FIGSIZE_L[2]*0.6))
_let = [L"\textbf{a}",L"\textbf{b}", L"\textbf{c}", L"\textbf{d}", L"\textbf{e}", L"\textbf{f}"]
gxs = vcat( [LightGraphs.grid([2,2,2]) for _ in 1:3],# seed = 2
    [LightGraphs.grid([8,1],periodic=true) for _ in 1:2], # seed = 2
    [SimpleGraph(Edge.([(1,2), (2,3),(2,4),(4,5),(5,8),(5,6),(6,7),(7,8),(5,7),(6,8)]))],
    # SimpleGraph(Edge.([(1,2), (2,3),(2,4),(4,5),(5,8),(5,6),(6,7),(7,8),(5,7),(6,8)]))]
    )
# pos = nx.nx_agraph.graphviz_layout(g, prog="neato")
# pos = nx.planar_layout(gx)
poss = vcat([nx.spring_layout(to_nx(gx),seed=2) for gx in gxs[1:5]],
    Dict(0 => [-2,-0.3], 1 => [-1.5,0], 2=> [-2,0.3], 3=> [-0.8,0],4 => [0.8,0], 6=> [1.5,0],7=>[2,-0.3],5 => [2,0.3]))
θs = [[ -1, 1, -1 , 1, -1 , 1, -1 , 1],
    [ -1, 1, 1 , -1, 1 , -1, -1 , 1],
    [ -1, -1, 1 , -1, 1 , -1, 1 , 1],
    [ -1, -1, 1 , 1, -1 , -1, 1 , 1],
    [ -1, 1, -1 , 1, -1 , 1, -1 , 1],
    [ 1, 1, 1 , 1, -1 , -1, -1 , -1]]
for i in 1:length(θs)
    ax = axx[i]
    g = gxs[i]
    gx = to_nx(gxs[i])
    θ = θs[i]
    pos = poss[i]
    soptim = Dict(collect(1:length(θ)) .=> θ)
    for h in [-1,1]
        _c = h == -1 ? "tab:blue" : "tab:red"
        node_list = [k-1 for (k,v) in enumerate(θ) if v == h]
        nx.draw_networkx_nodes(gx,pos,
                        nodelist = node_list,
                        edgecolors = _c,
                        node_size = 100.,
                        node_color = _c,
                        # linewidths = 3.,
                        # horizontalalignment = "right",
                        # verticalalignment = "baseline",
                        # alpha = 0.
                        # options,
                        # with_labels = false,
                        ax = ax,
                        )
    end
    # nx.draw_networkx_labels(gx,pos,ax=ax)
    nx.draw_networkx_edges(gx, pos, alpha=0.5, width=4,ax=ax)
    ax.margins(0.11)
    ax.axis("off")
    ax.set_title(L"r_\Theta = "*(@sprintf "%1.2f" soptim_correlation(g,soptim)),
                y = -0.1)
    ax.text(0, 1, labels[i], va="bottom", ha="left",
        fontsize=12,
        transform=ax.transAxes ,
        )
end
fig.tight_layout()
gcf()
fig.savefig("r_theta_for_different_graphs.pdf",
            dpi=1200,
            bbox_inches = "tight",
            )
