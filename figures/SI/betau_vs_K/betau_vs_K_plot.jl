# Adaptive landscape with 9 patches along a line
# /!\ t
using Random
cd(@__DIR__)
using Dates
using JLD2
using Revise
using EvoId,LightGraphs,UnPack
using LaTeXStrings, Printf
using IDEvol
using PyPlot
include("../../format.jl")
Line2D = matplotlib.lines.Line2D

@load "betau_vs_K.jld2" df

sort!(df,:K1)
###############
###plotting####
###############
df[!,"Q_ST_u"] = similar(df.betau)
df[!,"Q_ST_s"] = similar(df.betas)

[r.Q_ST_u = r.betau ./ r.gammau for r in eachrow(df)]
[r.Q_ST_s = r.betas ./ r.gammas for r in eachrow(df)]

plt.clf()
fig, axs = plt.subplots(1,3,figsize=(FIGSIZE_M[1]*2,FIGSIZE_M[2]*0.8))

# first row corresponds to K=10 for which there is too much stochasticity
axs[1].boxplot(df.Q_ST_u)
axs[1].set_ylabel( L"Q_{ST,u}")
axs[2].boxplot(df.Q_ST_s)
axs[2].set_ylabel( L"Q_{ST,s}");
axs[3].boxplot(df.N)
axs[3].set_ylabel( L"N");
[ax.set_xlabel("K") for ax in axs]

for ax in axs
    ax.set_xticks(1:(size(df,1)));
    ax.set_xticklabels(df.K1 .|> Int)
end
# custom legend
# fig.legend(handles=[Line2D([0], [0], color="grey", linestyle=":", label="star graph"),
#             Line2D([0], [0], color="grey", label="complete graph")],loc="upper left") #colors
# fig.legend(handles=[Line2D([0], [0], color="tab:blue", linestyle=":", label="m = 0.1"),
#                     Line2D([0], [0], color="tab:red", label="m=0.3")],loc="upper right") #linestyle
# other legend style
# fig.legend(bbox_to_anchor=[1.0,0.5], bbox_transform = axs[3].transAxes)

_let = ["a","b","c"]
_offx = - 0.2
_offy =  1.05
for (i,ax) in enumerate(axs)
    ax.text(_offx, _offy, _let[i],
        fontsize=12,
        fontweight="bold",
        va="bottom",
        ha="left",
        transform=ax.transAxes ,
    )
end

fig.tight_layout()

fig.savefig("betau_vs_K.pdf",
            dpi=1200,
            bbox_inches = "tight",
            )

gcf()
