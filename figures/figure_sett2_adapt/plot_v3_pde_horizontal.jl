cd(@__DIR__)
using LightGraphs
using DataFrames
using CSV
using Statistics
using JLD2
using EvoId
using Printf
using Polynomials
using KernelDensity
using ColorSchemes
using PyPlot
cm_eth = ColorMap([c for c in eth_grad_std.colors]);
include("../../code/graphs_utils/src/graphs_utils.jl")
include("../format.jl")

@load "../../code/simulations/setting_2/explo_r&m-pde/pde_meanfield_data.jld2" df_explo # data for heatmap
@load "../../code/analytics/adaptive_dynamics_PDE_data.jld2" dfg S # data for distribution

# processing data for heatmap
x = unique(df_explo.m)
y = unique(df_explo.rθ)
Q_ST_s = zeros(length(x),length(y))
pop = zeros(length(x),length(y))

for i in 1:length(x)
    for j in 1:length(y)
        k = findfirst(r -> (r.m == x[i]) && (r.rθ == y[j]),eachrow(df_explo) )
        Q_ST_s[i,j] = df_explo.Q_ST_s[k]
        # Q_ST_s[i,j] = df_explo.γ[k]
        pop[i,j] = df_explo.npop[k]
    end
end
df_explo_g = groupby(df_explo, :rθ, sort = true)
rt = sort!(unique(df_explo.rθ))
pde_data = Dict("m"=>x,"rθ"=>y,"βs"=>Q_ST_s,"N"=>pop)

# function to obtain arg of local maximum
function get_types(uend)
    du1 = uend[1,2:end]-uend[1,1:end-1]
    du10 = ((du1[2:end] .* du1[1:end-1]) .< 0)
    c = count(du10)
    println("derivative cancels at ", S[2:end-1][du10])
    if c == 3
        # the middle point where derivative is 0 corresponds to a local minimum
        return S[2:end-1][du10][[1,3]]
    elseif c == 1
        return S[2:end-1][du10]
    else 
        println("Cannot handle uend where derivative cancels $c times")
        return S[2:end-1][du10]
    end
end


#################
### PDE RESULT###
#################
fig = plt.figure(
            constrained_layout=true,
            figsize = (FIGSIZE_S[1]*2,FIGSIZE_S[2]*1.2)
            )
gs = fig.add_gridspec(2, 2,)

ax3 = fig.add_subplot(py"$(gs)[0:2,1]")
ax1 = fig.add_subplot(py"$(gs)[0,0]")
ax2 = fig.add_subplot(py"$(gs)[1,0]")

cbbox1 = mplt.inset_axes(ax3,width = "20%", height = "32%",loc=4)
for cbox in [cbbox1]
    cbox.tick_params(axis="both", left=false, top=false, right=false, bottom=false, labelleft=false, labeltop=false, labelright=false, labelbottom=false)
    # [cbox.spines[k].set_visible(false) for k in keys(cbox.spines)]
    # cbox.axis("off")
    cbox.set_facecolor([1,1,1,0.7])
end
cbar1 = mplt.inset_axes(cbbox1, "8%", "50%", loc = 7)

pc = ax3.pcolormesh( pde_data["rθ"],pde_data["m"], pde_data["βs"],cmap = cm_eth,
    # extent=[1,tend,S[1],S[end]]
    # vmax = K1
    # norm = matplotlib.colors.LogNorm()
    )
_cb = fig.colorbar(pc,cax =cbar1)
_cb.ax.set_xlabel(L"Q_{ST,s}",labelpad = 5, loc="right")
_cb.ax.xaxis.set_label_position("top")
# _cb.ax.tick_params(axis="x",direction = "in", labeltop = true)
# cbar1.xaxis.set_label_position("top")
cbar1.yaxis.set_ticks_position("left")
# _cb.outline.set_edgecolor("white")
# _cb.outline.set_linewidth(2)

ax3.set_ylabel(L"m")
# ax3.set_xlabel(L"r_\theta")

gcf()

###############
## drawing m* 
###############

rθ = pde_data["rθ"]
p = 1; θ=0.5
mstar = 1. ./ (1 .- rθ ) .* 4 * p * θ^2/(1 + 3*p * θ^2) #critical threshold
ax3.plot(rθ,mstar,c = "red",label = L"m^*")
ax3.legend(bbox_to_anchor=(0.66, 1.), bbox_transform=ax3.transAxes)
ax3.set_ylim(0.,1.)
ax3.set_xlabel(L"r_\theta")


##############################
##### plotting distribution ##
##############################
## ax1
ax = [ax1,ax2]
[ax[i].set_ylim([0.,4.]) for i in 1:2]

cols = ["tab:blue", "tab:orange", "tab:red"]
for (i,df) in enumerate(dfg)
    rθ = df.rθ[1]
    p = 1; θ=0.5
    mstar = 1. ./ (1 .- rθ ) .* 4 * p * θ^2/(1 + 3*p * θ^2) #critical threshold
    for (j,r) in enumerate(eachrow(df)[[1,7,10]])
        # plotting trait distribution
        ax[i].fill_between(S,
                        r.uend[1,:], 
                        alpha=0.4,
                        color = cols[j],
                        label = (i==2 ? "m = "*@sprintf("%.1f",r.m) : nothing),
                        )
        ax[i].plot(S,
                    r.uend[1,:], 
                    color = cols[j])
        s = get_types(r.uend)
        # plotting maximum corresponding to adaptive dynamics prediction
        # for some reason when curves are unimodal they return whether 1 or 5 derivatives which cancel out.
        # only the first one is meaningful
        if length(s) == 1 || length(s) > 2
            ax[i].vlines(s[1], ax[i].get_ylim()..., linestyles = "dashed", color = cols[j], linewidth=2)
        elseif length(s) == 2
            #calculating mean trait :
            _pop = (r.uend[1,s[1] .== S][] + r.uend[1,s[2] .== S][])
            _meantrait = (r.uend[1, s[1] .== S][] * s[1] + r.uend[1,s[2] .== S][] * s[2]) / _pop
            ax[i].vlines(_meantrait, ax[i].get_ylim()..., linestyles = "dashed", color = cols[j], linewidth=2)
        end
    end
    ax[i].set_title(
            latexstring("r_\\theta = ",(df.rθ[1]), ", m^\\star = ", @sprintf("%.1f",mstar) ),
            fontsize=12,
            fontweight="bold",
            )
    ax[i].set_ylabel(latexstring("\\bar{n}^{\\bcirc}"))
end
ax[2].set_xlabel("Adaptive trait")
ax[2].legend()

_let = [L"\textbf{a}",L"\textbf{b}", L"\textbf{c}", L"\textbf{d}"]
for (i,ax) in enumerate([ax1,ax2,ax3])
    _x = i < 3 ? -1.4 : -0.2
    _y = i == 2 ? 0.45 : 1.1
    ax.text(_x, 
        _y,
        _let[i],
        fontsize=12,
        va="bottom",
        ha="left",
        transform=ax3.transAxes,
    )
    ax.set_xticks(-1.:0.5:1.)
end


fig.tight_layout()

gcf()

fig.savefig("pde_v3_horizontal.pdf",
            dpi=1200,
            # bbox_inches = "tight",
            )
gcf()