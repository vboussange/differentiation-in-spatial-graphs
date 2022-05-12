#=
This script plots the last figure of the main text

habitat considered is forest

Results of simulations are plotted in contrast to the other script

v2 : bar plot
=#

# /!\ this script uses PyCall, that relies on python
# Make sure to install Python, as well as the packages
# "cartopy", "rasterio", "networkx", "matplotlib", "shapely".
# Then uncomment the two following lines, making sure to direct
# ENV["PYTHON"] your your python installation (see https://github.com/JuliaPy/PyCall.jl)
# for more explanations.
# Those two lines can be uncommented after the first usepackage

# ENV["PYTHON"] = "/usr/local/anaconda3/envs/land2graph/bin/python"
# using Pkg; Pkg.build("PyCall")

using PyCall
using Printf
using JLD2
ccrs = pyimport("cartopy.crs")
cf = pyimport("cartopy.feature")
rasterio = pyimport("rasterio")
using PyPlot
include("../../format.jl")
geometry = pyimport("shapely.geometry")
include("../../../code/graphs_utils/src/graphs_utils.jl")
cd(@__DIR__)
# this is needed for PyCall v1.92.5 
# cm_eth = ColorMap([c for c in eth_grad_std.colors])
# cmap = cm_eth
using ColorSchemes
cmap = ColorMap(ColorSchemes.tableau_red_blue.colors)
rivers_50m = cf.NaturalEarthFeature("physical", "rivers_lake_centerlines", "50m")

idx_m_toplot = [1,4] #idx to display results

#################################
##### Reading habitat raster ####
#################################
filename_temp = "../../../code/graphs_utils/Land2Graph/data/CHELSA_bio1_reprojected_nearest.tif"
filename_hab = "../../../code/graphs_utils/Land2Graph/data/iucn_habitatclassification_fraction_lvl1__100_Forest__ver004.tif"
# loading rasters for (A) and (B)
data_src_temp = rasterio.open(filename_temp, "r")
data_src_hab = rasterio.open(filename_hab, "r")

xmin = 90.78596; xmax = 105.8718; ymin = 22.91373; ymax = 34.47213
# pixel coordinate for longitude / latitudes for top left corner
row, col = data_src_hab.index(xmin, ymax)
row2, col2 = data_src_hab.index(xmax, ymin)
# window_size = 10
area_threshold = 500

mytemp = data_src_temp.read(1, window = rasterio.windows.Window( col, row, col2-col, row2-row)) /10. .- 273.15 .|> Float64
myhab =  data_src_hab.read(1, window = rasterio.windows.Window( col, row, col2-col, row2-row)) .|> Float64

real_graph_df = load("../../../code/graphs_utils/real_graphs/realistic_graph_hengduans.jld2", "sub_df_interesting_graphs")

##########################################################
# small window within the considered region for fig. (A) #
##########################################################

# create figure
# define cartopy crs for the raster
crs = ccrs.PlateCarree()

# fig,axs = plt.subplots(1, 2,            
#                     figsize = FIGSIZE_L,
#                     subplot_kw = Dict("projection"=>crs))
# ax1,ax2 = axs

# create figure
fig = plt.figure(
            # constrained_layout=true,
            figsize = (FIGSIZE_L[1],FIGSIZE_L[1])
            )
gs = fig.add_gridspec(2, 2)
gs.update(wspace=0.3) # set the spacing between axes.
ax1 = fig.add_subplot(py"$gs[0, 1]", projection=crs)
ax2 = fig.add_subplot(py"$gs[0, 0]", projection=crs)
axa = fig.add_subplot(py"$gs[1, 0]")
axb = fig.add_subplot(py"$gs[1, 1]")

################
#### ax1 #######
################
# ax1.set_title("Hengduan region")
ax1.set_xmargin(0.05)
ax1.set_ymargin(0.10)


################
#### (B) #######
################
# plot raster
pos1 = ax1.imshow(
                myhab/10,
                interpolation="nearest",
                cmap = :Purples,
                vmin = 0.0,
                vmax = 100.0,
                extent=[xmin, xmax, ymin, ymax],
                transform=crs,
                )
gl = ax1.gridlines(draw_labels=true)
gl.xlabels_top = false
gl.ylabels_right = false
gl.xlines = true
gl.ylines = true
gl.rotate_labels = false
gcf()
# plot features
ax1.coastlines(resolution="auto", color="red")
ax1.add_feature(cf.BORDERS, linestyle = "--", edgecolor = "grey")
ax1.add_feature(rivers_50m, facecolor="None", edgecolor="deepskyblue")
gcf()


# adding crosses on areas of the graph selected
# we add a cross on central pixel of the sub window
for r in eachrow(real_graph_df)
    # coordinates in the array
    y = col + r.geographic_coord[2]*r.window_size - (r.window_size/2-1)
    x = row + r.geographic_coord[1]*r.window_size - (r.window_size/2-1)
    # latitude / longitude
    lat, long = data_src_hab.xy(x, y)
    ax1.scatter(lat, long, c = "r", marker="+", transform=crs, zorder = 10)
end
ax1.axes.xaxis.set_visible(false)
ax1.axes.yaxis.set_visible(false)
gcf()
############################
############ (A) ###########
############################
ax2.set_xmargin(-0.01)
ax2.set_ymargin(-0.01)
[i.set_linewidth(3) for i in ax2.spines.values()]
[i.set_color("r") for i in ax2.spines.values()]

# plotting sub raster 
idx = 25 # choosing first graph considered {21,22,35}
wc = real_graph_df[idx,"geographic_coord"] # window coordinate (top left)
window_size = real_graph_df[idx,"window_size"]
xmin, ymax = data_src_hab.xy(row + wc[1]*window_size - (window_size-1), col + wc[2]*window_size - (window_size-1))
xmax, ymin = data_src_hab.xy(row + wc[1]*window_size + 1, col + wc[2]*window_size + 1) 
# notice that extent is given by coordinates of 0 and window_size + 1 pixel
pos2 = ax2.imshow(myhab[(wc[1]*window_size - (window_size-1)):(wc[1]*window_size), (wc[2]*window_size - (window_size-1)):(wc[2]*window_size)]/10,
                cmap = :Purples,
                origin="upper",
                vmin = 0.0,
                vmax = 100.0,
                extent=[xmin, xmax, ymin, ymax],
                transform=crs,
                )
gl = ax2.gridlines(draw_labels=true)
gl.xlabels_top = false
gl.ylabels_right = false
gl.xlines = false
gl.ylines = false
gl.rotate_labels = false
gcf()
# cb = colorbar(pos2, ax=ax2, label = "Grassland coverage", shrink=0.5)
# vals = cb.ax.get_yticks()
# cb.ax.set_yticklabels([@sprintf("%15.0f", x)*L"\%" for x in vals])

# Plotting graph on top
gx = to_nx(real_graph_df[idx, "graphs"])
coord_graph = real_graph_df[idx,"pos_interesting_graph_unorganised"]
temp = real_graph_df[idx, "temps"]
# assiging coordinates to the graph nodes
pos = Dict{Int,Array}()
row2 = row + wc[1]*window_size - (window_size-1)
col2 = col + wc[2]*window_size - (window_size-1)
for i in keys(coord_graph)
    xs, ys = data_src_hab.xy(row2 + coord_graph[i][2] + 0.5, col2 + coord_graph[i][1] + 0.5, offset="center")
    pos[i] = [xs[], ys[]]
end
nx.draw_networkx_nodes(gx,
                        pos,
                        node_size = 25.,
                        node_color = temp,
                        cmap=cmap,
                        ax = ax2,
                        )
nx.draw_networkx_edges(gx, pos, alpha=0.9, width=1,ax=ax2)
gcf()



# colorbar grassland
# fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.96, 0.745, 0.02, 0.1])
cb = fig.colorbar(pos1, cax=cbar_ax, label = "Grassland\ncoverage",shrink=10)
vals = cb.ax.get_yticks()
cb.ax.set_yticklabels([@sprintf("%15.0f", x)*L"\%" for x in vals])

## Colorbar vertices
# colorbar grassland
cbar_ax2 = fig.add_axes([0.96, 0.585, 0.02, 0.1])
sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=matplotlib.colors.Normalize(vmin=minimum(temp), vmax=maximum(temp)))
cb = fig.colorbar(sm, cax=cbar_ax2, label = "Standardized\ntemperature")
# vals = cb.ax.get_yticks()
# cb.ax.set_yticklabels([@sprintf("%15.0f", x)*L"\%" for x in vals])


## draw nice lines between (B) and (A)
x = row + wc[1]*window_size - (window_size-1)
y = col + wc[2]*window_size - (window_size-1)
xyB = data_src_hab.xy(x, y)
xyA = (1, 1)
con = matplotlib.patches.ConnectionPatch(
    xyA=xyA,
    xyB=xyB,
    coordsA=ax2.transAxes,
    coordsB="data",
    axesA=ax2,
    axesB=ax1,
    color = "r")
ax2.add_artist(con)
# low
xyA = (1, 0)
con = matplotlib.patches.ConnectionPatch(
    xyA=xyA,
    xyB=xyB,
    coordsA=ax2.transAxes,
    coordsB="data",
    axesA=ax2,
    axesB=ax1,
    color = "r")
ax2.add_artist(con)


###################################
### plotting simulation results ###
###################################
using GLM, Interpolations, Polynomials
using DataFrames

date_sim = "2022-01-15"
@load "../../../code/simulations/setting_2/realistic_graphs/setting_2_mu_01_realistic_graphs/setting_2_mu_01_realistic_graphs_$(date_sim)_aggreg.jld2" df_aggreg
df_aggreg_g = groupby(df_aggreg,:m,sort=true)
totsize = size(df_aggreg_g[1],1)

# calculating GLM for each m and extracting coefficients for Q_ST_u
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

axa.axes.xaxis.set_visible(true)
axa.axes.yaxis.set_visible(true)
axa.axes.yaxis.set_visible(true)
axb.axes.yaxis.set_visible(true)

_let = [L"\textbf{a}",L"\textbf{b}", L"\textbf{c}", L"\textbf{d}"]
for (i,ax) in enumerate([ax2,ax1, axa, axb])
    _x = -0.2
    _y = ax == ax1 ? 1.22 : 1.05
    ax.text(_x, _y, _let[i],
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
axa.set_facecolor("None")
axb.set_facecolor("None")
fig.set_facecolor("None")
###################
### annotating ####
##################
fig.tight_layout()


fig.savefig("land2graphs_with_simu.pdf",
            dpi=1200,
            bbox_inches = "tight",
            )

gcf()