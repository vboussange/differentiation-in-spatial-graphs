using Base: between
cd(@__DIR__)
using Pkg; Pkg.activate(".")
# using PyPlot
using ArchGDAL; const AG = ArchGDAL
using Statistics
using ProgressMeter
using Glob, JLD2
using Dates
verbose = true
test = false
plotting = false
simulation = true

include("src/extract_graph_from_raster.jl")
savedir = "results/lvl2_frac_1km_ver004"
isdir(savedir) ? nothing : mkpath(savedir)


if simulation
    datalist = glob("./data/lvl2_frac_1km_ver004/*.tif")
    println("Starting computation")
    println(now())
    Threads.@sync Threads.@threads for f in datalist
        dataset = AG.read(f)
        fsplit = split(f,"_")
        window_size = 100#km
        habitat = parse(UInt16,fsplit[9]) #lvl1 forest
        raster_i = 1:window_size:AG.width(dataset)
        raster_j = 1:window_size:AG.height(dataset)

        raster = Dict(
                        "mean_betw_centrality" => zeros(length(raster_i), length(raster_j)),
                        "var_betw_centrality" => zeros(length(raster_i), length(raster_j)),
                        "mean_degree" => zeros(length(raster_i), length(raster_j)),
                        )
        raster_g = Array{SimpleGraph{Int16}}(undef,length(raster_i), length(raster_j))

        # preallocating
        data = AG.read(dataset, 1, 1, 1, window_size, window_size) #x, y , window size

        # main loop
        for (ii,i) in enumerate(raster_i[1:end-1]), (jj,j) in enumerate(raster_j[1:end-1])
            data .= AG.read(dataset, 1, i, j, window_size, window_size) #x, y , window size
            ncells = count(data .> 0)
            if ncells > 0
                verbose && println(ncells, " cells of habitats ", habitat, " were found")
                g = extract_graph_1km(data)

                # calculating metrics
                bet = betweenness_centrality(g)
                raster["mean_betw_centrality"][ii,jj] = mean(bet)
                raster["var_betw_centrality"][ii,jj] = var(bet)
                raster["mean_degree"][ii,jj] = mean(degree(g))
                # storing graph
                raster_g[ii,jj] = g

                if verbose
                    println("mean_betw_centrality = ", raster["mean_betw_centrality"][ii,jj])
                    println("var_betw_centrality = ", raster["var_betw_centrality"][ii,jj])
                    println("mean_degree = ", raster["mean_degree"][ii,jj])
                end
            end

            if test
                if i > 1
                    break
                end
            end
        end
        savename = split(f,"/")[end]
        savename = string(split(savename,".")[1],".jld2")
        @save joinpath(savedir,savename) raster raster_g
    end
    println(now())
    println("Computation over")
end

# plotting
if plotting
    datalist = glob(savedir*"/*.jld2")
    for _dat in datalist
        @load _dat raster
        savename = split(_dat,"/")[end]
        savenamfig = string(split(savename,".")[1],".pdf")
        fig, ax = plt.subplots();
        ax.imshow(raster')
        fig.savefig(joinpath(savedir,savenamfig))
        plt.close(fig)
    end
end
