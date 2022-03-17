cd(@__DIR__)
using Pkg; Pkg.activate(".")
using ArchGDAL; const AG = ArchGDAL
using Statistics
using ProgressMeter
using Glob, JLD2
using Dates
verbose = true
test = false
plotting = false
simulation = true
area_threshold = 500
bio = "bio12" #bio12 for precipitations, bio1 for temperature
habitat = "forest"

_savename = "$(habitat)_lvl1_sqrtk_rtheta_$(bio)_nearest"
include("src/extract_graph_from_raster.jl")
include("src/graph_metrics.jl")

savedir = "results/$habitat"
isdir(savedir) ? nothing : mkpath(savedir)


@time if simulation
    dataset_hab = AG.read("./data/lvl1_frac_1km_ver004/iucn_habitatclassification_fraction_lvl1__100_Forest__ver004.tif")
    dataset_temp = AG.read("./data/chelsa/CHELSA_$(bio)_reprojected_nearest.tif")

    window_size = 100#km
    raster_i = 1:window_size:AG.width(dataset_hab)
    raster_j = 1:window_size:AG.height(dataset_hab)

    raster = Dict(
                    "sqrtk" => zeros(length(raster_i), length(raster_j)),
                    "assortativity" => zeros(length(raster_i), length(raster_j)),
                    )
    raster_g = Array{SimpleGraph{Int16}}(undef,length(raster_i), length(raster_j))

    # preallocating
    data_hab = AG.read(dataset_hab, 1) #x, y , window size
    data_temp = AG.read(dataset_temp, 1)

    # main loop
    Threads.@threads for ii in 1:length(raster_i)-1
        i = raster_i[ii]
        for (jj,j) in enumerate(raster_j[1:end-1])
            datahab_ij = view(data_hab, i:i+window_size,j:j+window_size)
            ncells = count(datahab_ij .> 0)
            if ncells > 0
                datatemp_ij = view(data_temp, i:i+window_size,j:j+window_size)
                verbose && println(ncells, " cells of habitats ", habitat, " were found")
                g, B = extract_graph_1km(datahab_ij, area_threshold)
                # calculating metrics
                raster["sqrtk"][ii,jj] = sqrtk(g)
                raster["assortativity"][ii,jj] = assortativity(g, datatemp_ij[B])
                # storing graph
                raster_g[ii,jj] = g

                if verbose
                    println("sqrtk = ", raster["sqrtk"][ii,jj])
                    println("assortativity = ", raster["assortativity"][ii,jj])
                end
            end
        end

        if test
            if i > 1
                break
            end
        end
    end
    _savename = string(_savename,".jld2")
    @save joinpath(savedir,_savename) raster raster_g
end
println(now())
println("Computation over")

if plotting
    using PyPlot
    datalist = glob(savedir*"/*.jld2")
    for _dat in datalist
        @load _dat raster
        savename = split(_dat,"/")[end]
        for k in keys(raster)
            savenamfig = string(split(savename,".")[1],k,"area",area_threshold,".pdf")
            fig, ax = plt.subplots();
            im = ax.imshow(raster[k]')
            plt.colorbar(im, label=k)
            fig.savefig(joinpath(savedir,savenamfig))
            plt.close(fig)
        end
    end
end
