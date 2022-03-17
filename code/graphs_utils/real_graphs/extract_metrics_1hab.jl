#=
This script contains utilities to 
extract metrics from a spatial window of habitats.

utilities are used in `extract_graphs_for_simu_hengduans_v2.jl``

=#

using ProgressMeter
function extract_metrics_1hab(myim, mytemp, myprec, habitat, myrich, window_size, area_threshold; radius = 1.5)
    println("Computing for ", habitat, " window_size = " , window_size, " area_threshold = ", area_threshold)
    # dividing the raster in coarser raster
    raster_i = 1:window_size:size(myim,1)
    raster_j = 1:window_size:size(myim,2)
    raster_sqrtk = fill(NaN, length(raster_i), length(raster_j))
    # raster_cl = fill(NaN, length(raster_i), length(raster_j))
    raster_assortativity_temp = fill(zeros(1,1), length(raster_i), length(raster_j))
    # raster_assortativity_prec = fill(NaN, length(raster_i), length(raster_j))
    raster_assortativity_prec = fill(zeros(1,1), length(raster_i), length(raster_j))
    raster_rich = fill(NaN, length(raster_i), length(raster_j))
    raster_temp = fill(NaN, length(raster_i), length(raster_j))
    raster_prec = fill(NaN, length(raster_i), length(raster_j))
    raster_g = fill(SimpleGraph(Int16(0)), length(raster_i), length(raster_j))
    raster_B = fill(BitArray(undef,1,1), length(raster_i), length(raster_j))

    # looping over the windows of the coarser raster
    for (ii,i) in enumerate(raster_i)
        for (jj,j) in enumerate(raster_j)
            datahab_ij = view(myim, i:i+window_size - 1, j:j+window_size - 1)
            ncells = count(datahab_ij .> 0) # is there the considered habitat in the window?
            if ncells > 0
                datatemp_ij = view(mytemp, i:i+window_size - 1,j:j+window_size -1 ) # needed for assortativity
                dataprec_ij = view(myprec, i:i+window_size - 1,j:j+window_size -1 ) # needed for assortativity
                myrich_ij = view(myrich, i:i+window_size - 1, j:j+window_size - 1 )
                verbose && println(ncells, " cells of habitats ", habitat, " were found")
                g, B = extract_graph_1km(datahab_ij, area = area_threshold, radius = radius)
                raster_g[ii,jj] = g
                raster_B[ii,jj] = B
                # calculating metrics
                raster_sqrtk[ii,jj] = sqrtk(g)
                # raster_cl[ii,jj] = (mean∘betweenness_centrality)(g)
                raster_assortativity_temp[ii,jj] = datatemp_ij
                raster_assortativity_prec[ii,jj] = dataprec_ij
                raster_rich[ii,jj] = mean(myrich_ij[B])
                raster_temp[ii,jj] = mean(datatemp_ij[B])
                raster_prec[ii,jj] = mean(dataprec_ij[B])
                if verbose
                    println("sqrtk = ", raster_sqrtk[ii,jj])
                    # println("cl = ", raster_cl[ii,jj])
                    println("assortativity temp = ", raster_assortativity_temp[ii,jj])
                    println("assortativity prec = ", raster_assortativity_prec[ii,jj])
                end
            end
        end
    end
    return raster_g,raster_B, raster_sqrtk, raster_assortativity_temp, raster_assortativity_prec, raster_rich, raster_temp, raster_prec
end
# data_src.close()

function extract_metrics_1hab(myim, mytemp, window_size, area_threshold; radius = 1.5)
    println("Computing for window_size = " , window_size, " area_threshold = ", area_threshold)
    # dividing the raster in coarser raster
    raster_i = 1:window_size:size(myim,1)
    raster_j = 1:window_size:size(myim,2)
    raster_sqrtk = fill(NaN, length(raster_i), length(raster_j))
    # raster_cl = fill(NaN, length(raster_i), length(raster_j))
    raster_assortativity_temp = fill(zeros(1,1), length(raster_i), length(raster_j))
    raster_temp = fill(NaN, length(raster_i), length(raster_j))
    raster_g = fill(SimpleGraph(Int16(0)), length(raster_i), length(raster_j))
    raster_B = fill(BitArray(undef,1,1), length(raster_i), length(raster_j))

    # looping over the windows of the coarser raster
    @showprogress for (ii,i) in enumerate(raster_i)
        for (jj,j) in enumerate(raster_j)
            datahab_ij = view(myim, i:i+window_size - 1, j:j+window_size - 1)
            ncells = count(datahab_ij .> 0) # is there the considered habitat in the window?
            if ncells > 0
                datatemp_ij = view(mytemp, i:i+window_size - 1,j:j+window_size -1 ) # needed for assortativity
                verbose && println(ncells, " cells of habitats were found")
                g, B = extract_graph_1km(datahab_ij, area = area_threshold, radius = radius)
                raster_g[ii,jj] = g
                raster_B[ii,jj] = B
                # calculating metrics
                raster_sqrtk[ii,jj] = sqrtk(g)
                raster_assortativity_temp[ii,jj] = datatemp_ij
                raster_temp[ii,jj] = mean(datatemp_ij[B])
            end
        end
    end
    return raster_g, raster_B, raster_sqrtk, raster_assortativity_temp, raster_temp
end
# data_src.close()