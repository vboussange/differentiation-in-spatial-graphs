cd(@__DIR__)
name_sim = "setting_2_mu_01_M=7_hetero_2_trait_comp_exp"
date_sim = "2022-03-17"
tend = 2000.
using EvoId,JLD2
using DataFrames
using Glob
using LightGraphs

df_arrays = []
flist = glob("$name_sim/*$date_sim*.jld2")
for f in flist
    try
        push!(df_arrays,jldopen(f,"r")["df"]);
    catch e
        println(e)
    end
end

# checking if all .jld2 with different seeds have same number of simulations
all([ s == size(df_arrays[1], 1) for s in size.(df_arrays, 1)])

# variables that do not need statistics
vars_nostat = ["m", "μ", "D", "graph", "rθ", "soptim", "σ_α"]
#variables needing statistics
vars = ["N", "alphas", "betas", "alphau", "betau", "Q_ST_u", "Q_ST_s"]
[df[!,"Q_ST_u"] = df.betau ./ (df.betau .+ df.alphau) for df in df_arrays]
[df[!,"Q_ST_s"] = df.betas ./ (df.betas .+ df.alphas) for df in df_arrays]

#filling df_aggreg
[sort!(df,["m", "μ", "D","rθ","σ_α"]) for df in df_arrays]
# making sure that all df correspond.
@assert all([all(df_arrays[1][:,vars_nostat] == df[:,vars_nostat]) for df in df_arrays])

mtemp = [df[:,vars] |> Matrix for df in df_arrays]
_stats = hcat(df_arrays[1][:,vars_nostat] |> Matrix, mean(mtemp), std(mtemp))
df_aggreg = DataFrame(_stats,[vars_nostat; vars.*"_mean"; vars.*"_std"])

# filtering nan etcs
for v in vars.*"_mean"
    println("we found ", count(isnan.(df_aggreg[:,v])), " NaN for $v")
    filter!(v => x -> !(ismissing(x) || isnothing(x) || isnan(x)), df_aggreg)
end

#################################
### matching df with graphs_df ##
#################################
@load "../../../graphs_utils/M=7/graph_prop_M=7.jld2" graphs_df
println("matching df and graphs_df")
metrics =  names(graphs_df[:,2:end])



for _n in metrics
    df_aggreg[!,_n] = zeros(size(df_aggreg,1))
end

global i = 0
for r in eachrow(df_aggreg)
    for _i in 1:size(graphs_df,1)
        if graphs_df.graph[_i] == r.graph
            r[metrics] = graphs_df[_i,2:end]
            global i+= 1
            break
        end
    end
end
println("df matched $i / $(size(df_aggreg,1))")

@save "$name_sim/$(name_sim)_$(date_sim)_aggreg.jld2" df_aggreg
using CSV
CSV.write("$name_sim/$(name_sim)_$(date_sim)_aggreg.csv", df_aggreg[:,Not("graph")])