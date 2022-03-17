#=
This script explores statistical properties of the bunch of graphs 
generated with habitat dataset.
=#

cd(@__DIR__)
using LightGraphs
using DataFrames
using CSV
using GLM
using Statistics
using MLBase
using JLD2
using Printf

@load "graph_prop_realistic_graphs.jld2" graphs_df

metrics = names(graphs_df)[2:end]
# we remove mean betweenness centrality as its correlation is one with cl
metrics = vcat(metrics[.!(metrics .∈ Ref([
                                        "mean_betweenness_centrality",
                                        "heat_hetero",
                                        "sqrtk_inv",
                                        "kk"]))])
                        
######################################
##### working clustering with R ######
######################################
dist = 1 .- cor(graphs_df[:,metrics] |> Array)

using RCall
R"""
dist_mat <- $(dist)
rownames(dist_mat) <- $metrics
colnames(dist_mat) <- $metrics
dist <- as.dist(dist_mat)
hc <- hclust(dist)
phc <- plot(hc)
"""