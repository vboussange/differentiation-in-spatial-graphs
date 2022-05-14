#=
Contains functions to calculate population size and differentiation metrics
from PDE solutions.
=#
using StatsBase
function αdiv(sol,p)
    ws = FrequencyWeights.(eachrow(sol))
    mygrid = p["S"]
    mean([var(mygrid,w,corrected=false) for w in ws])
end

function βdiv(sol,p)
    ws = FrequencyWeights.(eachrow(sol))
    mygrid = p["S"]
    var([mean(mygrid,w) for w in ws],corrected=false)
end

function γdiv(sol,p)
    ws = FrequencyWeights(vcat((eachrow(sol))...))
    mygrid = vcat(eachrow(repeat(p["S"]',p["M"],1))...)
    var(mygrid[:],ws[:],corrected=false)
end

function popsize(sol,p,loc = false)
    if size(sol,2) > 1
        dS = p["dS"]
        _s = 0.5 .* (sol[:,1] .+ sol[:,end] )
        _s .+= sum(sol[:,2:end-1],dims=2)[:] .* dS
        if loc
             return _s
         else
             return sum(_s)
         end
    else
        return sum(sol)
    end
end