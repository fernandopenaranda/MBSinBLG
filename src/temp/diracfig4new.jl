#understanding the low energy low field transition

using Distributed
addprocs(3)

@everywhere include("model.jl")
@everywhere include("random_interface.jl")
@everywhere include("plots.jl")
@everywhere include("Plot_functions.jl")


p = Params(Ln = 100, Ls = 0, scale = 40, λ= 5, α = 0, 
    EZ = SA[0, 0.6, 0], μN = 1.192, Δ = 0.3, d = 0, τ = 1);


 sp = spectrumsweepb(p, 0:.25:2, true)

 