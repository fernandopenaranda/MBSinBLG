module MBSinBLG
# Use README as the docstring of the module:
@doc read(joinpath(dirname(@__DIR__), "README.md"), String) Quantica

using Quantica, StaticArrays, Parameters, LinearAlgebra, StatsBase
using Baselet, Arpack, StaticArrays

using Requires
# function __init__()
#       @require GLMakie = "e9467ef8-e4e7-5192-8a1a-b1aee30e663a" include("plot_makie.jl")
#       @require VegaLite = "112f6efa-9a02-5b7d-90c0-432ed331239a" include("plot_vegalite.jl")
# end

include("model.jl")
include("random_interface.jl")

end
