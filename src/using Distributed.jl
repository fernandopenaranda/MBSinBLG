using Distributed 
# addprocs(3)
@everywhere using MBSinBLG
# using MKL
# BLAS.set_num_threads(2)
@everywhere using MKL, SharedArrays
@everywhere using Quantica, StaticArrays, Parameters, LinearAlgebra, StatsBase
@everywhere using Baselet, Arpack, StaticArrays, CSV, DataFrames, Dates
@everywhere using Optim, ProgressMeter, Distributed
@everywhere using PhysicalConstants.CODATA2018: ustrip, @u_str, ħ, k_B, m_e, e, μ_B


using MBSinBLG, Quantica, StaticArrays, Parameters, LinearAlgebra, StatsBase, SharedArrays
using Baselet, Arpack, StaticArrays, CSV, DataFrames, Dates
using PhysicalConstants.CODATA2018: ustrip, @u_str, ħ, k_B, m_e, e, μ_B
using Optim, ProgressMeter