using Distributed 
addprocs(2)
@everywhere using MBSinBLG
# using MKL
# BLAS.set_num_threads(2)
@everuwhere using MKL
@everywhere using Quantica, StaticArrays, Parameters, LinearAlgebra, StatsBase
@everywhere using Baselet, Arpack, StaticArrays, CSV, DataFrames, Dates
@everywhere using Optim, ProgressMeter, Distributed
@everywhere using PhysicalConstants.CODATA2018: ustrip, @u_str, ħ, k_B, m_e, e, μ_B