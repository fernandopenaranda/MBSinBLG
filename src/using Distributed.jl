using Distributed 
BLAS.set_num_threads(2)
@everywhere using Quantica, StaticArrays, Parameters, LinearAlgebra, StatsBase
@everywhere using Baselet, Arpack, StaticArrays, CSV, DataFrames, Dates
@everywhere using Optim, ProgressMeter, Distributed
@everywhere using PhysicalConstants.CODATA2018: ustrip, @u_str, ħ, k_B, m_e, e, μ_B
@everywhere using MKL