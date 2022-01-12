module MBSinBLG
    # Use README as the docstring of the module:
    @doc read(joinpath(dirname(@__DIR__), "README.md"), String) MBSinBLG


    using Requires

    function __init__()
        @require CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0" begin
            @require  LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f" begin 
                @require  VegaLite = "112f6efa-9a02-5b7d-90c0-432ed331239a" begin 
                    @require  Colors = "5ae59095-9a9b-59fe-a467-6f913c188581" include("figplots.jl")
                end
            end
        end
    end

    using SharedArrays, Distributed
    using Quantica, StaticArrays, Parameters, LinearAlgebra, StatsBase
    using Baselet, Arpack, StaticArrays, CSV, DataFrames, Dates
    using Optim, ProgressMeter
    
    using PhysicalConstants.CODATA2018: ustrip, @u_str, ħ, k_B, m_e, e, μ_B

    export runfig2, runfig3, runfig4, runfig5, runfig6, runfig7 
    export fig2plot, fig3plot, fig4plot, fig5plot, fig6plot, fig7plot, ldosonlattice

    export nanoribbonS, nanoribbonSA, nanoribbonSZ, Params,  modelS, rectangle_weaklink, 
        rectangle_randombounds_sc, ldosonlattice_averaged_sc

    include("model.jl")
    include("nanoribbon.jl")
    include("bounded_sys.jl")
    include("ldos.jl")
    include("spectrum.jl")
    include("squid.jl")
    include("save.jl")
    include("fig2.jl")
    include("fig3.jl")
    include("fig4.jl")
    include("fig5.jl")
    include("fig6.jl")
    include("fig7.jl")
end