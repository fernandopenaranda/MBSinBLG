module MBSinBLG
    # Use README as the docstring of the module:
    @doc read(joinpath(dirname(@__DIR__), "README.md"), String) MBSinBLG

    using Quantica, StaticArrays, Parameters, LinearAlgebra, StatsBase
    using Baselet, Arpack, StaticArrays, CSV, DataFrames, Dates
    using Optim, ProgressMeter, Distributed, SharedArrays
    using PhysicalConstants.CODATA2018: ustrip, @u_str, ħ, k_B, m_e, e, μ_B
    
    using Colors, CairoMakie, ElectronDisplay, LaTeXStrings, VegaLite

    export nanoribbonS, nanoribbonSA, nanoribbonSZ, Params, rectangle_squid, 
        rectangle_randombounds_sc, splittingvsrotation, savespectrum, modelS,
        cpr, savepsi, ldosonlattice_averaged_sc, savecprs

    export fig2run, fig3run, fig4run, fig5run, fig6run, fig7run, ldosonlattice
    export fig2plot, fig3plot, fig4plot, fig5plot, fig6plot, fig7plot

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
    include("figplots.jl")
end