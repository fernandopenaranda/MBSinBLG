using Quantica, ArnoldiMethod, Arpack

include("model.jl")
include("plots.jl")
include("random_interface.jl")

using CSV, Dates, DataFrames

function runfig4() #optional: reduce W to make the calculation faster
    p = Params(Ln = 2000, Ls = 0, scale = 40, λ= 5, α = 0, 
        EZ = SA[0, 0.6, 0], μN = 1.192, Δ = 0.3, d = 0, τ = 1);

    presets_fig = Fig3_presets(0.0, 0.0, 0, 4, false)
    data = ldosonlattice_averaged_sc(p, presets_fig)
    savepsi("fig4", p, data[1])
        
    presets_fig = Fig3_presets(0.0, 0.0, 0, 4, true)
    data = ldosonlattice_averaged_sc(p, presets_fig)
    savepsi("fig4", p, data[1])
        
    sp = spectrumsweepb(p, 0:.1:2, true)
    savespectrum("fig4", p, sp[1], sp[2])

    sp = spectrumsweepb(p, 0:.1:2, false)
    savespectrum("fig4", p, sp[1], sp[2])
end
# Data
patha = "/Users/fernandopenarandadelrio/Desktop/diracfuse/BLG/data/fig4/2021-12-05T14:11:53.482"


function figure_4(path)
    Ez, e = readspectrum(path)
    return figure_4(Ez.EZ, Matrix(e))
end
function figure_4(Ez, ea,  ylims = (-.15,.15))
    fig = Figure(resolution = (500, 250), font = "Times New Roman") 
    axa = Axis(fig[1, 1],  xlabel = L"$E_Z$ [meV]", xaxisposition = :top)
    lines!(axa, Ez, ea[1,:])
    mean = size(ea, 1) ÷ 2
    println(mean)
    for i in 2:size(ea, 1)
        lines!(axa, Ez, ea[i,:], color = ifelse(in(i, mean-1:mean+2), 
            (:dodgerblue3, 1), :gray), opacity = .5)
    end
    vlines!(axa, [0.6], color = :black, linewidth = .9, linestyle = :dash)
    axa.ylabel = "E [meV]"
    ylims!(axa, ylims)
    xlims!(axa, low = 0, high = 2)
    fig
end



#########################################
# # LDOS with KPM (deprecated)

# p = Params(Ln = 2000, Ls = 0, scale = 40, λ= 5, α = 2, 
# EZ = SA[0, 0.6, 0], μN = 1.192, Δ = 0.3, d = 0, τ = 1);

# index = 33755
# presets_fig = Fig3_presets(0.0, 0.0, 0, 4, false)
# ldosa = ldos(p, index, presets_fig.η, presets_fig.angle)
# presets_fig = Fig3_presets(0.0, 0.0, 0, 4, true)
# ldosb = ldos(p, index, presets_fig.η, presets_fig.angle)

# # store
# adata = DataFrame(x = ldosa[1], y = ldosa[2])
# bdata = DataFrame(x = ldosb[1], y = ldosb[2])
# CSV.write("data/fig3_ldosa", adata; delim = '\t')
# CSV.write("data/fig3_ldosb", bdata; delim = '\t')

# # plot
# data = CSV.read(joinpath(datapath), DataFrame, delim='\t')
# panelrightfig3(q[1], adata.y, bdata.y, cdata.y, ddata.y)