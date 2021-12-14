function runfig6()
    p = Params(Ln = 1440, Ls = 0, scale = 40, λ = 5, α = 0,   
    μN = 0.6, Δ = 1, d = 0, τ = 1)

    # p = reconstruct(p, EZ = SA[0, 2, 0])
    # presets_fig = Fig3_presets(0.0, 0.0, 0, 8, true)
    # data = ldosonlattice_averaged_sc(p, presets_fig)
    # savepsi("fig6", p, data[1])
    # savepsi("fig6", p, data[2])

    # p = reconstruct(p, EZ = SA[0, 4, 0])
    # presets_fig = Fig3_presets(0.0, 0.0, 0, 4, true)
    # data = ldosonlattice_averaged_sc(p, presets_fig)
    # savepsi("fig6", p, data[1])
    
    # p = reconstruct(p, EZ = SA[0, 2, 0])
    # presets_fig = Fig4_presets(0.00, π/180, 5, 8, true)
    # data = ldosonlattice_averaged_sc(p, presets_fig)
    # savepsi("fig6", p, data[1])
    # savepsi("fig6", p, data[2])

    # p = reconstruct(p, EZ = SA[0, 4, 0])
    # presets_fig = Fig3_presets(0.005, π/180, 5, 4, true)
    # data = ldosonlattice_averaged_sc(p, presets_fig)
    # savepsi("fig6", p, data[1])

    sp = splittingvsrotation(reconstruct(p, EZ = SA[0, 2, 0]), 0:.25:5, true)
    savespectrum("fig6", p, sp[1], sp[2])

    sp = splittingvsrotation(reconstruct(p, EZ = SA[0, 4, 0]), 0:.25:5, true)
    savespectrum("fig6", p, sp[1], sp[2])

    
end

function figure_6(patha, pathb)
    Ez, ea = readspectrum(patha)
    _, eb = readspectrum(pathb)
    return figure_6(Ez.EZ, Matrix(ea), Matrix(eb))
end

function figure_6(Ez, ea, eb,  ylims = (-0.4, 0.4))
    fig = Figure(resolution = (500, 300), font = "Times New Roman") 
    axa = Axis(fig[1,1])
    axb = Axis(fig[2,1])
    mean = size(ea, 1) ÷ 2
    lines!(axa, 180/π .* Ez, ea[1,:])
    lines!(axb, 180/π .* Ez, eb[1,:])
    for i in 2:size(ea, 1)
        lines!(axa, 180/π .*Ez, ea[i,:], color = ifelse(in(i, mean-3:mean+4), 
        ifelse(in(i,[mean-1,mean+2]),  (:dodgerblue3, 1), (:firebrick3, 0.5)), :gray), opacity = .5)
        lines!(axb, 180/π .*Ez, eb[i,:], color = ifelse(in(i, mean-3:mean+4), 
        ifelse(in(i,[mean-1,mean+2]),  (:dodgerblue3, 1), (:firebrick3, 0.5)), :gray), opacity = .5)
    end
    axb.xlabel = " θ  [degrees]"
    axa.ylabel = "E [meV]"
    axb.ylabel = "E [meV]"
    hidexdecorations!(axa, grid = false)
    ylims!(axa, ylims)
    fig
end

# ## load and plot
# psi = readfig3psi(path);
