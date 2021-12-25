function runfig6()
    p = Params(Ln = 1440, W = 1440, Ls = 0, scale = 40, λ = 5, α = 0,   
    μN = 0.6, Δ = 1, d = 0, τ = 1)
     # p = Params(Ln = 1440, W = 1440, Ls = 0, scale = 40, λ = 5, α = 0,   
    # μN = 0.6, Δ = 1, d = 0, τ = 1)

    p = reconstruct(p, EZ = SA[0, 2, 0])
    presets_fig = Fig4_presets(0.0, 0.0, 0, 8, true)
    data = ldosonlattice_averaged_sc(p, presets_fig)
    savepsi("fig6", p, data[1])
    savepsi("fig6", p, data[2])

    p = reconstruct(p, EZ = SA[0, 4, 0])
    presets_fig = Fig4_presets(0.0, 0.0, 0, 4, true)
    data = ldosonlattice_averaged_sc(p, presets_fig)
    savepsi("fig6", p, data[1])
    
    # p = reconstruct(p, EZ = SA[0, 2, 0])
    # presets_fig = Fig4_presets(0.005, π/180, 0, 8, true)
    # data1 = ldosonlattice_averaged_sc(p, presets_fig)
    # data2 = ldosonlattice_averaged_sc(p, presets_fig)
    # data3 = ldosonlattice_averaged_sc(p, presets_fig)
    # data4 = ldosonlattice_averaged_sc(p, presets_fig)
    # data5 = ldosonlattice_averaged_sc(p, presets_fig)

    p = reconstruct(p, EZ = SA[0, 2, 0])
    presets_fig = Fig4_presets(0.005, π/180, 19, 8, true)
    data = ldosonlattice_averaged_sc(p, presets_fig)
    savepsi("fig6", p, data[1])
    savepsi("fig6", p, data[2])

    p = reconstruct(p, EZ = SA[0, 4, 0])
    presets_fig = Fig4_presets(0.005, π/180, 19, 4, true)
    data = ldosonlattice_averaged_sc(p, presets_fig)
    savepsi("fig6", p, data[1])

    
    sp = splittingvsrotation(reconstruct(p, Ln = 800, W = 2000, EZ = SA[0, 2, 0]), 0:.25:5, true, ϕ0 = pi)
    savespectrum("fig6", p, sp[1], sp[2])


    sp = splittingvsrotation(reconstruct(p, Ln = 800, W = 2000, EZ = SA[0, 5, 0]), 0:.25:5, true, ϕ0 = pi)
    savespectrum("fig6", p, sp[1], sp[2])


    sp = splittingvsdisorder(reconstruct(p, Ln = 400, W = 2000, EZ = SA[0, 2, 0]), 0:0.02:0.5, true, ϕ0 = pi)
    savespectrum("fig6testdisorder2", p, sp[1], sp[2])

    sp = splittingvsdisorder(reconstruct(p, Ln = 400, W = 2000, EZ = SA[0, 5, 0]), 0:0.02:0.5, true, ϕ0 = pi)
    savespectrum("fig6testdisorder5", p, sp[1], sp[2])
   
end


# ## load and plot
# psi = readfig3psi(path);
# ldosonlattice
