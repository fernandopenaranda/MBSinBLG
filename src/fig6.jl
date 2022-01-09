function runfig6()
    p = Params(Ln = 2000, W = 2000, Ls = 20, Ws = 20, scale = 40, λ = 5, α = 0,   
    Δ = 1, d = 0, τ = 1)
    
    # PANELS a-b
    p = reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6)
    presets_fig = Fig4_presets(0.000, 0π/180, 0, 8, true)
    data = ldosonlattice_averaged_sc(p, presets_fig)
    savepsi("fig6psi", p, data[1])
    savepsi("fig6psi", p, data[2])

    p = reconstruct(p, EZ = SA[0, 4, 0], μN = 1.192)
    presets_fig = Fig4_presets(0, 0π/180, 0, 4, true)
    data = ldosonlattice_averaged_sc(p, presets_fig)
    savepsi("fig6psi", p, data[1])

   # PANELS c-d
    p = reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6)
    presets_fig = Fig4_presets(0.002, 2π/180, 9, 8, true)
    data = ldosonlattice_averaged_sc(p, presets_fig)
    savepsi("fig6psi", p, data[1])
    savepsi("fig6psi", p, data[2])

    p = reconstruct(p, EZ = SA[0, 4, 0], μN = 1.192)
    presets_fig = Fig4_presets(0.002, 2π/180, 9, 4, true)
    data = ldosonlattice_averaged_sc(p, presets_fig)
    savepsi("fig6psi", p, data[1])

#############

    p = Params(Ln = 1500, W = 1500, Ls = 20, Ws = 20, scale = 40, λ = 5, α = 0,   
    Δ = 1, d = 0, τ = 1)

    p = reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6)
    presets_fig = Fig4_presets(0.01, 2π/180, 0, 8, true)
    data = ldosonlattice_averaged_sc(p, presets_fig)
    savepsi("fig6psi3", p, data[1])
    savepsi("fig6psi3", p, data[2])

    p = reconstruct(p, EZ = SA[0, 4, 0], μN = 1.192)
    presets_fig = Fig4_presets(0.01, 2π/180, 0, 4, true)
    data = ldosonlattice_averaged_sc(p, presets_fig)
    savepsi("fig6psi3", p, data[1])


    # PANELS e-f
    sp = splittingvsrotation(reconstruct(p, Ln = 400, EZ = SA[0, 2, 0], μN = 0.6), 0:.25:10, true, ϕ0 = pi)
    savespectrum("fig6", p, sp[1], sp[2])

    sp = splittingvsrotation(reconstruct(p, Ln = 400, EZ = SA[0, 4, 0], μN = 1.192), 0:.25:10, true, ϕ0 = pi)
    savespectrum("fig6", p, sp[1], sp[2])

    sp = splittingvsrotation(reconstruct(p, Ln = 400, EZ = SA[0, 4, 0], μN = 0.6), 0:.5:10, true, ϕ0 = pi)
    savespectrum("fig6", p, sp[1], sp[2])

    # Study vs spectrum vs disorder (not included in the manuscript)
    # sp = splittingvsdisorder(reconstruct(p, Ln = 800, W = 2000, EZ = SA[0, 2, 0]), 0:0.02:0.5, true, ϕ0 = pi)
    # savespectrum("fig6testdisorder2", p, sp[1], sp[2])

    # sp = splittingvsdisorder(reconstruct(p, Ln = 800, W = 2000, EZ = SA[0, 5, 0]), 0:0.02:0.5, true, ϕ0 = pi)
    # savespectrum("fig6testdisorder5", p, sp[1], sp[2])
end


function sims()
    p = Params(Ln = 1000, W = 2000, Ls = 20, Ws = 20, scale = 40, λ = 5, α = 0,   
    Δ = 1, d = 0, τ = 1)

    # p = reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6)
    # presets_fig = Fig4_presets(0.0, 0*π/180, 0, 8, true)
    # data = ldosonlattice_averaged_sc(p, presets_fig)
    # savepsi("fig6psieta0rot2_BDI_red_1000x2000_f32", p, data[1])
    # savepsi("fig6psieta0rot2_BDI_blue_1000x2000_f32", p, data[2])

    # p = reconstruct(p, EZ = SA[0, 4, 0], μN = 0.6)
    # presets_fig = Fig4_presets(0, 0π/180, 0, 4, true)
    # data = ldosonlattice_averaged_sc(p, presets_fig)
    # savepsi("fig6psieta0rot2_D_1000x2000_f32", p, data[1])

    for i in 1:10
        println(i)
        p = reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6)
        presets_fig = Fig4_presets(0.005, 2π/180, 0, 8, true)
        data = ldosonlattice_averaged_sc(p, presets_fig)
        savepsi("fig6psieta0.005rot2_BDI_red_1000x2000_f64", p, data[1])
        savepsi("fig6psieta0.005rot2_BDI_blue_1000x2000_f64", p, data[2])

        p = reconstruct(p, EZ = SA[0, 4, 0], μN = 0.6)
        presets_fig = Fig4_presets(0.005, 2π/180, 0, 4, true)
        data = ldosonlattice_averaged_sc(p, presets_fig)
        savepsi("fig6psieta0.005rot2_D_1000x2000_f64", p, data[1])

        # println(i)
        # p = reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6)
        # presets_fig = Fig4_presets(0.01, 2π/180, 0, 8, true)
        # data = ldosonlattice_averaged_sc(p, presets_fig)
        # savepsi("fig6psieta0.01rot2_BDI_red_1500x1500_f32", p, data[1])
        # savepsi("fig6psieta0.01rot2_BDI_blue_1500x1500_f32", p, data[2])
        # p = reconstruct(p, EZ = SA[0, 4, 0], μN = 1.192)
        # presets_fig = Fig4_presets(0.01, 2π/180, 0, 4, true)
        # data = ldosonlattice_averaged_sc(p, presets_fig)
        # savepsi("fig6psieta0.01rot2_D_1500x1500_f32", p, data[1])
        # p = reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6)
        # presets_fig = Fig4_presets(0.005, 2π/180, 0, 8, true)
        # data = ldosonlattice_averaged_sc(p, presets_fig)
        # savepsi("fig6psieta0.005rot2_BDI_red_1500x1500_f32", p, data[1])
        # savepsi("fig6psieta0.005rot2_BDI_blue_1500x1500_f32", p, data[2])
        # p = reconstruct(p, EZ = SA[0, 4, 0], μN = 1.192)
        # presets_fig = Fig4_presets(0.005, 2π/180, 0, 4, true)
        # data = ldosonlattice_averaged_sc(p, presets_fig)
        # savepsi("fig6psieta0.005rot2_D_1500x1500_f32", p, data[1])
        # p = reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6)
        # presets_fig = Fig4_presets(0.001, 2π/180, 0, 8, true)
        # data = ldosonlattice_averaged_sc(p, presets_fig)
        # savepsi("fig6psieta0.001rot2_BDI_red_1500x1500_f32", p, data[1])
        # savepsi("fig6psieta0.001rot2_BDI_blue_1500x1500_f32", p, data[2])
        # p = reconstruct(p, EZ = SA[0, 4, 0], μN = 1.192)
        # presets_fig = Fig4_presets(0.001, 2π/180, 0, 4, true)
        # data = ldosonlattice_averaged_sc(p, presets_fig)
        # savepsi("fig6psieta0.001rot2_D_1500x1500_f32", p, data[1])
    end
end