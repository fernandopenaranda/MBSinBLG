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