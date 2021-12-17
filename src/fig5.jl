function runfig5()
    p = Params(Ln = 1440, W = 1440, Ls = 0, scale = 40, 
        λ = 5, α = 0, EZ = SA[0, 0, 0],Δ = 1, d = 0, τ = 1);

    p = reconstruct(p, μN = 0.000001)
    @time sp = spectrumsweepb(p, 0:.25:5, true)
    savespectrum("fig5", p, sp[1], sp[2])

    p = reconstruct(p, μN = 0.6)
    sp = spectrumsweepb(p, 0:.25:5, true)
    savespectrum("fig5", p, sp[1], sp[2])

    p = reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6)
    sp = spectrumsweepα(p, 0:.25:4, true)
    savespectrum("fig5", p, sp[1], sp[2])

    p = reconstruct(p, EZ = SA[0, 4, 0], μN = 0.6)
    sp = spectrumsweepα(p, 0:.25:4, true)
    savespectrum("fig5", p, sp[1], sp[2])
end