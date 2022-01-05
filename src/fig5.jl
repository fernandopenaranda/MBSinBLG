function runfig5()
    p = Params(Ln = 2000, W = 2000, Ls = 20, Ws = 20, scale = 40, λ = 5, α = 0,   
    Δ = 1, d = 0, τ = 1)

    p = reconstruct(p, μN = 0.000001)
    sp = spectrumsweepb(p, 0:.25:5, true)
    savespectrum("fig5", p, sp[1], sp[2])

    p = reconstruct(p, μN = 0.6)
    sp = spectrumsweepb(p, 0:.25:5, true)
    savespectrum("fig5", p, sp[1], sp[2])

    p = reconstruct(p, μN = 1.192)
    sp = spectrumsweepb(p, 0:.25:5, true)
    savespectrum("fig5", p, sp[1], sp[2])

    p = reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6)
    sp = spectrumsweepα(p, 0:.25:4, true)
    savespectrum("fig5", p, sp[1], sp[2])

    p = reconstruct(p, EZ = SA[0, 4, 0], μN = 1.192)
    sp = spectrumsweepα(p, 0:.25:4, true)
    savespectrum("fig5", p, sp[1], sp[2])
end
