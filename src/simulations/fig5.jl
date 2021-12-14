function runfig5()
    p = Params(Ln = 1440, Ls = 0, scale = 40, 
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


function figure_5(patha, pathb, pathc, pathd)
    Ez, ea = readspectrum(patha)
    _, eb = readspectrum(pathb)
    α, ec = readspectrum(pathc)
    _, ed = readspectrum(pathd)
    return figure_5(Ez.EZ, α.EZ, Matrix(ea), Matrix(eb), Matrix(ec),  Matrix(ed))
end


function figure_5(Ez, α, ea, eb, ec, ed,  ylims = (-.4,.4))
    fig = Figure(resolution = (500, 300), font = "Times New Roman") 
    # gb = fig[1:2, 2] = GridLayout
    g = fig[1:2, 1:3] = GridLayout()
    # axa = Axis(fig[1, 1])
    # axb = Axis(fig[2, 2])
    # axc = Axis(fig[1, 3],  xlabel = L"$\alpha$", yaxisposition = :right)
    # axd = Axis(fig[2, 3],  xlabel = L"$\alpha$", yaxisposition = :right)
    axa = Axis(fig[1:2, 1],  xlabel = L"$E_Z$ [meV]")
    axb = Axis(fig[1:2, 2],  xlabel = L"$E_Z$ [meV]")
    axc = Axis(fig[1, 3],  xlabel = L"$\alpha$ [meV]", yaxisposition = :right)
    axd = Axis(fig[2, 3],  xlabel = L"$\alpha$ [meV]", yaxisposition = :right)
    lines!(axa, Ez, ea[1,:])
    lines!(axb, Ez, eb[1,:])
    mean = size(ea, 1) ÷ 2
    for i in 2:size(ea, 1)
        lines!(axa, Ez, ea[i,:], color = ifelse(in(i, mean-3:mean+4), 
            ifelse(in(i,[mean-1,mean+2]),  (:dodgerblue3, 1), (:firebrick3, 0.5)), :gray), opacity = .5)
        lines!(axb, Ez, eb[i,:], color = ifelse(in(i, mean-3:mean+4), 
            ifelse(in(i,[mean-1,mean+2]),  (:dodgerblue3, 1), (:firebrick3, 0.5)), :gray), opacity = .5)
    end
    vlines!(axb, [2, 4], color = :black, linewidth = .9, linestyle = :dash)

    lines!(axc, α, ec[1,:])
    lines!(axd, α, ed[1,:])
    mean = size(ec, 1) ÷ 2
    for i in 2:size(ec, 1)
        lines!(axc, α, ec[i,:], color = ifelse(in(i, mean-3:mean+4), 
            ifelse(in(i,[mean-1,mean+2]),  (:dodgerblue3, 1), (:firebrick3, 0.5)), :gray))
        lines!(axd, α, ed[i,:], color = ifelse(in(i, mean-3:mean+4), 
            ifelse(in(i,[mean-1,mean+2]),  (:dodgerblue3, 1), (:firebrick3, 0.5)), :gray))
    end
    
    axa.ylabel = "E [meV]"
    axb.ylabel = " "
    axc.ylabel = "E [meV]"
    axd.ylabel = "E [meV]"
    ylims!(axa, (-0.5,0.5))
    ylims!(axb,  (-0.5,0.5))
    ylims!(axc, ylims)
    ylims!(axd, (-0.15, 0.15))
    xlims!(axa, low = 0)
    xlims!(axb, low = 0)
    xlims!(axc, (0, 4))
    xlims!(axd, low = 0)
    xlims!(axd, (0, 4))
    colgap!(g, 20)
    rowgap!(g, 10)
    hidexdecorations!(axc, grid = false)
    hideydecorations!(axb, grid = false)
    # Label(g[1, 1, TopLeft()], "a", textsize = 16, font = "Times New Roman", 
    #     padding = (0, 5, 5, 0), halign = :right)
    # Label(g[1, 2, TopRight()], "b", textsize = 16, font = "Times New Roman", 
    #     padding = (0, 40, 0, 0), halign = :right)
    # Label(g[2, 2, TopRight()], "c", textsize = 16, font = "Times New Roman", 
    #     padding = (-0, 40, 0, 0), halign = :right)
    # colsize!(fig.layout, 1, Auto(0.5))
    fig
end
### QUICK PLOT
# # path = ""
# x, y = readspectrum(path)
# spectrumvsb(Matrix(y), x.EZ, (-0.4, 0.4))