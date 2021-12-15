
function runfig7()
    p = Params(Ln = 400, Ls = 0, scale = 40, λ = 5, α = 0,   
        μN = 0.6, Δ = 1, d = 0, τ = 1, EZ = SA[0, 4, 0]);

    ϕlist = -2π:π/6:2π
    println("1")
    _, esa = cpr(ϕlist, reconstruct(p, EZ = SA[0, 2, 0], α = 0), nev = 16, W = 1440)
    println("2")
    _, esb = cpr(ϕlist, reconstruct(p, EZ = SA[0, 4, 0], α = 0), nev = 16, W = 1440)
    println("3")
    _, esc = cpr(ϕlist, reconstruct(p, EZ = SA[0, 2, 0], α = 1), nev = 16, W = 1440)
    println("4")
    ϕs, esd = cpr(ϕlist, reconstruct(p, EZ = SA[0, 4, 0], α = 1), nev = 16, W = 1440)
    savecprs("fig7", p, ϕs, esa, esb, esc, esd)
end

function figure_7(path; kw...)
    ϕs, ea, eb, ec, ed = readcpr(path)
    return figure_7(ϕs.phi, Matrix(ea), Matrix(eb), Matrix(ec), Matrix(ed); kw...)
end

function figure_7(ϕs, ea, eb, ec, ed; ylims = (-0.1, 0.1))
    
    fig = Figure(resolution = (500, 500), font = "Times New Roman") 
    axa = Axis(fig[1,1])
    axb = Axis(fig[1,2])
    axc = Axis(fig[2,1])
    axd = Axis(fig[2,2])
    mean = size(ea, 1) ÷ 2
    lines!(axa, ϕs/π, ea[1,:], color = :gray)
    lines!(axb, ϕs/π, eb[1,:], color = :gray)
    lines!(axc, ϕs/π, ec[1,:], color = :gray)
    lines!(axd, ϕs/π, ed[1,:], color = :gray)
    for i in 2:size(ea, 1)
        lines!(axa, ϕs/π, ea[i,:], color = ifelse(in(i, mean-3:mean+4),
        ifelse(in(i,[mean-1,mean+2]),  (:dodgerblue3, 1), (:firebrick3, 0.5)), :gray), opacity = .5)
        lines!(axb, ϕs/π, eb[i,:], color = ifelse(in(i, mean-3:mean+4),
        ifelse(in(i,[mean-1,mean+2]),  (:dodgerblue3, 1), (:firebrick3, 0.5)), :gray), opacity = .5)
        lines!(axc, ϕs/π, ec[i,:], color = ifelse(in(i, mean-3:mean+4),
        ifelse(in(i,[mean-1,mean+2]),  (:dodgerblue3, 1), (:firebrick3, 0.5)), :gray), opacity = .5)
        lines!(axd, ϕs/π, ed[i,:], color = ifelse(in(i, mean-3:mean+4),
        ifelse(in(i,[mean-1,mean+2]),  (:dodgerblue3, 1), (:firebrick3, 0.5)), :gray), opacity = .5)
    end
    axa.ylabel = "E [meV]"
    axa.xlabel = "Δϕ/π"
    axb.ylabel = "E [meV]"
    axb.xlabel = "Δϕ/π"
    axc.ylabel = "E [meV]"
    axc.xlabel = "Δϕ/π"
    axd.ylabel = "E [meV]"
    axd.xlabel = "Δϕ/π"
    ylims!(axa, ylims)
    ylims!(axb, ylims)
    ylims!(axc, ylims)
    ylims!(axd, ylims)

    xlims!(axa, (-2,2))
    xlims!(axb,  (-2,2))
    xlims!(axc,  (-2,2))
    xlims!(axd,  (-2,2))
    hidexdecorations!(axa, grid = false)
    hidexdecorations!(axb, grid = false)
    hideydecorations!(axb, grid = false)
    hideydecorations!(axd, grid = false)
    fig
end
