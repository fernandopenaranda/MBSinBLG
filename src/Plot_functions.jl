include("model.jl")
using VegaLite, CairoMakie
#___________________________________________________________________________________________
# LDOS on lattice Plot

ldosonlattice(psi, h0) = vlplot(h0, psi, 
    sitesize = DensityShader(), sitecolor = DensityShader(), siteopacity = 0.5,
    mindiameter = 1e-1, size = (300), colorscheme = "blues", sitestroke = nothing, 
    maxthickness = 0.005, plotlinks = false, discretecolorscheme = :gray,
    labels = ("x (nm)", "y (nm)"))

latticeplot(psi, h0) = vlplot(h0, psi, 
    sitesize = DensityShader(), sitecolor = DensityShader(), siteopacity = 0.5,
    mindiameter = 1e-1, size = (300), colorscheme = "greys", sitestroke = nothing,
    maxthickness = 0.001, plotlinks = false, discretecolorscheme = :gray, 
    labels = ("x (nm)", "y (nm)"))
#___________________________________________________________________________________________
# Bands

function normalbands(p, which = "armchair")
    if which == "armchair"
        kmin, kmax = (-2., 2.)
        axis = (0, 1)
        h = nanoribbonN(p, axis);
        b = bandstructure(h, cuboid((-2.,2.), subticks = 201),
            method = ArpackPackage(sigma = -0.00001, nev = 16));
        # a = vlplot(b, ylims = (-10, 10), size = (300))
    else
        kmin, kmax = (0., 2pi)
        axis = (1, 0)
        h = nanoribbonN(p, axis);
        b = bandstructure(h, cuboid((0, 2pi), subticks = 101), 
            method = ArpackPackage(sigma = -0.00001, nev = 16));
    end
    return vlplot(b, ylims = (-10, 10), size = (300))#, vlplot(b, ylims = (-10, 10), size = (300))
end

function scbands(p, which = "armchair")
    if which == "armchair"
        axis = (0,1)
        h = nanoribbonS(p, axis);
        b = bandstructure(h, cuboid((0,2), subticks = 101),
            method = ArpackPackage(sigma = -0.00001, nev = 16));
        a = vlplot(b, ylims = (-2, 2), size = (300))
    else
        axis = (1,0)
        h = nanoribbonS(p, axis);
        b = bandstructure(h, cuboid((0, 2pi), subticks = 301), 
            method = ArpackPackage(sigma = -0.00001, nev = 16));
        a = vlplot(b, ylims = (-10, 10), size = (300))
    end
    return a#, a2
end

#vlplot(h, size = (300), maxdiameter = 1, mindiameter = 0, maxthickness = 0.4, 
#   colorscheme = "greys", plotlinks = true, sitestroke = nothing, plotsites = true,
#   linkcolor = -0.3)

########
function spectrumvsmu(e, mu, ylims = (-.2,.2))
    fig = Figure(resolution = (500, 300), font = "CMU Serif") 
    ax = Axis(fig[1,1])
    lines!(ax, mu, e[1,:])
    for i in 2:size(e, 1)
        lines!(ax, mu, e[i,:])
    end
    ax.xlabel = "mu [meV]"
    ax.ylabel = "E [meV]"
    ylims!(ax, ylims)
    fig
end

function spectrumvsλ(e, λ, ylims = (-.2,.2))
    fig = Figure(resolution = (500, 300), font = "CMU Serif") 
    ax = Axis(fig[1,1])
    lines!(ax, λ, e[1,:])
    for i in 2:size(e, 1)
        lines!(ax, λ, e[i,:])
    end
    ax.xlabel = "λ [meV]"
    ax.ylabel = "E [meV]"
    ylims!(ax, ylims)
    fig
end

function spectrumvsb(e, b, ylims = (-.2,.2))
    fig = Figure(resolution = (500, 300), font = "CMU Serif") 
    ax = Axis(fig[1,1])
    lines!(ax, b, e[1,:])
    for i in 2:size(e, 1)
        lines!(ax, b, e[i,:])
    end
    ax.xlabel = "EZ [meV]"
    ax.ylabel = "E [meV]"
    ylims!(ax, ylims)
    fig
end


function spectrumvsα(e, b, ylims = (-.2,.2))
    fig = Figure(resolution = (500, 300), font = "CMU Serif") 
    ax = Axis(fig[1,1])
    lines!(ax, b, e[1,:])
    for i in 2:size(e, 1)
        lines!(ax, b, e[i,:])
    end
    ax.xlabel = "α [meV]"
    ax.ylabel = "E [meV]"
    ylims!(ax, ylims)
    fig
end

function spectrumvsangle(e, b, ylims = (-.2,.2))
    fig = Figure(resolution = (500, 300), font = "CMU Serif") 
    ax = Axis(fig[1,1])
    b = b .*180/π
    lines!(ax, b, e[1,:])
    for i in 2:size(e, 1)
        lines!(ax, b, e[i,:])
    end
    ax.xlabel = " θ  [degrees]"
    ax.ylabel = "E [meV]"
    ylims!(ax, ylims)
    fig
end

#___________________________________________________________________________________________
# Figure 3
#___________________________________________________________________________________________
function plotldos(x, y, xrange = (-200,200))
    fig = Figure(resolution = (500, 300), font = "CMU Serif") 
    ax = Axis(fig[1,1])
    lines!(ax, x, y)
    xlims!(ax, xrange)
    fig
end

function panelrightfig3(x, y1, y2, y3, y4, xrange = (-1.9,1.9))
    ynorm1 = maximum([maximum(y1), maximum(y2)])
    ynorm2 = maximum([maximum(y3), maximum(y4)])
  
    scene, layout = layoutscene(30, size = (300,100))
    sub_c = GridLayout()
    axc = sub_c[1:2, 1] = [Axis(scene) for i in 1:2]
    layout[1, 1] = sub_c

    axc[1].xlabel = "LDOS [arb. units]"
    axc[1].ylabel = "E [meV]"
    axc[2].xlabel = "LDOS [arb. units]"
    axc[2].ylabel = "E [meV]"

    lines!(axc[1], x, y1./ynorm1, color = :blue)#, markersize = 5, strokewidth = 0)
    lines!(axc[1], x, y2./ynorm1, color = :black)
    lines!(axc[2], x, y3./ynorm2, color = :purple)
    lines!(axc[2], x, y4./ynorm2, color = :orange)  
    
    CairoMakie.rotate!(axc[1].scene, -π/2) 
    CairoMakie.rotate!(axc[2].scene, -π/2) 
    CairoMakie.ylims!(axc[1], xrange)
    CairoMakie.ylims!(axc[2], xrange)
    CairoMakie.xlims!(axc[1], (0,1.1))
    CairoMakie.xlims!(axc[2], (0,1.1))
    hidexdecorations!(axc[1], grid = false)   
    # rowgap!(sub_r, 1, Relative(0.0))
    colsize!(layout, 1, Relative(1/5))
    elem_1 = LineElement(color = :blue, strokecolor = :blue, strokewidth = 3)
    elem_2 = LineElement(color = :black, strokecolor = :black, strokewidth = 3)
    elem_3 = LineElement(color = :purple, strokecolor = :purple, strokewidth = 3)
    elem_4 = LineElement(color = :orange, strokecolor = :orange, strokewidth = 3)
    leg1 = Legend(axc[1].scene, [elem_1, elem_2], ["a)", "b)"],
      orientation = :horizontal, halign = 1.85, valign = 4.2, width =60, height = Auto(), 
      tellwidth = false, tellheight = false)
    leg2 = Legend(axc[2].scene, [elem_3, elem_4], ["c)", "d)"],
      orientation = :horizontal, halign = 1.85, valign = 4.2, width =60, height = Auto(), 
      tellwidth = false, tellheight = false)
     
    leg1.nbanks = 2
    leg2.nbanks = 2
    # axe.yreversed = true
    return scene
end
