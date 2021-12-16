############################################################################################
# Figure 2
############################################################################################

figure2(datapath::String; kw...) = 
    figure2(CSV.read(datapath, DataFrame,  delim = '\t'); kw...)
function figure2(data; ylims = (-10, 10))
     #remark: if you change the chemical potential in fig2 change the hline
    noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
    noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")
    scene, layout = layoutscene(30, resolution = (700, 700),
        font = noto_sans, fontsize = 22) 
    tag_list = unique(data.tag)
    println(tag_list)
    numbands = maximum(data.psi)
    sub = GridLayout()
    axes = sub[1:2, 1:2] = [Axis(scene, xlabel = L"k_y ⋅ a_0", ylabel = L"E \quad  [meV]", 
        xlabelsize = 24, ylabelsize = 24) for i in 1:2 for j in 1:2]
    tightlimits!.(axes)
    layout[1, 1] = sub
    #left
    customgray = Colors.RGB(.5, .5, .5)
    for (j, ax) in enumerate(axes)
        ax.xlabel = #L"k_y \dot a_0"
        ax.ylabel = "E (meV)"
        tag_indexes = findall(i -> i == tag_list[j], data.tag)
        for i in 1:numbands
            ψ_indexes = findall(x -> x == i, data.psi)
            selection = intersect(tag_indexes, ψ_indexes)

            if j == 1 || j == 3
                kcros = 0.046
                datak = data.karray[selection]
                lines!(ax, data.karray[selection][datak.<=-kcros+.001],
                      data.ϵs[selection][datak.<=-kcros+.001], 
                    color = ifelse(j != 1000 && 
                    intersect(i,(-1+Int(numbands/2):2+Int(numbands/2)))!=[],
                    :red, customgray),fxaa = true)

                lines!(ax, data.karray[selection][datak.>=kcros-.001],
                     data.ϵs[selection][datak.>=kcros-.001], 
                    color = ifelse(j != 1000 && 
                    intersect(i,(-1+Int(numbands/2):2+Int(numbands/2)))!=[],
                    :red, customgray),fxaa = true)
                
                lines!(ax, datak[abs.(datak) .< kcros], 
                     (data.ϵs[selection])[abs.(datak).<kcros], 
                    color = ifelse(j != 1000 && 
                    intersect(i,(-3+Int(numbands/2): 4+Int(numbands/2)))!=[],
                    :blue, customgray),fxaa = true)
                auxvect=collect(0:-2+Int(numbands/2))
                append!(auxvect, collect(3+Int(numbands/2):numbands))
                lines!(ax, datak[datak .> kcros], 
                    (data.ϵs[selection])[datak.>kcros], 
                    color = ifelse(intersect(i, auxvect)!=[], 
                    (customgray, 1), (customgray, 0)), fxaa = true)
                lines!(ax, datak[datak .< -kcros],
                     (data.ϵs[selection])[datak.<-kcros], 
                    color = ifelse(intersect(i, auxvect)!=[], 
                    (customgray, 1), (customgray, 0)), fxaa = true)           
            else
                if j == 2
                    kcrosr = 0.565
                    kcrosl = 0.435
                    ϵ = -10
                else
                    kcrosr = 0
                    kcrosl = 0
                    ϵ = 0
                end

                datak = range(0.2, stop = 0.8, length = length(data.ϵs[selection]))
                lines!(ax, datak[datak.<= kcrosl], 
                     data.ϵs[selection][datak.<= kcrosl], color = ifelse( 
                    intersect(i,(-1+Int(numbands/2):2+Int(numbands/2)))!=[],
                    :red, (customgray,0)), fxaa = true)
                lines!(ax, datak[datak.>= kcrosr], 
                     data.ϵs[selection][datak.>= kcrosr], color = ifelse( 
                    intersect(i,(-1+Int(numbands/2):2+Int(numbands/2)))!=[],
                    :red, (customgray,0)), fxaa = true)
                lines!(ax, datak[kcrosl+ ϵ .<= datak .<= kcrosr- ϵ], 
                    data.ϵs[selection][kcrosl + ϵ .<=datak.<= kcrosr- ϵ],
                    color = ifelse(intersect(i,(-3+Int(numbands/2):4+Int(numbands/2)))!=[],
                        :blue, customgray), fxaa = true)

                auxvect = [-3,-2,3,4] .+ Int(numbands/2)
                inter = intersect(i, auxvect)
                lines!(ax, datak[datak .>= kcrosr ], 
                     (data.ϵs[selection])[datak.>=kcrosr ], 
                    color = ifelse(inter!=[], (:blue, 1), 
                        ifelse(findall(x -> x == i, [-1,0,1,2] .+ Int(numbands/2)) == [], 
                        (customgray, 1), (:red,1))), fxaa = true)
                lines!(ax, datak[datak .<= kcrosl],
                    (data.ϵs[selection])[datak .<= kcrosl], 
                    color = ifelse(inter!=[], (:blue, 1), 
                        ifelse(findall(x -> x == i, [-1,0,1,2] .+ Int(numbands/2)) == [],
                        (customgray, 1), (:red,1))), fxaa = true)
            end
                if j == 3 && i == 1 || j == 4 && i == 1
                    hlines!(ax, 2.4, xmax = 10, color = :purple, linestyle= :dash, linewidth=1.8)
                else nothing end
        end
        ax.aspect = AxisAspect(1)
        if j == 1 || j == 3
             xlims!(ax, (-0.11,0.125)) 
             ax.xticks = ([-0.1, 0., 0.1], ["-0.1","0.0", "0.1"])
             ax.yticks = ([-5,0,5])
        else 
            xlims!(ax, (0.22, 0.807)) 
            ax.yticks = ([-5,0,5])
            ax.xticks = ([0.3, 0.5, 0.7])
        end
        ylims!(ax, ylims) 
  

    end
     hidexdecorations!(axes[1], grid = false)
     hidexdecorations!(axes[2], grid = false)
     hideydecorations!(axes[2], grid = false)
     hideydecorations!(axes[4], grid = false)
     axes[1].title = "Armchair"
     axes[2].title = "Zigzag"
     colgap!(sub, 1, Relative(-0.01))
    return scene
end

############################################################################################
# Figure 3
############################################################################################
function figure3(path)
    f = Figure()
    ax = Axis(f[1,1], xlabel = "Ez [meV]", ylabel = "μN [mev]")
    lines!(ax, boundarySZ, dμs)
    lines!(ax, boundarySA, dμs)
    return f
end
############################################################################################
# Figure 4
############################################################################################

function figure4(path)
    Ez, e = readspectrum(path)
    return figure4(Ez.EZ, Matrix(e))
end
function figure4(Ez, ea,  ylims = (-.15,.15))
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

#___________________________________________________________________________________________
# LDOS on lattice Plot 

ldosonlattice(psi, h0, scheme = "blues") = vlplot(h0, psi, 
    sitesize = DensityShader(), sitecolor = DensityShader(), siteopacity = 0.5,
    mindiameter = 5e-2, size = (300), colorscheme = scheme, sitestroke = nothing, #reds
    maxdiameter =30, plotlinks = false, discretecolorscheme = :gray, 
    labels = ("x (nm)", "y (nm)")) 

# logldosonlattice(psi, h0) = vlplot(h0, -log.(abs.(psi)), 
#     sitesize = DensityShader(), sitecolor = DensityShader(), siteopacity = 0.5,
#     mindiameter = maximum(-log.(abs.(psi)))*4/5, size = (300), colorscheme = "blues", sitestroke = nothing,
#     maxdiameter =maximum(-log.(abs.(psi))), plotlinks = false, discretecolorscheme = :gray, 
#     labels = ("x (nm)", "y (nm)")) 

latticeplot(psi, h0) = vlplot(h0, psi, 
    sitesize = DensityShader(), sitecolor = DensityShader(), siteopacity = 0.5,
    mindiameter = 1e-1, size = (300), colorscheme = "greys", sitestroke = nothing,
    maxthickness = 0.001, plotlinks = false, discretecolorscheme = :gray, 
    labels = ("x (nm)", "y (nm)"))
#___________________________________________________________________________________________



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
############################################################################################
# Figure 5
############################################################################################

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

############################################################################################
# Figure 6
############################################################################################

function figure6(patha, pathb)
    Ez, ea = readspectrum(patha)
    _, eb = readspectrum(pathb)
    return figure6(Ez.EZ, Matrix(ea), Matrix(eb))
end

function figure6(Ez, ea, eb,  ylims = (-0.4, 0.4))
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

############################################################################################
# Figure 7
############################################################################################

function figure7(path; kw...)
    ϕs, ea, eb, ec, ed = readcpr(path)
    return figure7(ϕs.phi, Matrix(ea), Matrix(eb), Matrix(ec), Matrix(ed); kw...)
end

function figure7(ϕs, ea, eb, ec, ed; ylims = (-0.1, 0.1))
    
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




#########
# Extra
########

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
    return vlplot(b, ylims = (-10, 10), size = (300))
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
    return a
end