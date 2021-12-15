############################################################################################
# Figure 2
############################################################################################

fig2(datapath::String; kw...) = 
    fig2(CSV.read(datapath, DataFrame,  delim = '\t'); kw...)
function fig2(data; ylims = (-10, 10))
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
