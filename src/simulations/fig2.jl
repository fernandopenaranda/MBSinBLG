using Colors

function runfig2()
    p = Params(Ln = 1440, Ls = 0, scale = 40, U = 1e-7, λ = 10, α = 2, EZ = SA[0, 0, 0], 
        Δ = 0., μN = 0.);
    p = reconstruct(p, EZ = SA[0, 0, 0])
    p_EZ = reconstruct(p, EZ = SA[0, 2.8, 0])
    smat = fig2bands([p, p, p_EZ, p_EZ], 
        kmin = -0.15, kmax = .15, kpoints = 301, numbands = 56)
    writecsv(smat, [p, p, p_EZ, p_EZ], "fig1", ["1", "2", "3", "4"])
end


fig2bands(listp::Array{Params,1}; kw...) =
    [fig2bands(listp[i]; kw...) for i in 1:length(listp)]

function fig2bands(p, which = "armchair")
    kpoints = 301
    numbands = 56
    if which == "armchair"
        kmin, kmax = (-2., 2.)
        axis = (0, 1)
        h = nanoribbonN(p, axis);
        b = bandstructure(h, cuboid((-2.,2.), subticks = kpoints),
            method = ArpackPackage(sigma = -0.00001, nev = numbands));
    else
        kmin, kmax = (0., 2pi)
        axis = (1, 0)
        h = nanoribbonN(p, axis);
        b = bandstructure(h, cuboid((0, 2pi), subticks = kpoints), 
            method = ArpackPackage(sigma = -0.00001, nev = numbands));
    end
    karray = collect(range(kmin, kmax, length = kpoints))
    ϵs = zeros(ComplexF64, numbands, kpoints)
    for i in 1:length(karray)
        ϵs[:,i], _ = b.diag((karray[i],))
    end
    return karray, ϵs, b
end

floattostring(vec::Array{Float64,1}) = [floattostring(vec[i]) for i in 1:length(vec)]
floattostring(i::Float64) = "$(Int(i))"
function writecsv(vec, p, name, taglist)
    karray = collect(Iterators.flatten(
    [   collect(Iterators.flatten([append!([],vec[1][1]) for i in 1:size(vec[1][2],1)])) for j in 1:size(vec,1)]))
    ϵs = real(collect(Iterators.flatten(
        [collect(Iterators.flatten([append!([], vec[j][2][i,:]) for i in 1:size(vec[j][2], 1)])) for j in 1:size(vec,1)])))
    psi = collect(Iterators.flatten(
        [floattostring([ceil(i/size((vec[1][2]),2)) for i in 1:length(vec[1][2])]) for j in 1:size(vec,1)]))
    tag = [taglist[Int(ceil(i/length(vec[1][2])))] for i in 1:size(vec,1) * length(vec[1][2])]
    datafig = DataFrame(tag = tag, karray = karray, ϵs = ϵs, obs = obs, psi = psi)
    savecsv(name, p, datafig)
end

function savecsv(figname, p, data::DataFrame,xmin, xlen, xstep, ymin, ylen, ystep, kpoints) 
    #check and/or create folder structure
    #r= range(xmin, length = xlen, step = xstep) ok
    time_str = string(now())
    isdir("data") ? nothing : mkdir(data)
    isdir("data/$(figname)") ? nothing : mkdir("data/$(figname)")
    pathstr = join([pwd(), "/data/$(figname)/", time_str])
    mkdir(pathstr) 
    #build the csv.files
    dict = type2dict(p)
    push!(dict, :xlen => xlen, :xstep => xstep, :xmin => xmin, :ylen => ylen, :ystep => ystep, :ymin => ymin, :kpoints => kpoints)
    CSV.write(join([pathstr,"/presets.csv"]), dict)
    CSV.write(join([pathstr,"/$(figname).csv"]), data; delim = '\t')
    println(join([pathstr,"/$(figname).csv"]))
end 

fig2(datapath::String; kw...) = 
    cairofig1_v4(CSV.read(datapath, DataFrame,  delim = '\t'); kw...)
function fig2(data; ylims = (-0.015, 0.015))
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
                     1000 .* data.ϵs[selection][datak.<=-kcros+.001], 
                    color = ifelse(j != 1000 && 
                    intersect(i,(-1+Int(numbands/2):2+Int(numbands/2)))!=[],
                    :red, customgray),fxaa = true)

                lines!(ax, data.karray[selection][datak.>=kcros-.001],
                    1000 .* data.ϵs[selection][datak.>=kcros-.001], 
                    color = ifelse(j != 1000 && 
                    intersect(i,(-1+Int(numbands/2):2+Int(numbands/2)))!=[],
                    :red, customgray),fxaa = true)
                
                lines!(ax, datak[abs.(datak) .< kcros], 
                    1000 .* (data.ϵs[selection])[abs.(datak).<kcros], 
                    color = ifelse(j != 1000 && 
                    intersect(i,(-3+Int(numbands/2): 4+Int(numbands/2)))!=[],
                    :blue, customgray),fxaa = true)
                auxvect=collect(0:-2+Int(numbands/2))
                append!(auxvect, collect(3+Int(numbands/2):numbands))
                lines!(ax, datak[datak .> kcros], 
                    1000 .* (data.ϵs[selection])[datak.>kcros], 
                    color = ifelse(intersect(i, auxvect)!=[], 
                    (customgray, 1), (customgray, 0)), fxaa = true)
                lines!(ax, datak[datak .< -kcros],
                    1000 .* (data.ϵs[selection])[datak.<-kcros], 
                    color = ifelse(intersect(i, auxvect)!=[], 
                    (customgray, 1), (customgray, 0)), fxaa = true)           
            else
                if j == 2
                    kcrosr = 0.565
                    kcrosl = 0.435
                    ϵ = -0.01
                else
                    kcrosr = 0
                    kcrosl = 0
                    ϵ = 0
                end

                datak = range(0.2, stop = 0.8, length = length(data.ϵs[selection]))
                lines!(ax, datak[datak.<= kcrosl], 
                    1000 .* data.ϵs[selection][datak.<= kcrosl], color = ifelse( 
                    intersect(i,(-1+Int(numbands/2):2+Int(numbands/2)))!=[],
                    :red, (customgray,0)), fxaa = true)
                lines!(ax, datak[datak.>= kcrosr], 
                    1000 .* data.ϵs[selection][datak.>= kcrosr], color = ifelse( 
                    intersect(i,(-1+Int(numbands/2):2+Int(numbands/2)))!=[],
                    :red, (customgray,0)), fxaa = true)
                lines!(ax, datak[kcrosl+ ϵ .<= datak .<= kcrosr- ϵ], 
                    1000 .* data.ϵs[selection][kcrosl + ϵ .<=datak.<= kcrosr- ϵ],
                    color = ifelse(intersect(i,(-3+Int(numbands/2):4+Int(numbands/2)))!=[],
                        :blue, customgray), fxaa = true)

                auxvect = [-3,-2,3,4] .+ Int(numbands/2)
                inter = intersect(i, auxvect)
                lines!(ax, datak[datak .>= kcrosr ], 
                    1000 .* (data.ϵs[selection])[datak.>=kcrosr ], 
                    color = ifelse(inter!=[], (:blue, 1), 
                        ifelse(findall(x -> x == i, [-1,0,1,2] .+ Int(numbands/2)) == [], 
                        (customgray, 1), (:red,1))), fxaa = true)
                lines!(ax, datak[datak .<= kcrosl],
                    1000 .* (data.ϵs[selection])[datak .<= kcrosl], 
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
