function spectrumsweep(p, μlist, selfy = false)
    pn = reconstruct(p, μN = μlist[1])
    h = rectangle_randombounds_sc(pn, 0, 0.0, sidecontacts = true, selfy = selfy)
    sp = spectrum(h, method = ArpackPackage(nev = 16, sigma = 0.0im))
    energieslist = zeros(Float64, length(sp.energies), length(μlist))
    energieslist[:,1] = sp.energies
    for i in 2:length(μlist)
        println(i/length(μlist))
        pn = reconstruct(p, μN = μlist[i])
        h = rectangle_randombounds_sc(pn, 0, 0.0, sidecontacts = true, selfy = selfy)
        sp = spectrum(h,
            method = ArpackPackage(nev = 16, sigma = 0.0im))
        energieslist[:,i] = sp.energies
    end
    return energieslist, μlist
end

function spectrumsweepλ(p, λlist, selfy = false)
    pn = reconstruct(p, λ = λlist[1])
    h = rectangle_randombounds_sc(pn, 0, 0.0, sidecontacts = true, selfy = selfy)
    sp = spectrum(h, method = ArpackPackage(nev = 16, sigma = 0.0im))
    energieslist = zeros(Float64, length(sp.energies), length(λlist))
    energieslist[:,1] = sp.energies
    for i in 2:length(λlist)
        println(i/length(λlist))
        pn = reconstruct(p, λ = λlist[i])
        h = rectangle_randombounds_sc(pn, 0, 0.0, sidecontacts = true)
        sp = spectrum(h,
            method = ArpackPackage(nev = 16, sigma = 0.0im))
        energieslist[:,i] = sp.energies
    end
    return energieslist, λlist
end

function spectrumsweepb(p, blist, selfy = false)
    pn = reconstruct(p, EZ = SA[0, blist[1], 0])
    h = rectangle_randombounds_sc(pn, 0., 0.0, sidecontacts = true, selfy = selfy)
    sp = spectrum(h, method = ArpackPackage(nev = 32, sigma = 0.0im))
    energieslist = zeros(Float64, length(sp.energies), length(blist))
    energieslist[:,1] = sp.energies
    for i in 2:length(blist)
        println(i/length(blist))
        pn = reconstruct(p, EZ = SA[0, blist[i], 0])
        h = rectangle_randombounds_sc(pn, 0, 0.0, sidecontacts = true, selfy = selfy)
        sp = spectrum(h,
            method = ArpackPackage(nev = 32, sigma = 0.0im))
        energieslist[:,i] = sp.energies
    end
    return energieslist, blist
end

function spectrumsweepα(p, list, selfy = false)
    pn = reconstruct(p, α = list[1])
    h = rectangle_randombounds_sc(pn, 0, 0.0, sidecontacts = true, selfy = selfy)
    @time sp = spectrum(h, method = ArpackPackage(nev = 32, sigma = 0.0im))
    energieslist = SharedArray(zeros(Float64, length(sp.energies), length(list)))
    energieslist[:,1] = sp.energies
    for i in 2:length(list) #@sync @distributed 
        println(i/length(list))
        pn = reconstruct(p, α = list[i])
        h = rectangle_randombounds_sc(pn, 0, 0.0, sidecontacts = true, selfy = selfy)
        sp = spectrum(h,
            method = ArpackPackage(nev = 32, sigma = 0.0im))
        energieslist[:,i] = sp.energies
    end
    return energieslist, list
end

"given somre presets and a list of angles (in degrees) computes the spectrum of the 
lowest energy states as a funcion of θ"
function splittingvsrotation(p, list, selfy = true; kw...)
    list = list .* π/180
    println("first h...")
    h = rectangle_randombounds_sc(p, list[1], 0.0, sidecontacts = true, selfy = selfy; kw...)
    println("first sp...")
    sp = spectrum(h, method = ArpackPackage(nev = 16, sigma = 0.0im))
    energieslist = SharedArray(zeros(Float64, length(sp.energies), length(list)))
    energieslist[:,1] = sp.energies
    @sync @distributed for i in 2:length(list) 
        println(i/length(list))
        println("h ...")
        hn = rectangle_randombounds_sc(p, list[i], 0.0, sidecontacts = true, selfy = selfy; kw...)
        println("sp...")
        spn = spectrum(hn,
            method = ArpackPackage(nev = 16, sigma = 0.0im))
        energieslist[:,i] = spn.energies
    end
    return energieslist, list
end



function splittingvsdisorder(p, list, selfy = true; kw...)
    list = list .* π/180
    println("first h...")
    h = rectangle_randombounds_sc(p, 0.0, list[1], sidecontacts = true, selfy = selfy; kw...)
    println("first sp...")
    sp = spectrum(h, method = ArpackPackage(nev = 16, sigma = 0.0im))
    energieslist = SharedArray(zeros(Float64, length(sp.energies), length(list)))
    energieslist[:,1] = sp.energies
    @sync @distributed for i in 2:length(list) 
        println(i/length(list))
        println("list[i] ...")
        hn = rectangle_randombounds_sc(p, 0.0, list[i], sidecontacts = true, selfy = selfy; kw...)
        println("sp...")
        spn = spectrum(hn,
            method = ArpackPackage(nev = 16, sigma = 0.0im))
        energieslist[:,i] = spn.energies
    end
    return energieslist, list
end

# function merge_max2(a, b)
#     x, y = Baselet.sort((a..., b...); by = last, rev = true)
#     return (x, y)
# end
