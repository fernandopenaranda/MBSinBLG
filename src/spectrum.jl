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
    h = rectangle_randombounds_sc(pn, 0, 0.0, sidecontacts = true, selfy = selfy)
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
function splittingvsrotation(p, list, selfy = true)
    list = list .* π/180
    h = rectangle_randombounds_sc(p, list[1], 0.0, sidecontacts = true, selfy = selfy)
    sp = spectrum(h, method = ArpackPackage(nev = 32, sigma = 0.0im))
    energieslist = zeros(Float64, length(sp.energies), length(list))
    energieslist[:,1] = sp.energies
    for i in 2:length(list) #@sync @distributed 
        println(i/length(list))
        h = rectangle_randombounds_sc(p, list[i], 0.0, sidecontacts = true, selfy = selfy)
        sp = spectrum(h,
            method = ArpackPackage(nev = 32, sigma = 0.0im))
        energieslist[:,i] = sp.energies
    end
    return energieslist, list
end

# function merge_max2(a, b)
#     x, y = Baselet.sort((a..., b...); by = last, rev = true)
#     return (x, y)
# end

########
# Save
######

"save data"
function savepsi(figname, p::Params, psi)
    time_str = string(now())
    isdir("data") ? nothing : mkdir(data)
    isdir("data/$(figname)") ? nothing : mkdir("data/$(figname)")
    pathstr = join([pwd(), "/data/$(figname)/", time_str])
    mkdir(pathstr) 
    CSV.write(join([pathstr,"/presets.csv"]), DataFrame(type2dict(p)))
    psi_data = DataFrame(psi = psi)
    CSV.write(join([pathstr,"/psi.csv"]), psi_data; delim = '\t')
    return pathstr
end

function savespectrum(figname, p::Params, sp1, sp2)
    time_str = string(now())
    isdir("data") ? nothing : mkdir("data")
    isdir("data/$(figname)") ? nothing : mkdir("data/$(figname)")
    pathstr = join([pwd(), "/data/$(figname)/", time_str])
    mkdir(pathstr) 
    CSV.write(join([pathstr,"/presets.csv"]), DataFrame(type2dict(p)))
    CSV.write(join([pathstr,"/spectrum.csv"]),  DataFrame(sp1, :auto))
    psi_data = DataFrame(EZ = sp2)
    CSV.write(join([pathstr,"/EZ.csv"]), psi_data; delim = '\t')
    return pathstr
end

function savecprs(figname, p::Params, phi, ea, eb, ec, ed)
    time_str = string(now())
    isdir("data") ? nothing : mkdir("data")
    isdir("data/$(figname)") ? nothing : mkdir("data/$(figname)")
    pathstr = join([pwd(), "/data/$(figname)/", time_str])
    mkdir(pathstr)
    CSV.write(join([pathstr,"/presets.csv"]), DataFrame(type2dict(p)))
    phi_data = DataFrame(phi = collect(phi))
    CSV.write(join([pathstr,"/phi.csv"]), phi_data; delim = '\t')
    CSV.write(join([pathstr,"/cpra.csv"]),  DataFrame(ea, :auto))
    CSV.write(join([pathstr,"/cprb.csv"]),  DataFrame(eb, :auto))
    CSV.write(join([pathstr,"/cprc.csv"]),  DataFrame(ec, :auto))
    CSV.write(join([pathstr,"/cprd.csv"]),  DataFrame(ed, :auto))
    return pathstr
end
#read

function readfig3psi(datapath::String)
    presets = CSV.read(joinpath(datapath,"presets.csv"), DataFrame,  delim = '\t')
    data = CSV.read(joinpath(datapath,"psi.csv"), DataFrame, delim='\t')
    return  data.psi#map(x->parse(Float64,x), data.psi)
end


function readspectrum(datapath::String)
    y = CSV.read(join([datapath,"/spectrum.csv"]), DataFrame)
    x = CSV.read(join([datapath,"/EZ.csv"]), DataFrame)
    return x, y
end


function readcpr(datapath::String)
    phi = CSV.read(join([datapath,"/phi.csv"]), DataFrame)
    cpra = CSV.read(join([datapath,"/cpra.csv"]), DataFrame)
    cprb = CSV.read(join([datapath,"/cprb.csv"]), DataFrame)
    cprc = CSV.read(join([datapath,"/cprc.csv"]), DataFrame)
    cprd = CSV.read(join([datapath,"/cprd.csv"]), DataFrame)
    return phi, cpra, cprb, cprc, cprd
end
