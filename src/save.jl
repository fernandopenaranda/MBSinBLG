############################################################################################
# Saving functions to store the output of figXrun() where X ∈ {2, 3, 4, 5, 6, 7}

using CSV: write, read

"save data"
function savepsi(figname, p::Params, psi)
    path_str = filesstorage_settings(figname)
    write(join([path_str,"/presets.csv"]), DataFrame(type2dict(p)))
    psi_data = DataFrame(psi = psi)
    write(join([path_str,"/psi.csv"]), psi_data; delim = '\t')
    return path_str
end

function savespectrum(figname, p::Params, sp1, sp2)
    path_str = filesstorage_settings(figname)
    write(join([path_str,"/presets.csv"]), DataFrame(type2dict(p)))
    write(join([path_str,"/spectrum.csv"]),  DataFrame(sp1, :auto))
    psi_data = DataFrame(EZ = sp2)
    write(join([path_str,"/EZ.csv"]), psi_data; delim = '\t')
    return path_str
end

"save data"
function savespectrum(figname, p::Params, sp::Quantica.Spectrum)
    path_str = filesstorage_settings(figname)
    write(join([path_str,"/presets.csv"]), DataFrame(type2dict(p)))
    psi_data = DataFrame(psi = sp.states[:,1])
    write(join([path_str,"/psi.csv"]), psi_data; delim = '\t')
    return path_str
end

function savecprs(figname, p::Params, phi, ea, eb, ec, ed)
    path_str = filesstorage_settings(figname)
    write(join([path_str,"/presets.csv"]), DataFrame(type2dict(p)))
    phi_data = DataFrame(phi = collect(phi))
    write(join([path_str,"/phi.csv"]), phi_data; delim = '\t')
    write(join([path_str,"/cpra.csv"]),  DataFrame(ea, :auto))
    write(join([path_str,"/cprb.csv"]),  DataFrame(eb, :auto))
    write(join([path_str,"/cprc.csv"]),  DataFrame(ec, :auto))
    write(join([path_str,"/cprd.csv"]),  DataFrame(ed, :auto))
    return path_str
end

""" check and/or create folder structure"""
function filesstorage_settings(figname)
    time_str = string(now())
    isdir("data") ? nothing : mkdir("data")
    isdir("data/$(figname)") ? nothing : mkdir("data/$(figname)")
    path_str = join([pwd(), "/data/$(figname)/", time_str])
    mkdir(path_str)
    return path_str
end
#fig2

floattostring(vec::Array{Float64,1}) = [floattostring(vec[i]) for i in 1:length(vec)]
floattostring(i::Float64) = "$(Int(i))"
function writecsv(vec, p, name, taglist)
    karray = collect(Iterators.flatten(
    [collect(Iterators.flatten([append!([],vec[1][1]) 
        for i in 1:size(vec[1][2],1)])) for j in 1:size(vec,1)]))
    ϵs = real(collect(Iterators.flatten(
        [collect(Iterators.flatten([append!([], vec[j][2][i,:]) 
            for i in 1:size(vec[j][2], 1)])) for j in 1:size(vec,1)])))
    psi = collect(Iterators.flatten(
        [floattostring([ceil(i/size((vec[1][2]),2)) 
            for i in 1:length(vec[1][2])]) for j in 1:size(vec,1)]))
    tag = [taglist[Int(ceil(i/length(vec[1][2])))] 
        for i in 1:size(vec,1) * length(vec[1][2])]
    datafig = DataFrame(tag = tag, karray = karray, ϵs = ϵs, psi = psi)
    savecsv(name, p, datafig)
end

function savecsv(figname, presetsarray::AbstractArray, data::DataFrame) 
    path_str = filesstorage_settings()
    for i in 1:length(presetsarray)
        write(join([path_str,"/presets","$(i)",".csv"]),
            DataFrame(type2dict(presetsarray[i])))
    end
    write(join([path_str,"/$(figname).csv"]), data; delim = '\t')
    println(join([path_str,"/$(figname).csv"]))
end 


#read
function readfig4data(datapath::String, angle = 0*pi/180)
    println("hey")
    presets = CSV.read(joinpath(datapath,"presets.csv"), DataFrame,  delim = '\t')
    println(presets)
    data = CSV.read(joinpath(datapath,"psi.csv"), DataFrame, delim='\t')
    return  map(x->parse(ComplexF32,x), data.psi)
end

function readfig4psi(datapath::String)
    presets = CSV.read(joinpath(datapath,"presets.csv"), DataFrame,  delim = '\t')
    data = CSV.read(joinpath(datapath,"psi.csv"), DataFrame, delim='\t')
    return  data.psi
end

function readspectrum(datapath::String)
    y = CSV.read(join([datapath,"/spectrum.csv"]), DataFrame)
    println("hey")
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

