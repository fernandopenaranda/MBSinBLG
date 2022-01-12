############################################################################################
# Functions to write/read the output of figXrun() where X ∈ {2, 3, 4, 5, 6, 7} as .csv files
############################################################################################

# using CSV: write, read

#######
# WRITE
#######

function savepsi(figname, p::Params, psi)
    path_str = filesstorage_settings(figname)
    CSV.write(join([path_str,"/presets.csv"]), DataFrame(type2dict(p)))
    psi_data = DataFrame(psi = psi)
    CSV.write(join([path_str,"/psi.csv"]), psi_data; delim = '\t')
    return path_str
end

function savespectrum(figname, p::Params, sp1, sp2)
    path_str = filesstorage_settings(figname)
    CSV.write(join([path_str,"/presets.csv"]), DataFrame(type2dict(p)))
    CSV.write(join([path_str,"/spectrum.csv"]),  DataFrame(sp1, :auto))
    psi_data = DataFrame(EZ = sp2)
    CSV.write(join([path_str,"/EZ.csv"]), psi_data; delim = '\t')
    return path_str
end

"save data"
function savespectrum(figname, p::Params, sp::Quantica.Spectrum)
    path_str = filesstorage_settings(figname)
    CSV.write(join([path_str,"/presets.csv"]), DataFrame(type2dict(p)))
    psi_data = DataFrame(psi = sp.states[:,1])
    CSV.write(join([path_str,"/psi.csv"]), psi_data; delim = '\t')
    return path_str
end

function savespectrumvsphase(figname, p::Params, phi, ea, eb, ec, ed)
    path_str = filesstorage_settings(figname)
    CSV.write(join([path_str,"/presets.csv"]), DataFrame(type2dict(p)))
    phi_data = DataFrame(phi = collect(phi))
    CSV.write(join([path_str,"/phi.csv"]), phi_data; delim = '\t')
    CSV.write(join([path_str,"/cpra.csv"]),  DataFrame(ea, :auto))
    CSV.write(join([path_str,"/cprb.csv"]),  DataFrame(eb, :auto))
    CSV.write(join([path_str,"/cprc.csv"]),  DataFrame(ec, :auto))
    CSV.write(join([path_str,"/cprd.csv"]),  DataFrame(ed, :auto))
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
    path_str = filesstorage_settings(figname)
    for i in 1:length(presetsarray)
        CSV.write(join([path_str,"/presets","$(i)",".csv"]),
            DataFrame(type2dict(presetsarray[i])))
    end
    CSV.write(join([path_str,"/$(figname).csv"]), data; delim = '\t')
    println(join([path_str,"/$(figname).csv"]))
end 

#######
# READ
#######
function readfig4data(datapath::String, angle = 0*pi/180)
    presets = CSV.read(joinpath(datapath,"presets.csv"), DataFrame,  delim = '\t')
    println(presets)
    data = CSV.read(joinpath(datapath,"psi.csv"), DataFrame, delim='\t')
    return  map(x->parse(ComplexF32,x), data.psi)
end

function readfig4psi(datapath::String)
    data = CSV.read(joinpath(datapath,"psi.csv"), DataFrame, delim='\t')
    return  data.psi
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

