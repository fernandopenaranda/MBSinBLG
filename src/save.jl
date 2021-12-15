############################################################################################
# Saving functions to store the output of figXrun() where X ∈ {2, 3, 4, 5, 6, 7}

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

"save data"
function savespectrum(figname, p::Params, sp::Quantica.Spectrum)
    time_str = string(now())
    isdir("data") ? nothing : mkdir(data)
    isdir("data/$(figname)") ? nothing : mkdir("data/$(figname)")
    pathstr = join([pwd(), "/data/$(figname)/", time_str])
    mkdir(pathstr) 
    CSV.write(join([pathstr,"/presets.csv"]), DataFrame(type2dict(p)))
    psi_data = DataFrame(psi = sp.states[:,1])
    CSV.write(join([pathstr,"/psi.csv"]), psi_data; delim = '\t')
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
    #check and/or create folder structure
    time_str = string(now())
    isdir("data") ? nothing : mkdir(data)
    isdir("data/$(figname)") ? nothing : mkdir("data/$(figname)")
    pathstr = join([pwd(), "/data/$(figname)/", time_str])
    mkdir(pathstr) 
    #build the csv.files
    for i in 1:length(presetsarray)
        CSV.write(join([pathstr,"/presets","$(i)",".csv"]),
            DataFrame(type2dict(presetsarray[i])))
    end
    CSV.write(join([pathstr,"/$(figname).csv"]), data; delim = '\t')
    println(join([pathstr,"/$(figname).csv"]))
end 
#


#read
function readfig3data(datapath::String, angle = 0*pi/180)
    println("hey")
    presets = CSV.read(joinpath(datapath,"presets.csv"), DataFrame,  delim = '\t')
    println(presets)
    data = CSV.read(joinpath(datapath,"psi.csv"), DataFrame, delim='\t')
    return  map(x->parse(ComplexF64,x), data.psi)
end

function readfig3psi(datapath::String)
    presets = CSV.read(joinpath(datapath,"presets.csv"), DataFrame,  delim = '\t')
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

