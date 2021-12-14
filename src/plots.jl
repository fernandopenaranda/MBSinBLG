using LinearAlgebra, Arpack,  Dates, CSV, DataFrames, VegaLite, FFTW, MKLSparse
include("model.jl")
include("Plot_functions.jl")

#________________________________FUNCTIONS_FIGURE_3_________________________________________
function ldosdata_gen()
    fig3data(pa, "a")
    fig3data(pb, "b")
    fig3data(pc, "c")
    fig3data(pc, "d")
end

"Density of states for the lowest energy state on top of the lattice"
function fig3panels(presets_mask, whichpanel = "a"; load = true)
    sidecontacts = (ifelse(whichpanel == "a" || whichpanel == "b", true, false))
    angle = (ifelse(whichpanel == "d", 1*pi/180, 0.))
    presets_mask_0 = reconstruct(presets_mask, λ = 0, α = 0)
    h = rectangle(presets_mask, angle, sidecontacts = sidecontacts);
    h0 = rectangle(presets_mask_0, angle, sidecontacts = sidecontacts);
    s = spectrum(h, method = ArpackPackage(sigma = 1e-7, nev = 8))
    println(s.energies)
    println(s[2])
    if load == true 
        datapath = savespectrum("fig3", presets_mask, s)
    else nothing end
    return ldosonlattice(s[2], h0)
end

"KPM calculation of localdensity of states at a site ind for the bounded system"
function ldos(params, whichpanel, ind, angle = 0; order = 30000, bandrange = missing)
    sidecontacts = (ifelse(whichpanel == "a" || whichpanel == "b", true, false))
    h = rectangle(params, angle, sidecontacts =  sidecontacts);
    si = SMatrix{4,4}(Diagonal(SA[1,1,-1,-1]))
    km = ketmodel(si, indices = ind)
    dos = dosKPM(h, ket = km, bandrange = bandrange; order = order)
    return dos
end

"finds the position (index) where the maximum amplitude is to pass it to the 
KPM function for ldos vs w calculations"
function cornerfinder(presets_mask, whichpanel = "a"; load = true)
    sidecontacts = (ifelse(whichpanel == "a" || whichpanel == "b", true, false))
    angle = (ifelse(whichpanel == "d", 1*pi/180, 0.))
    presets_mask_0 = reconstruct(presets_mask, λ = 0, α = 0)
    h = rectangle(presets_mask, angle, sidecontacts = sidecontacts);
    h0 = rectangle(presets_mask_0, angle, sidecontacts = sidecontacts);
    s = spectrum(h, method = ArpackPackage(sigma = 1e-7, nev = 4))
    println(s.energies)
    if load == true 
        datapath = savespectrum("fig3", presets_mask, s)
    else nothing end
    return psi = abs.(s.states[:,1])
end

# function fig3data(presets_mask, whichpanel = "a"; load = true)
#     sidecontacts = (ifelse(whichpanel == "a" || whichpanel == "b", true, false))
#     angle = (ifelse(whichpanel == "d", 1*pi/180, 0.))
#     presets_mask_0 = reconstruct(presets_mask, λ = 0, α = 0)
#     h = rectangle(presets_mask, angle, sidecontacts = sidecontacts);
#     h0 = rectangle(presets_mask_0, angle, sidecontacts = sidecontacts);
#     s = spectrum(h, method = ArpackPackage(sigma = 1e-7, nev = 4))
#     println(s.energies)
#     if load == true 
#         datapath = savespectrum("fig3", presets_mask, s)
#     else nothing end
#     return nothing
# end
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

"import data"
function readfig3data(datapath::String, angle = 0*pi/180)
    println("hey")
    presets = CSV.read(joinpath(datapath,"presets.csv"), DataFrame,  delim = '\t')
    println(presets)
    data = CSV.read(joinpath(datapath,"psi.csv"), DataFrame, delim='\t')
    return  map(x->parse(ComplexF64,x), data.psi)
end

#________________________________COMPUTATION________________________________________________
## Optimized presets

palarge = Params(Ln = 800, Ls = 200, scale = 40, λ = 3, α = 2, EZ = SA[0, 0.1, 0],
                  E = [0.000000, 0, 0], μN = 0.85, Δ = 1.5, d = 1, τ = .8);

pblarge = Params(Ln = 800, Ls = 200, scale = 40, λ = 3, α = 2, EZ = SA[0, 1, 0],
    E = [0.000000, 0, 0], μN = 0.85, Δ = 1.5, d = 1, τ = .8);

pcdlarge = Params(Ln = 800, Ls = 200, scale = 40, λ = 3, α = 2, EZ = SA[0, 1, 0],
    E = [0.000000, 0, 0], μN = 0.85, Δ = 1.5, d = 1, τ = .1);

## indexes for the corners
indexa = 26199 #in the middle of the edge (not used)
indexab = 29959 
indexc = 44968
indexd = 41297

# Parameters
pcdtest = Params(Ln = 800, Ls = 100, scale = 40, λ = 3, α = 0, EZ = SA[0, 3, 0],
           E = [0.000000, 0, 0], μN = 0.65, Δ = 5, d = 0, τ = .2);

pcdtest = Params(Ln = 400, Ls = 50, scale = 40, λ = 3, α = 0, EZ = SA[0, 2, 0],
           E = [0.000000, 0, 0], μN = 0.65, Δ = 1, d = 1.5, τ = 1);

pbtest = Params(Ln = 200, Ls = 50, scale = 20, λ = 60, α = 0.0, d = 4, 
    τ = 0.55, EZ = SA[0, 12, 0], E = SA[0.0000, 0.0000, 0], μN = 20, 
    Δ = 20)


pbtest2 = Params(Ln = 800, Ls = 200, scale = 40, λ = 3, α = 2, EZ = SA[0, 2, 0],
    E = [0.000000, 0, 0], μN = 0.85, Δ = 1.5, d = 1, τ = 1);

# See normal bands
# bands= normalbands(pcdtest); bands[1]; bands[2]
# bands= scbands(pcdtest); bands[1]; bands[2]


#___________________________________________________________________________________________
## Computation and plot panels a-d

# fig = fig3panels(palarge, "a"; load = false)
# fig = fig3panels(pblarge, "b"; load = false)
# fig = fig3panels(pcdlarge, "c"; load = false)
# fig = fig3panels(pcdlarge, "d"; load = false)

#___________________________________________________________________________________________
## Computation panels e-f
# ldosa = ldos(palarge, "a", indexab, 0, order = 30000) 
# ldosb = ldos(pblarge, "b", indexab, 0, order = 30000)
# ldosc = ldos(pcdlarge, "c", indexc, 0, order = 30000)
# ldosd = ldos(pcdlarge, "d", indexd, pi/180, order = 30000)

# adata = DataFrame(x = ldosa[1], y = ldosa[2])
# bdata = DataFrame(x = ldosb[1], y = ldosb[2])
# cdata = DataFrame(x = ldosc[1], y = ldosc[2])
# ddata = DataFrame(x = ldosd[1], y = ldosd[2]) 

# CSV.write("data/fig3_ldosa", adata; delim = '\t')
# CSV.write("data/fig3_ldosb", bdata; delim = '\t')
# CSV.write("data/fig3_ldosc", cdata; delim = '\t')
# CSV.write("data/fig3_ldosd", ddata; delim = '\t')

## Plot panelrightfig3(q[1], dataa.y, datab.y, datac.y, datad.y)

#panelrightfig3(q[1], adata.y, bdata.y, cdata.y, ddata.y)

#___________________________________________________________________________________________
# Visualization tests

########
# Lattice plot
# pcdminimal = Params(Ln = 800, Ls = 200, scale = 200, λ = 0, α = 0, EZ = SA[0, 0, 0], 
#     E = [0.000000, 0, 0], μN = 2.7, Δ = 10, d = 4, τ = .5)
# h = rectangle(pcdminimal, 0/180, sidecontacts = false);

# lata = vlplot(h, size = (300), maxdiameter = .5, mindiameter = 0, maxthickness = 2, 
#   colorscheme = "greys", plotlinks = true, sitestroke = nothing, plotsites = true, 
#   linkcolor = -0.5)

# latb = vlplot(h, size = (300), maxdiameter = .5, mindiameter = 0, maxthickness = 2,
#   colorscheme = "lightgreyred", linkcolor = 2,  
#   plotlinks = true, sitestroke = nothing, plotsites = true, linkopacity = 0.1)
#######

