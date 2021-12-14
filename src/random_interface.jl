using Quantica, LinearAlgebra, StatsBase
using Baselet, Arpack
# include("model.jl")

############################################################################################
# Random model
############################################################################################

"builds the hamiltonian of BLG with a random model implemented 
on the edges of the system see `rectangle`, the degree of disorder
is parametrized by a value `η` which determines the average number
of atom vacancies in the edges. 
    `η = 0` -> crystalographic edges
    `η = 0.` is the default    
"
function rectangle_randombounds(p, θ = 0, η = 0.; mono = false, sidecontacts = false)
    (; Ln, Ls, Δ, a0, τ, d) = p
    (; model0, field!) = modelS(p)
    lat0 = mono ? latSLG(p) : latBLG(p)
    Quantica.transform!(r -> SA[cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1] * r, lat0)

    myreg = RP.rectangle( (Ln , Ln))
    lat = unitcell(lat0, region = myreg)
    regionS(r) = abs(r[1]) >= Ln/2 || abs(r[2]) >= Ln/2
    regionS(r, dr) = regionS(r)
    
    smooth_sides(r) = iszero(d) ? ifelse(abs(r[1])<Ln/2 && abs(r[2])<Ln/2, 0.0, 1.0) :
            1 - (1 + tanh((Ln/2-abs(r[1]))/d))/2

    smooth(r) = iszero(d) ? ifelse(abs(r[1])<Ln/2 && abs(r[2])<Ln/2, 0.0, 1.0) :
                1 - (1 + tanh((Ln/2-abs(r[1]))/d))*(1 + tanh((Ln/2-abs(r[2]))/d))/4

    smooth_method(r) = ifelse(sidecontacts == true , smooth_sides(r), smooth(r))
    
    #random model (we are randomly replacing sites by vacancies) within rand_region
    rand_region(r) = Ln/2 - 2a0 <= abs(r[1]) <= Ln/2 || Ln/2 - 2a0 <= abs(r[2]) <= Ln/2
    random_mat() = 
        @SMatrix[sample([0, 1], Weights([1-η, η])) 0 0 0;
            0 sample([0, 1], Weights([1-η, η])) 0 0;
            0 0 -sample([0, 1], Weights([1-η, η])) 0; 
            0 0 0 -sample([0, 1], Weights([1-η, η]))]
    so_rand! = @onsite!((o, r) -> o + rand_region(r) * ifelse(η != 0, 1e3, 0) * 
        random_mat())

    #building the hamiltonian
    ph = lat |> hamiltonian(model0; orbitals = Val(4)) |> parametric(field!, so_rand!)
    h = ph()
    return h
end

"builds the hamiltonian of BLG with a random model implemented 
on the edges of the system see `rectangle` we also add two possible SC configurations
(square mask not yet implemented),
the degree of disorder is parametrized by a value `η` which determines the average number
of atom vacancies in the edges. 
    `η = 0` -> crystalographic edges
    `η = 0.` is the default    
"
function rectangle_randombounds_sc(p, θ = 0, η = 0.; mono = false, sidecontacts = false, selfy = false)
    println(selfy)
    (; Ln, Ls, Δ, a0, τ, d) = p
    (; model0, field!) = modelS(p)
    lat0 = mono ? latSLG(p) : latBLG(p)
    Quantica.transform!(r -> SA[cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1] * r, lat0)
    # contact model (Josephson vs rectangular mask architectures)
    W = Ln
    #Ln = 1000
    if sidecontacts == true
        myreg = RP.rectangle( (Ln + 2Ls, W))
    else
        myreg = RP.rectangle( (Ln + 2Ls, Ln + 2Ls))
    end

    # lat = unitcell(lat0, region = myreg)
    lat = unitcell(lat0, region = (r) -> r[1]> -Ln/2 && r[1]< Ln/2-0*1/2*a0 && r[2] < W/2 && r[2]> -W/2)
    
    regionS(r) = abs(r[1]) >= Ln/2 || abs(r[2]) >= W/2
    regionS(r, dr) = regionS(r)
    
    smooth_sides(r) = iszero(d) ? ifelse(abs(r[1])<Ln/2 && abs(r[2])<W/2, 0.0, 1.0) :
            1 - (1 + tanh((Ln/2-abs(r[1]))/d))/2

    smooth(r) = iszero(d) ? ifelse(abs(r[1])<Ln/2 && abs(r[2])<W/2, 0.0, 1.0) :
                1 - (1 + tanh((Ln/2-abs(r[1]))/d))*(1 + tanh((W/2-abs(r[2]))/d))/4

    smooth_method(r) = ifelse(sidecontacts == true , smooth_sides(r), smooth(r))
    
    #random model (we are randomly replacing sites by vacancies) within rand_region
    rand_region(r) = Ln/2 - 2a0 <= abs(r[1]) <= Ln/2 #&& Ln/2 - 2a0 <= abs(r[2]) <= Ln/2
    random_mat() = 
        sample([0, 1], StatsBase.Weights([1-η, η])) * σ0τz
    so_rand! = @onsite!((o, r) -> o + rand_region(r) * ifelse(η != 0, 1e3, 0) *
        random_mat())
    # local self energy model
    self_region(r) = ifelse(selfy == false, abs(r[1]) ≥ (round(0.5*Ln/a0)-0.9)*a0,  
        abs(r[1]) ≥ (round(0.5*Ln/a0)-0.9)*a0|| abs(r[2]) ≥ (round(0.5*W/a0)-0.9)*a0)

    # phase = 
    # Superconductivity
    SCo! = @onsite!((o, r) -> o + smooth_method(r) * Δ * σyτy)
    sCself! = @onsite!((o, r) -> o + self_region(r) * Δ * σyτy)
    SCh! = @hopping!((t, r, dr) -> t * (1-smooth_method(r));
        sublats = (:At, :Bt ,:Ab, :Bb, :A, :B) .=> (:At, :Bt ,:Ab, :Bb, :A, :B))
    SCτ! = @hopping!((t, r, dr) -> t * τ; region = (r, dr) -> xor(regionS(r - dr/2),
        regionS(r + dr/2)))
    #building the hamiltonian
    ph = lat |> hamiltonian(model0; orbitals = Val(4)) |> 
        parametric(field!, so_rand!, SCo!, sCself!, SCh!, SCτ!)
    return ph()
end

"""
SQUID geometry. Self energy on the zigzag edge with a phase difference (-ϕ/2 for x<0 ϕ/2 for x>0)
Neither smoothness nor leads are implemented. Also disorder and rotation are not considered.
"""
function rectangle_squid(p; W = 1440)
    (; Ln, Ls, Δ, a0, τ, d) = p
    (; model0, field!) = modelS(p)
    lat0 = latBLG(p)
    lat = unitcell(lat0, 
        region = (r) -> r[1]> -Ln/2 && r[1]< Ln/2-0*1/2*a0 && r[2] < W/2 && r[2]> -W/2)
    
    # SC
    self_region(r) =  
        abs(r[1]) ≥ (round(0.5*Ln/a0)-0.9)*a0|| abs(r[2]) ≥ (round(0.5*W/a0)-0.9)*a0

    diagphi(ϕ) = Diagonal(SA[cis(ϕ), cis(ϕ), cis(-ϕ), cis(-ϕ)])
    sCself! = @onsite!((o, r; ϕ) -> o + 
        self_region(r) * Δ * diagphi(sign(r[1])* ϕ/4) * σyτy * diagphi(sign(r[1])*ϕ/4)')     
    
    #return the parametric hamiltonian in ϕ
    return lat |> hamiltonian(model0; orbitals = Val(4))|> parametric(sCself!)
end

"""
CPR
"""
function cpr(ϕlist, p; nev = 16, kw...)
    ph = rectangle_squid(p; kw...)
    elist = zeros(Float64, nev, length(ϕlist))
    for i in 1:length(ϕlist)
        s = spectrum(ph(ϕ = ϕlist[i]), method = ArpackPackage(sigma = 1e-7, nev = nev))
        elist[:,i] = s.energies
    end
    return ϕlist, elist
end


############################################################################################
# LDOS
############################################################################################

function ldosonlattice_averaged(p, η, whichpanel = "a", iterations = 0; kw...) 
    psi, h0 = ldosonlattice_rand(p, η, whichpanel; kw...)
    [psi .+= ldosonlattice_rand(p, η, whichpanel; kw...)[1] for i in 1:iterations]
    psi ./= iterations +1
    return ldosonlattice(psi, h0)
end

struct Fig4_presets
    η::Float64
    angle::Float64
    iterations::Int64
    nummodes::Int64
    selfy::Bool
end

function ldosonlattice_averaged_sc(p, fig_p::Fig4_presets; kw...) 
    ldosonlattice_averaged_sc(p, fig_p.η, fig_p.angle,
        fig_p.iterations; nummodes = fig_p.nummodes, selfy = fig_p.selfy, kw...) 
end

function ldosonlattice_averaged_sc(p, η, angle = 0, iterations = 5; nummodes = 4, selfy = false, kw...) 
    println(nummodes)
    psi1, psi2 = ldosonlattice_rand_sc(p, η, angle, nummodes = nummodes, selfy = selfy; kw...)
    for i in 1:iterations
         psis = ldosonlattice_rand_sc(p, η, angle, nummodes = nummodes, selfy = selfy; kw...)
         psi1 .+= psis[1]
         psi2 .+= psis[2]
    end
    # psi ./= iterations +1
    return psi1, psi2#, h0#ldosonlattice(psi, h0)
end

"KPM calculation of localdensity of states at a site ind for the bounded system"
function ldos(params, ind, η, angle = 0.0; order = 300000, bandrange = missing, kw...)
    sidecontacts = true
    h = rectangle_randombounds_sc(params, angle, η, sidecontacts = sidecontacts; kw...);
    si = SMatrix{4,4}(Diagonal(SA[1,1,-1,-1]))
    km = ketmodel(si, indices = ind)
    dos = dosKPM(h, ket = km, bandrange = bandrange; order = order)
    return dos
end


"returns the wavefunctions of nev eigenvalues (you must select them to be inside the 
helical gap manually and the h0 for plotting purposes for a random model on the edges
    see: rectangle_randombounds∂"
function ldosonlattice_rand_norm(p, η, whichpanel = "a"; kw...)
    p_0 = reconstruct(p, λ = 0, α = 0)
    h = rectangle_randombounds(p, angle, η; kw...);
    h0 = rectangle_randombounds(p_0, angle, η; kw...);
    s = spectrum(h, method = ArpackPackage(sigma = 1e-7, nev = 64))
    println(s.energies)
    psi = s[1].basis
    [psi .+= s[i].basis for i in 2:ceil(length(s.energies)÷2)]
    return psi, h0
end

function ldosonlattice_rand_sc(p, η, angle = 0; nummodes = 4, selfy = false, kw...)
    p_0 = reconstruct(p, λ = 0, α = 0)
    sidecontacts = true
    h = rectangle_randombounds_sc(p, angle, η, sidecontacts = sidecontacts, selfy = selfy; kw...);
    s = spectrum(h, method = ArpackPackage(sigma = 1e-7, nev = ifelse(nummodes == 4, 4, 8))) #nev = 64 for normal plots
    println(s.energies)
    if nummodes == 4
        psi =  Vector{Float64}(undef, Int(size(h,1)*4))
        psi .=  abs.(s.states[:,1])
        psi .+= abs.(s.states[:,2])
        psi .+= abs.(s.states[:,3])
        psi .+= abs.(s.states[:,4])
        return psi, psi
    else nummodes == 8
        psi1 =  Vector{Float64}(undef, Int(size(h,1)*4))
        psi2 =  Vector{Float64}(undef, Int(size(h,1)*4))
        psi1 .=  abs.(s.states[:,1])
        psi1 .+= abs.(s.states[:,2])
        psi1 .+= abs.(s.states[:,7])
        psi1 .+= abs.(s.states[:,8])
        psi2 .=  abs.(s.states[:,3])
        psi2 .+= abs.(s.states[:,4])
        psi2 .+= abs.(s.states[:,5])
        psi2 .+= abs.(s.states[:,6])
        return psi1, psi2
    end

end

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

function merge_max2(a, b)
    x, y = Baselet.sort((a..., b...); by = last, rev = true)
    return (x, y)
end

# ldosonlattice(psi, h0) = vlplot(h0, psi, 
#     sitesize = DensityShader(), sitecolor = DensityShader(), siteopacity = 0.5,
#     mindiameter = 2e-2, size = (300), colorscheme = "blues", sitestroke = nothing,
#     maxthickness = 0.005, plotlinks = false, discretecolorscheme = :gray, 
#     labels = ("x (nm)", "y (nm)")) 

ldosonlattice(psi, h0, scheme = "blues") = vlplot(h0, psi, 
    sitesize = DensityShader(), sitecolor = DensityShader(), siteopacity = 0.5,
    mindiameter = 5e-2, size = (300), colorscheme = scheme, sitestroke = nothing, #reds
    maxdiameter =30, plotlinks = false, discretecolorscheme = :gray, 
    labels = ("x (nm)", "y (nm)")) 


logldosonlattice(psi, h0) = vlplot(h0, -log.(abs.(psi)), 
    sitesize = DensityShader(), sitecolor = DensityShader(), siteopacity = 0.5,
    mindiameter = maximum(-log.(abs.(psi)))*4/5, size = (300), colorscheme = "blues", sitestroke = nothing,
    maxdiameter =maximum(-log.(abs.(psi))), plotlinks = false, discretecolorscheme = :gray, 
    labels = ("x (nm)", "y (nm)")) 

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
