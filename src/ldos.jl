############################################################################################
# LDOS
############################################################################################

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