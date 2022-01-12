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

"
    ldosonlattice_averaged_sc(p, fig_p::Fig4_presets; kw...) 
average ldos weights of a number: Fig4_presets.nummodes of inner eigenvalues of a system 
built with parameters p and  a given rotation angle Fig4_presets.angle and disorder 
model Fig4_presets.angle.  The output is averaged overa a number of disorder realizations
specified by Fig4_presets.iterations.
    see: Fig4_presets, ldosonlattice_rand_sc, rectangle_randombounds_sc "

function ldosonlattice_averaged_sc(p, fig_p::Fig4_presets; kw...) 
    ldosonlattice_averaged_sc(p, fig_p.η, fig_p.angle,
        fig_p.iterations; nummodes = fig_p.nummodes, selfy = fig_p.selfy, kw...) 
end

function ldosonlattice_averaged_sc(p, η, angle = 0, iterations = 5; nummodes = 4, 
    selfy = false, kw...) 
    println(nummodes)
    psi1, psi2 = ldosonlattice_rand_sc(p, η, angle, nummodes = nummodes, selfy = selfy; kw...)
    for i in 1:iterations
        psis = ldosonlattice_rand_sc(p, η, angle, nummodes = nummodes, selfy = selfy; kw...)
        psi1 .+= psis[1]
        psi2 .+= psis[2]
    end
    psi ./= iterations +1
    return psi1, psi2
end

"KPM calculation of localdensity of states at a site ind for the bounded system"
function ldos(params, ind, η, angle = 0.0; order = 300000, bandrange = missing, kw...)
    h = rectangle_randombounds_sc(params, angle, η; kw...);
    si = SMatrix{4,4}(Diagonal(SA[1,1,-1,-1]))
    km = ketmodel(si, indices = ind)
    dos = dosKPM(h, ket = km, bandrange = bandrange; order = order)
    return dos
end

"""
ldosonlattice_rand_sc(p, η, angle = 0; nummodes = 4, selfy = false, kw...)
returns the wavefunctions of nev inner eigenvalues for a system, created by 
       rectangle_randombounds_sc()
"""
function ldosonlattice_rand_sc(p, η, angle = 0; nummodes = 4, selfy = false, kw...)
    ho = rectangle_randombounds_sc(p, angle, 0, selfy = selfy);
    hd = rectangle_randombounds_sc(p, angle, η, selfy = selfy; kw...);
    posd = sitepositions(hd)
    positionsdisordered(r) = r in posd ? true : false
    disord_indices = siteindices(ho, region = positionsdisordered)
    s = spectrum(hd, method = ArpackPackage(sigma = 1e-7, nev = ifelse(nummodes == 4, 4, 8)))
    println(s.energies)
    if nummodes == 4
        psi = Vector{Float64}(undef, Int(size(hd,1)*4))
        psi .=  abs.(s.states[:,1])
        psi .+= abs.(s.states[:,2])
        psi .+= abs.(s.states[:,3])
        psi .+= abs.(s.states[:,4])
        psi_ord = zeros(Float64, Int(size(ho,1)*4))
        count = 1
        for i in disord_indices
            psi_ord[4*i-3:4*i] = psi[count:count+3]
            count += 4
        end
        return psi_ord, psi_ord
    else nummodes == 8
        psi1 =  Vector{Float64}(undef, Int(size(hd,1)*4))
        psi2 =  Vector{Float64}(undef, Int(size(hd,1)*4))
        psi1 .=  abs.(s.states[:,1])
        psi1 .+= abs.(s.states[:,2])
        psi1 .+= abs.(s.states[:,7])
        psi1 .+= abs.(s.states[:,8])
        psi2 .=  abs.(s.states[:,3])
        psi2 .+= abs.(s.states[:,4])
        psi2 .+= abs.(s.states[:,5])
        psi2 .+= abs.(s.states[:,6])
        psi_ord1 = zeros(Float64, Int(size(ho,1)*4))
        psi_ord2 = zeros(Float64, Int(size(ho,1)*4))
        count = 1
        for i in disord_indices
            psi_ord1[4*i-3:4*i] = psi1[count:count+3]
            psi_ord2[4*i-3:4*i] = psi2[count:count+3]
            count += 4
        end
        return psi_ord1, psi_ord2
    end
end