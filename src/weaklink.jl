"""
Implements configuration B, i.e. BLG in Bernal stacking with self energy on both the ZZ and
the AC edges. There is a different phase difference between left and right sides of the 
weak link (-ϕ/2 for x<0 ϕ/2 for x>0).
Neither smoothness nor leads are implemented. Also disorder and rotation are not considered.
"""
function rectangle_weaklink(p)
    (; Ls, Δ, Ws) = p
    (; model0, field!, modelinter) = modelS(p)
    lat_top, lat_bot = MBSinBLG.latBLG(p, 0) #0 angle rotation
    
    # LOCAL SELF-ENERGY MODEL
    xmin, xmax = extrema(r->r[1], sitepositions(Quantica.combine(lat_top, lat_bot)))
    ymin, ymax = extrema(r->r[2], sitepositions(Quantica.combine(lat_top, lat_bot)))
    self_region(r) = !(xmin+Ls <= r[1] <= xmax - Ls) || !(ymin+Ws <= r[2] <= ymax - Ws)
    diagphi(phi) = Diagonal(SA[cis(phi), cis(phi), cis(-phi), cis(-phi)])
    sCself! = @onsite!((o, r; phi) -> o + 
        self_region(r) * Δ * diagphi(sign(r[1])* phi/4) * σyτy * diagphi(sign(r[1])*phi/4)') 
                                                              
    # HAMILTONIAN BUILD
    h_top = lat_top |> hamiltonian(model0; orbitals = Val(4)) |> 
        unitcell(mincoordination = 5)    
    h_bot = lat_bot |> hamiltonian(model0; orbitals = Val(4)) |>
        unitcell(mincoordination = 5)
    ph = Quantica.combine(h_top, h_bot; coupling = modelinter) |> 
        parametric(sCself!) 
            return ph
end

"""
    spectrumvsphase(philist, p; nev = 16, kw...)
computes a number nev of eigenvalues as a function of the phase difference across the weak
link. 
    see: rectangle_weakling(::Params()), Params(), modelS()
"""
function spectrumvsphase(philist, p; nev = 16, kw...)
    ph = rectangle_weaklink(p; kw...)
    elist = zeros(Float64, nev, length(philist))
    for i in 1:length(philist)
        s = spectrum(ph(phi = philist[i]), method = ArpackPackage(sigma = 1e-7, nev = nev))
        elist[:,i] = s.energies
    end
    return philist, elist
end
