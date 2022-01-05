"""
SQUID geometry. Self energy on the zigzag edge with a phase difference (-ϕ/2 for x<0 ϕ/2 for x>0)
Neither smoothness nor leads are implemented. Also disorder and rotation are not considered.
"""
function rectangle_squids(p)
    (; Ls, Δ, Ws) = p
    (; model0, field!, modelinter) = modelS(p)
    lat_top, lat_bot = latBLG(p, 0) #0 angle rotation
    
    # LOCAL SELF-ENERGY MODEL
    xmin, xmax = extrema(r->r[1], sitepositions(Quantica.combine(lat_top, lat_bot)))
    ymin, ymax = extrema(r->r[2], sitepositions(Quantica.combine(lat_top, lat_bot)))
    self_region(r) = !(xmin+Ls <= r[1] <= xmax - Ls) || !(ymin+Ws <= r[2] <= ymax - Ws)
    diagphi(ϕ) = Diagonal(SA[cis(ϕ), cis(ϕ), cis(-ϕ), cis(-ϕ)])
    sCself! = @onsite!((o, r; ϕ) -> o + 
        self_region(r) * Δ * diagphi(sign(r[1])* ϕ/4) * σyτy * diagphi(sign(r[1])*ϕ/4)') 
                                                              
    # HAMILTONIAN BUILD
    h_top = lat_top |> hamiltonian(model0; orbitals = Val(4)) |> 
        unitcell(mincoordination = 5, region = !isvacancy)    
    h_bot = lat_bot |> hamiltonian(model0; orbitals = Val(4)) |>
        unitcell(mincoordination = 5, region = !isvacancy)
    ph = Quantica.combine(h_top, h_bot; coupling = modelinter) |> 
        parametric(sCself!) 
            return ph
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

