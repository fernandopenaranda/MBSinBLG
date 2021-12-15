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

