"builds the hamiltonian of BLG with a random model implemented 
on the edges of the system see `rectangle` we also add two possible SC configurations: 
A and B in the manuscript.
The amount of disorder is parametrized by a value `η` which determines the average number
of atom vacancies in the edges. 
    `η = 0` -> crystalographic edges
we can control the rotation of the lattices w.r.t the SC leads using the keyword `θ`"

function rectangle_randombounds_sc(p, θ = 0, η = 0.; selfy = false, ϕ0 = 0)
    (; Ls, Δ, Ws) = p
    (; model0, field!, modelinter) = modelS(p)
    lat_top, lat_bot = latBLG(p, θ)
             
    # LOCAL SELF-ENERGY MODEL
    xmin, xmax = extrema(r->r[1], sitepositions(Quantica.combine(lat_top, lat_bot)))
    ymin, ymax = extrema(r->r[2], sitepositions(Quantica.combine(lat_top, lat_bot)))
    self_region(r) = ifelse(selfy == false, !(xmin+Ls <= r[1] <= xmax - Ls),  
        !(xmin+Ls <= r[1] <= xmax - Ls) || !(ymin+Ws <= r[2] <= ymax - Ws))
    diagphi(ϕ) = Diagonal(SA[cis(ϕ), cis(ϕ), cis(-ϕ), cis(-ϕ)])
    sCself! = @onsite!((o, r; ϕ) -> o + 
        self_region(r) * Δ * diagphi(sign(r[1])* ϕ/4) * σyτy * diagphi(sign(r[1])*ϕ/4)')     
    
    # RANDOM MODEL (we randomly replace sites within self_energy by vacancies) 
    vacancies(r) = sample([false, true], StatsBase.Weights([1-η, η])) * self_region(r)
    
    # HAMILTONIAN BUILD
    h_top = lat_top |> hamiltonian(model0; orbitals = Val(4)) |>
        unitcell(mincoordination = 5, region = !vacancies)    
    h_bot = lat_bot |> hamiltonian(model0; orbitals = Val(4)) |>
        unitcell(mincoordination = 5, region = !vacancies)
     ph = Quantica.combine(h_top, h_bot; coupling = modelinter) |> 
        parametric(field!, sCself!)
    ## checking vacancy probability
    if η !=0
        sitesinhdis = length(collect(siteindices(ph(ϕ=0))))
        println(sitesinhdis)
        normalsites, allsites = 
            sitesinhord(lat_top, lat_bot, self_region, model0, modelinter)
        println("vacancy probability: ", 1-(sitesinhdis-normalsites)/(allsites-normalsites))
    else nothing 
    end
    return ph(ϕ = ϕ0)
end

"""
    sitesinhord(latt, latb, self_region, model, modelinter)
returns the number of normal sites and the total number of sites in the non-disordered lattice 
    see: rectangle_randombounds_sc()
"""
function sitesinhord(latt, latb, self_region, model, modelinter)
    h_top = latt |> hamiltonian(model; orbitals = Val(4))
    h_bot = latb |> hamiltonian(model; orbitals = Val(4))
    normalsites = length(collect(siteindices(Quantica.combine(h_top |> 
        unitcell(mincoordination = 5, region = !self_region), 
            h_bot |> unitcell(mincoordination = 5, region = !self_region); 
                coupling = modelinter))))
    allsites = length(collect(siteindices(Quantica.combine(h_top |> 
        unitcell(mincoordination = 5), h_bot |> 
            unitcell(mincoordination = 5); coupling = modelinter))))
    return normalsites, allsites
end

# Deprecated terms in H (explicit construction of the leads)
# ONSITE SC MODEL (transparencies... disabled by default)
# SCo! = @onsite!((o, r) -> o + smooth_method(r) * Δ * σyτy)
# SCh! = @hopping!((t, r, dr) -> t * (1-smooth_method(r));
#     sublats = (:At, :Bt ,:Ab, :Bb, :A, :B) .=> (:At, :Bt ,:Ab, :Bb, :A, :B))
# SCτ! = @hopping!((t, r, dr) -> t * τ; region = (r, dr) -> xor(regionS(r - dr/2),
#     regionS(r + dr/2)))
    # contact model (Josephson vs rectangular mask architectures)
    # smooth_sides(r) = iszero(d) ? ifelse(abs(r[1])<Ln/2 && abs(r[2])<W/2, 0.0, 1.0) :
    #         1 - (1 + tanh((Ln/2-abs(r[1]))/d))/2
    # smooth_rectangle(r) = iszero(d) ? ifelse(abs(r[1])<Ln/2 && abs(r[2])<W/2, 0.0, 1.0) :
    #             1 - (1 + tanh((Ln/2-abs(r[1]))/d))*(1 + tanh((W/2-abs(r[2]))/d))/4
    # smooth_method(r) = ifelse(sidecontacts == true , smooth_sides(r), smooth_rectangle(r))