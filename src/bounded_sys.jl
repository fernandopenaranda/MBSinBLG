
"builds the hamiltonian of BLG with a random model implemented 
on the edges of the system see `rectangle` we also add two possible SC configurations
(square mask not yet implemented),
the degree of disorder is parametrized by a value `η` which determines the average number
of atom vacancies in the edges. 
    `η = 0` -> crystalographic edges
    `η = 0.` is the default 
we can control the rotation of the lattices w.r.t the SC using the keyword `θ`"

function rectangle_randombounds_sc(p, θ = 0, η = 0.; sidecontacts = false, 
        selfy = false, ϕ0 = 0)
    (; Ln, Ls, Δ, a0, τ, d, W, Ws) = p
    (; model0, field!, modelinter) = modelS(p)
    lat_top, lat_bot = latBLG(p, θ)
       
    # contact model (Josephson vs rectangular mask architectures)
    smooth_sides(r) = iszero(d) ? ifelse(abs(r[1])<Ln/2 && abs(r[2])<W/2, 0.0, 1.0) :
            1 - (1 + tanh((Ln/2-abs(r[1]))/d))/2
    smooth_rectangle(r) = iszero(d) ? ifelse(abs(r[1])<Ln/2 && abs(r[2])<W/2, 0.0, 1.0) :
                1 - (1 + tanh((Ln/2-abs(r[1]))/d))*(1 + tanh((W/2-abs(r[2]))/d))/4
    smooth_method(r) = ifelse(sidecontacts == true , smooth_sides(r), smooth_rectangle(r))
    
    # RANDOM MODEL (we randomly replace sites within rand_region by vacancies) 
    rand_region(r) = Ln/2 - 0.7a0 <= abs(r[1]) <= Ln/2 #&& Ln/2 - 2a0 <= abs(r[2]) <= Ln/2
    random_mat() = sample([0, 1], StatsBase.Weights([1-η, η])) * σ0τz
    so_rand! = @onsite!((o, r) -> o + rand_region(r) * ifelse(η != 0, 1e3, 0) * random_mat())
    
    # LOCAL SELF-ENERGY MODEL

    xmin, xmax = extrema(r->r[1], sitepositions(Quantica.combine(lat_top, lat_bot)))
    ymin, ymax = extrema(r->r[2], sitepositions(Quantica.combine(lat_top, lat_bot)))

    self_region(r) = ifelse(selfy == false, abs(r[1]) >= Ln/2 - Ls,  
        !(xmin+Ls <= r[1] <= xmax - Ls) || !(ymin+Ws <= r[2] <= ymax - Ws))
    diagphi(ϕ) = Diagonal(SA[cis(ϕ), cis(ϕ), cis(-ϕ), cis(-ϕ)])
    sCself! = @onsite!((o, r; ϕ) -> o + 
            self_region(r) * Δ * diagphi(sign(r[1])* ϕ/4) * σyτy * diagphi(sign(r[1])*ϕ/4)')     
    
    # ONSITE SC MODEL (transparencies... disabled by default)
    # SCo! = @onsite!((o, r) -> o + smooth_method(r) * Δ * σyτy)
    # SCh! = @hopping!((t, r, dr) -> t * (1-smooth_method(r));
    #     sublats = (:At, :Bt ,:Ab, :Bb, :A, :B) .=> (:At, :Bt ,:Ab, :Bb, :A, :B))
    # SCτ! = @hopping!((t, r, dr) -> t * τ; region = (r, dr) -> xor(regionS(r - dr/2),
    #     regionS(r + dr/2)))
    
    # HAMILTONIAN BUILD
    h_top = lat_top |> hamiltonian(model0; orbitals = Val(4)) |> unitcell(mincoordination = 5)    
    h_bot = lat_bot |> hamiltonian(model0; orbitals = Val(4)) |> unitcell(mincoordination = 5)
     ph = Quantica.combine(h_top, h_bot; coupling = modelinter) |> 
        parametric(field!, so_rand!, sCself!) #, SCh!, SCτ!)
    return ph(ϕ = ϕ0)
end
