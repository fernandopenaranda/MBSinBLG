
"builds the hamiltonian of BLG with a random model implemented 
on the edges of the system see `rectangle` we also add two possible SC configurations
(square mask not yet implemented),
the degree of disorder is parametrized by a value `η` which determines the average number
of atom vacancies in the edges. 
    `η = 0` -> crystalographic edges
    `η = 0.` is the default 
we can control the rotation of the lattices w.r.t the SC using the keyword `θ`"

function rectangle_randombounds_sc(p, θ = 0, η = 0.; sidecontacts = false, 
        selfy = false, mono = false, ϕ0 = 0)
    (; Ln, Ls, Δ, a0, τ, d, W) = p
    (; model0, field!) = modelS(p)
    lat0 = mono ? latSLG(p) : latBLG(p)
    rot(r) = SA[cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1] * r
    Quantica.transform!(r -> rot(r), lat0)
    # contact model (Josephson vs rectangular mask architectures)
    if sidecontacts == true
        myreg = RP.rectangle( (Ln + 2Ls, W))
    else
        myreg = RP.rectangle( (Ln + 2Ls, Ln + 2Ls))
    end
    # lat = unitcell(lat0, region = myreg)
    lat = unitcell(lat0, 
        region = (r) -> r[1]> -Ln/2 && r[1]< Ln/2 && r[2] < W/2 && r[2]> -W/2)
    
    regionS(r) = abs(r[1]) >= Ln/2 || abs(r[2]) >= W/2
    regionS(r, dr) = regionS(r)
    
    smooth_sides(r) = iszero(d) ? ifelse(abs(r[1])<Ln/2 && abs(r[2])<W/2, 0.0, 1.0) :
            1 - (1 + tanh((Ln/2-abs(r[1]))/d))/2

    smooth(r) = iszero(d) ? ifelse(abs(r[1])<Ln/2 && abs(r[2])<W/2, 0.0, 1.0) :
                1 - (1 + tanh((Ln/2-abs(r[1]))/d))*(1 + tanh((W/2-abs(r[2]))/d))/4

    smooth_method(r) = ifelse(sidecontacts == true , smooth_sides(r), smooth(r))
    
    #random model (we are randomly replacing sites by vacancies) within rand_region
    rand_region(r) = Ln/2 - 3a0 <= abs(r[1]) <= Ln/2 #&& Ln/2 - 2a0 <= abs(r[2]) <= Ln/2
    random_mat() = 
        sample([0, 1], StatsBase.Weights([1-η, η])) * σ0τz
    so_rand! = @onsite!((o, r) -> o + rand_region(r) * ifelse(η != 0, 1e3, 0) *
        random_mat())
    # local self energy model
    self_region(r) = ifelse(selfy == false, abs(r[1]) ≥ (round(0.5*Ln/a0)-0.9)*a0,  
        # r[1]> -Ln/2 && r[1]< Ln/2 && r[2] < W/2 && r[2]> -W/2)
        abs(r[1]) ≥ (round(0.5*Ln/a0)-2)*a0 || abs(r[2]) ≥ (round(0.5*W/a0)-2)*a0)

    # self_region(r) = ifelse(selfy == false, abs(r[1]) ≥ (round(0.5*Ln/a0)-0.9)*a0,  
    # r[1] ≤ -tan(θ)*r[2]-(round(0.5*Ln/a0)-2)*a0 ||  
    # r[1] ≥ -tan(θ)*r[2]+(round(0.5*Ln/a0)-2)*a0 || r[2] ≥ tan(θ)*(r[1]-0.5*Ln/a0) +
    #     (round(0.5*W/a0)-2)*a0 || r[2] ≤ tan(θ)*(r[1]-0.5*Ln/a0) -(round(0.5*W/a0)-2)*a0 ) 


    diagphi(ϕ) = Diagonal(SA[cis(ϕ), cis(ϕ), cis(-ϕ), cis(-ϕ)])
    sCself! = @onsite!((o, r; ϕ) -> o + 
            self_region(r) * Δ * diagphi(sign(r[1])* ϕ/4) * σyτy * diagphi(sign(r[1])*ϕ/4)')     
    # Superconductivity
    SCo! = @onsite!((o, r) -> o + smooth_method(r) * Δ * σyτy)
    #sCself! = @onsite!((o, r) -> o + self_region(r) * Δ * σyτy)
    SCh! = @hopping!((t, r, dr) -> t * (1-smooth_method(r));
        sublats = (:At, :Bt ,:Ab, :Bb, :A, :B) .=> (:At, :Bt ,:Ab, :Bb, :A, :B))
    SCτ! = @hopping!((t, r, dr) -> t * τ; region = (r, dr) -> xor(regionS(r - dr/2),
        regionS(r + dr/2)))
    #building the hamiltonian
    ph = lat |> hamiltonian(model0; orbitals = Val(4)) |> 
        parametric(field!, so_rand!, SCo!, sCself!, SCh!, SCτ!)
    return ph(ϕ = ϕ0)# ph()
end
