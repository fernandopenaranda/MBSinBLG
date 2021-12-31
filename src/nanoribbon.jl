function nanoribbonN(p, (nx, ny))
    (; Ln, Ls) = p
    (; model0, field!) = modelN(p)
    lat_top, lat_bot = latBLG(p)
    perpaxis = Quantica.normalize(bravais(lat) * SA[1 1; -1 1] * SA[-ny, nx])

    h = lat |> hamiltonian(model0; orbitals = Val(2)) |>
               unitcell((1, -1), (1, 1)) |>
               unitcell((nx, ny), region = r -> abs(dot(perpaxis, r)) < Ln/2 + Ls, 
                modifiers = (field!,))
    return h
end

function nanoribbonS(p, θ)
    (; Ln, Ls, Δ, a0, τ, d, Ws) = p
    (; model0, field!, modelinter) = modelS(p)
    lat_top, lat_bot = latBLG(p, θ)
   
    # regionS(r) = abs(dot(perpaxis, r)) > ifelse(iszero(Ls), (round(0.5*Ln/a0)-WS)*a0, Ln/2) #0.9
    # regionS(r, dr) = regionS(r)
    regionS(r) = abs(r[1]) >= Ln/2 - Ls || abs(r[2]) >= W/2 - Ws
    regionS(r, dr) = regionS(r)

    smooth(r) = iszero(d) ? ifelse(abs(dot(perpaxis, r)) < Ln/2 - dxS * a0, 0.0, 1.0) :
                            1-(1+tanh((Ln/2-abs(dot(perpaxis, r)))/d))/2

    SCo! = @onsite!((o, r) -> o + smooth(r) * Δ * σyτy)
    SCh! = @hopping!((t, r, dr) -> t * ifelse(iszero(d),  1.0, (1.0 - smooth(r)));
        sublats = (:At, :Bt ,:Ab, :Bb, :A, :B) .=> (:At, :Bt ,:Ab, :Bb, :A, :B))
    SCτ! = @hopping!(t -> t * τ; region = (r, dr) -> xor(regionS(r - dr/2), 
        regionS(r + dr/2)))

    h_top = lat_top |> hamiltonian(model0; orbitals = Val(4)) 
    h_bot = lat_bot |> hamiltonian(model0; orbitals = Val(4)) 
    h = Quantica.combine(h_top, h_bot, coupling = modelinter)# |> unitcell((1, -1), (1, 1), modifiers = (field!, SCo!, SCh!, SCτ!))

    # h = lat |> hamiltonian(model0; orbitals = Val(4)) |>
    #            unitcell((1, -1), (1, 1)) |>
    #            unitcell((nx, ny), region = r -> abs(dot(perpaxis, r)) < Ln/2 + Ls, 
    #            modifiers = (field!, SCo!, SCh!, SCτ!))
    return h
end

nanoribbonNA(p = Params(); kw...) = nanoribbonN(p, 0; kw...)
nanoribbonNZ(p = Params(); kw...) = nanoribbonN(p, π/2; kw...)
nanoribbonSA(p = Params(); kw...) = nanoribbonS(p, 0; kw...)
nanoribbonSZ(p = Params(); kw...) = nanoribbonS(p, π/2; kw...)