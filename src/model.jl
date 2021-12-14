using Quantica, StaticArrays, Parameters#, VegaLite
using PhysicalConstants.CODATA2018: ustrip, @u_str, ħ, k_B, m_e, e, μ_B

const ħoec = ustrip(u"T*nm^2",ħ/e)
const μB = ustrip(u"T^-1*meV",μ_B)
const σ0τz = @SMatrix[1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 -1]
const σ0τ0 = @SMatrix[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
const σzτ0 = @SMatrix[1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 -1]
const σzτz = @SMatrix[1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1]
const σyτy = @SMatrix[0 0 0 -1; 0 0 1 0; 0 1 0 0; -1 0 0 0]
const σyτz = @SMatrix[0 -im 0 0; im 0 0 0; 0 0 0 im; 0 0 -im 0]
const σyτ0 = @SMatrix[0 -im 0 0; im 0 0 0; 0 0 0 -im; 0 0 im 0]
const σxτz = @SMatrix[0 1 0 0; 1 0 0 0; 0 0 0 -1; 0 0 -1 0]
const σxτ0 = @SMatrix[0 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 1 0]
const σ0 = SA[1 0; 0 1]
const σx = SA[0 1; 1 0]
const σy = SA[0 -im; im 0]
const σz = SA[1 0; 0 -1]

@with_kw struct Params @deftype Float64
    # Units: nm, meV
    # Params AB bernal bilayer obtained from DFT calculations: 
    # ref -> https://arxiv.org/pdf/1511.06706.pdf (flag <-:dft)
    scale = 1                  #scaling parameter
    a0 = 0.246 * scale
    dinter = 1.36 * a0/scale  #Note that we leave dinter unaltered by scale
    dintra = a0/sqrt(3)
    dinter_warp = sqrt(dintra^2 + dinter^2)
    Ls = 0                      #units of a0
    Ln = 30                     #units of a0
    W = 30                      #units of a0
    d = 0                       # transition length of NS contact

    t0 =  1e3 * 2.6 / scale     #in meV
    t1 = 1e3 * .36              #in meV
    t3 = 1e3 * .28 / scale      #in meV
    μN = 0 * t0
    μS = t0 / 3.25              # t0/3 > μS > 2 * t1
    τ = 1                       # NS transparency
    barrier = 0                 # NS barrier

    U_dimer = 0 * 0.015         #Energy difference between dimer and non dimer sites: 
                                # https://arxiv.org/pdf/1901.01332.pdf
    U = 0                       #Interlayer bias
    Δ = 1.3                     #That of SC MoRe: ref ->  https://arxiv.org/abs/2006.05522
    EZ::SVector{3,Float64} = [0,20,0]
    E::SVector{3,Float64} = [0,0,0]
    λ = 30                      #Between [-2.5,2.5] meV: 
                                #ref -> https://arxiv.org/abs/1905.08688
    α = 20                      #Rashba coef
    g = 2
end

function latBLG(p = Params())
    @unpack a0, dinter = p
    sAbot = sublat((0.0,-1.0a0/sqrt(3.0), - dinter/2); name = :Ab)
    sBbot = sublat((0.0, 0.0a0/sqrt(3.0), - dinter/2); name = :Bb)
    sAtop = sublat((0.0, 0.0a0/sqrt(3.0), + dinter/2); name = :At)
    sBtop = sublat((0.0, 1.0a0/sqrt(3.0), + dinter/2); name = :Bt)
    br = a0 * SA[cos(pi/3) sin(pi/3) 0; -cos(pi/3) sin(pi/3) 0]'
    lat = lattice(sAtop, sBtop, sAbot, sBbot; bravais = br)
    return lat
end

latSLG(p = Params()) = LP.honeycomb(; a0 = p.a0, dim = Val(3))

function modelN(p = Params())
    @unpack a0, dinter, dintra, dinter_warp, Ls, Ln, W, μN,μS , t0, t1,
            t3, Δ, λ, g, U, U_dimer, α, scale, EZ, E = p

    modelintra = hopping(-t0 * σ0, range = dintra,
            sublats = (:At=>:Bt,:Ab=>:Bb,:Bb=>:Ab,:Bt=>:At, :A=>:B, :B=>:A))

    modelinter = hopping(t1 * σ0, range = dinter, sublats = (:At=>:Bb,:Bb=>:At)) +
                 hopping(t3 * σ0, range = dinter_warp, sublats = (:Bt=>:Ab, :Ab=>:Bt))

    hopI(λ, dr) = λ/(2*3*sqrt(3)) * im * cos(3*atan(dr[2],dr[1])) * σz

    modelIsing = hopping((r, dr) -> hopI(λ,dr) * sign(r[3]),
        sublats = (:At, :Bt ,:Ab, :Bb) .=> (:At, :Bt ,:Ab, :Bb), range = a0)

    modelKaneMele =
        hopping((r, dr) -> hopI(λ, dr), sublats = :A=>:A, range = a0) -
        hopping((r, dr) -> hopI(λ, dr), sublats = :B=>:B, range = a0)

    hopR(α, dr) = 2α/3 * im * (dr[2]/dintra * σx - dr[1]/dintra * σy)

    modelRashba =
        hopping((r,dr) -> hopR(α, dr), range = a0/√3, 
            sublats = (:At=>:Bt, :Bt=>:At, :A=>:B, :B=>:A)) -
        hopping((r,dr) -> hopR(α, dr), range = a0/√3, sublats = (:Ab=>:Bb, :Bb=>:Ab))

    modelonsiteN = onsite(-μN * σ0 + EZ'*SA[σx, σy, σz])

    model0 = modelonsiteN + modelintra + modelinter + modelIsing + 
        modelKaneMele + modelRashba

    field! = @onsite!((o, r) -> o + E'r * σ0)

    return (; model0, field!)
end

function modelS(p = Params())
    @unpack a0, dinter, dintra, dinter_warp, Ls, Ln, W, μN,μS , t0, t1,
            t3, Δ, λ, g, U, U_dimer, α, scale, EZ, E = p

    modelintra = hopping(-t0 * σ0τz, range = dintra,
            sublats = (:At=>:Bt,:Ab=>:Bb,:Bb=>:Ab,:Bt=>:At, :A=>:B, :B=>:A))

    modelinter = hopping(t1 * σ0τz, range = dinter, sublats = (:At=>:Bb,:Bb=>:At)) +
                 hopping(t3 * σ0τz, range = dinter_warp, sublats = (:Bt=>:Ab, :Ab=>:Bt))

    hopI(λ, dr) = λ/(2*3*sqrt(3)) * im * cos(3*atan(dr[2],dr[1])) * σzτ0

    modelIsing = hopping((r, dr) -> hopI(λ,dr) * sign(r[3]),
        sublats = (:At, :Bt ,:Ab, :Bb) .=> (:At, :Bt ,:Ab, :Bb), range = a0)

    modelKaneMele =
        hopping((r, dr) -> hopI(λ, dr), sublats = :A=>:A, range = a0) -
        hopping((r, dr) -> hopI(λ, dr), sublats = :B=>:B, range = a0)

    hopR(α, dr) = 2α/3 * im * (dr[2]/dintra * σxτ0 - dr[1]/dintra * σyτz)

    modelRashba =
        hopping((r,dr) -> hopR(α, dr), range = a0/√3, 
            sublats = (:At=>:Bt, :Bt=>:At, :A=>:B, :B=>:A)) -
        hopping((r,dr) -> hopR(α, dr), range = a0/√3, sublats = (:Ab=>:Bb, :Bb=>:Ab))

    modelonsiteN = onsite(-μN * σ0τz + EZ'*SA[σxτz, σyτ0, σzτz])

    model0 = modelonsiteN + modelintra + modelinter + modelIsing + 
        modelKaneMele + modelRashba

    field! = @onsite!((o, r) -> o + E'r * σ0τz)

    return (; model0, field!)
end

function nanoribbonN(p, (nx, ny); mono = false)
    (; Ln, Ls, U, a0) = p
    (; model0, field!) = modelN(p)
    lat = mono ? latSLG(p) : latBLG(p)
    perpaxis = Quantica.normalize(bravais(lat) * SA[1 1; -1 1] * SA[-ny, nx])

    regionS(r) = abs(dot(perpaxis, r)) > ifelse(iszero(Ls), (round(0.5*Ln/a0)-0.9)*a0, Ln/2)
    smooth(r, β = 5) = iszero(Ls) ? 1.0 : tanh(β*(abs(dot(perpaxis, r)) - Ln/2) / Ls)

    # U! = @onsite!((o, r) -> o + U * sign(r[3]) * smooth(r) * σ0; region = regionS)

    h = lat |> hamiltonian(model0; orbitals = Val(2)) |>
               unitcell((1, -1), (1, 1)) |>
               unitcell((nx, ny), region = r -> abs(dot(perpaxis, r)) < Ln/2 + Ls, 
                modifiers = (field!,))
    return h
end

function nanoribbonS(p, (nx, ny); mono = false, dxS = 1.)
    (; Ln, Ls, Δ, a0, τ, d) = p
    (; model0, field!) = modelS(p)
    lat = mono ? latSLG(p) : latBLG(p)
    perpaxis = Quantica.normalize(bravais(lat) * SA[1 1; -1 1] * SA[-ny, nx])

    regionS(r) = abs(dot(perpaxis, r)) > ifelse(iszero(Ls), (round(0.5*Ln/a0)-0.9)*a0, Ln/2)
    regionS(r, dr) = regionS(r)

    smooth(r) = iszero(d) ? ifelse(abs(dot(perpaxis, r)) < Ln/2 - dxS * a0, 0.0, 1.0) :
                            1-(1+tanh((Ln/2-abs(dot(perpaxis, r)))/d))/2

    SCo! = @onsite!((o, r) -> o + smooth(r) * Δ * σyτy)
    SCh! = @hopping!((t, r, dr) -> t * ifelse(iszero(d),  1.0, (1.0 - smooth(r)));
        sublats = (:At, :Bt ,:Ab, :Bb, :A, :B) .=> (:At, :Bt ,:Ab, :Bb, :A, :B))
    SCτ! = @hopping!(t -> t * τ; region = (r, dr) -> xor(regionS(r - dr/2), 
        regionS(r + dr/2)))

    h = lat |> hamiltonian(model0; orbitals = Val(4)) |>
               unitcell((1, -1), (1, 1)) |>
               unitcell((nx, ny), region = r -> abs(dot(perpaxis, r)) < Ln/2 + Ls, 
               modifiers = (field!, SCo!, SCh!, SCτ!))
    return h
end

nanoribbonNA(p = Params(); kw...) = nanoribbonN(p, (0, 1); kw...)
nanoribbonNZ(p = Params(); kw...) = nanoribbonN(p, (1, 0); kw...)
nanoribbonSA(p = Params(); kw...) = nanoribbonS(p, (0, 1); kw...)
nanoribbonSZ(p = Params(); kw...) = nanoribbonS(p, (1, 0); kw...)

function square(p, θ = 0; mono = false)
    (; Ln, Ls, Δ, a0, τ, d) = p
    (; model0, field!) = modelS(p)
    lat0 = mono ? latSLG(p) : latBLG(p)
    Quantica.transform!(r -> SA[cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1] * r, lat0)
    lat = unitcell(lat0, region = RP.square(Ln + 2Ls))

    regionS(r) = abs(r[1]) >= Ln/2 || abs(r[2]) >= Ln/2
    regionS(r, dr) = regionS(r)

    smooth(r) = iszero(d) ? ifelse(abs(r[1])<Ln/2 && abs(r[2])<Ln/2, 0.0, 1.0) :
        1 - (1 + tanh((Ln/2-abs(r[1]))/d))*(1 + tanh((Ln/2-abs(r[2]))/d))/4

    SCo! = @onsite!((o, r) -> o + smooth(r) * Δ * σyτy)
    SCh! = @hopping!((t, r, dr) -> t * (1-smooth(r));
        sublats = (:At, :Bt ,:Ab, :Bb, :A, :B) .=> (:At, :Bt ,:Ab, :Bb, :A, :B))
    SCτ! = @hopping!((t, r, dr) -> t * τ; region = (r, dr) -> xor(regionS(r - dr/2),
        regionS(r + dr/2)))

    ph = lat |> hamiltonian(model0; 
        orbitals = Val(4)) |> parametric(field!, SCo!, SCh!, SCτ!)
    h = ph()

    return h
end

function rectangle(p, θ = 0; mono = false, sidecontacts = false, dims = "homog")
    (; Ln, Ls, Δ, a0, τ, d) = p
    (; model0, field!) = modelS(p)
    lat0 = mono ? latSLG(p) : latBLG(p)
    Quantica.transform!(r -> SA[cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1] * r, lat0)

    if sidecontacts == true
        myreg = RP.rectangle( (Ln + 2Ls, Ln))
    else
        myreg = RP.rectangle( (Ln + 2Ls, Ln + 2Ls))
    end

    lat = unitcell(lat0, region = myreg)
    regionS(r) = abs(r[1]) >= Ln/2 || abs(r[2]) >= Ln/2
    regionS(r, dr) = regionS(r)
    
        smooth_sides(r) = iszero(d) ? ifelse(abs(r[1])<Ln/2 && abs(r[2])<Ln/2, 0.0, 1.0) :
                1 - (1 + tanh((Ln/2-abs(r[1]))/d))/2
    
        smooth(r) = iszero(d) ? ifelse(abs(r[1])<Ln/2 && abs(r[2])<Ln/2, 0.0, 1.0) :
                 1 - (1 + tanh((Ln/2-abs(r[1]))/d))*(1 + tanh((Ln/2-abs(r[2]))/d))/4
    
    smooth_method(r) = ifelse(sidecontacts == true , smooth_sides(r), smooth(r))
    SCo! = @onsite!((o, r) -> o + smooth_method(r) * Δ * σyτy)
    SCh! = @hopping!((t, r, dr) -> t * (1-smooth_method(r));
        sublats = (:At, :Bt ,:Ab, :Bb, :A, :B) .=> (:At, :Bt ,:Ab, :Bb, :A, :B))
    SCτ! = @hopping!((t, r, dr) -> t * τ; region = (r, dr) -> xor(regionS(r - dr/2),
        regionS(r + dr/2)))

    ph = lat |> hamiltonian(model0; 
        orbitals = Val(4)) |> parametric(field!, SCo!, SCh!, SCτ!)
    h = ph()
    return h
end