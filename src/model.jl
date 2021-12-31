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
    Ws = 0
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

function latBLG_unbounded(p = Params())
    @unpack a0, dinter = p
    sAbot = sublat((0.0,-1.0a0/sqrt(3.0), - dinter/2); name = :Ab)
    sBbot = sublat((0.0, 0.0a0/sqrt(3.0), - dinter/2); name = :Bb)
    sAtop = sublat((0.0, 0.0a0/sqrt(3.0), + dinter/2); name = :At)
    sBtop = sublat((0.0, 1.0a0/sqrt(3.0), + dinter/2); name = :Bt)
    br = a0 * SA[cos(pi/3) sin(pi/3) 0; -cos(pi/3) sin(pi/3) 0]'
    lat = lattice(sAtop, sBtop, sAbot, sBbot; bravais = br)
    return lat
end


function latBLG(p, θ)
    (; a0, dinter, Ln, W) = p
    lat0_slg = LP.honeycomb(; a0, dim = 3)
    function rotated_region(r)
        r´ = SA[r[1], floor(abs(r[2])/(0.5*√3 * a0)) * sign(r[2])*(0.5*√3 * a0), r[3]]
        r´ = SA[cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1] * r´
        return abs(r´[1]) <= Ln/2 && abs(r´[2]) <= W/2
    end
    lat_slg = unitcell(lat0_slg, region = rotated_region)
    Quantica.transform!(r -> SA[cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1] * r, lat_slg)
    lat_bot = lattice(Quantica.transform!(r -> r + SA[cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1] *
         SA[0, -0.5*a0/sqrt(3.0), -dinter/2], copy(lat_slg)); names = (:Ab, :Bb))
    lat_top = lattice(Quantica.transform!(r -> r + SA[cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1] * 
        SA[0, 0.5*a0/sqrt(3.0), dinter/2], copy(lat_slg)); names = (:At, :Bt))
    return lat_top, lat_bot
end

latSLG(p = Params()) = LP.honeycomb(; a0 = p.a0, dim = Val(3))

function modelN(p = Params())
    @unpack a0, dinter, dintra, dinter_warp, Ls, Ln, W, μN, μS, t0, t1,
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

    modelinter = hopping(t1 * σ0τz, range = dinter, sublats = (:At=>:Bb, :Bb=>:At)) +
                 hopping(t3 * σ0τz, range = dinter_warp, sublats = (:Bt=>:Ab, :Ab=>:Bt))

    hopI(λ, dr) = λ/(2*3*sqrt(3)) * im * cos(3*atan(dr[2],dr[1])) * σzτ0

    modelIsing = hopping((r, dr) -> hopI(λ,dr) * sign(r[3]),
        sublats = (:At, :Bt ,:Ab, :Bb) .=> (:At, :Bt ,:Ab, :Bb), range = a0)

    modelKaneMele =
        hopping((r, dr) -> 0*hopI(λ, dr), sublats = :A=> :A, range = a0) -
        hopping((r, dr) -> 0*hopI(λ, dr), sublats = :B=> :B, range = a0)

    hopR(α, dr) = 2α/3 * im * (dr[2]/dintra * σxτ0 - dr[1]/dintra * σyτz)

    modelRashba =
        hopping((r,dr) -> hopR(α, dr), range = a0/√3, 
            sublats = (:At=>:Bt, :Bt=>:At, :A=>:B, :B=>:A)) -
        hopping((r,dr) -> hopR(α, dr), range = a0/√3, sublats = (:Ab=>:Bb, :Bb=>:Ab))

    modelonsiteN = onsite(-μN * σ0τz + EZ' * SA[σxτz, σyτ0, σzτz])

    model0 = modelonsiteN + modelintra + modelinter + modelIsing + 
        modelKaneMele + modelRashba

    field! = @onsite!((o, r) -> o + E'r * σ0τz)

    return (; model0, field!, modelinter)
end