## Phase diagram (superconducting, armchair, optim)

function runfig3()
    p = Params(Ln = 1440, Ls = 0, scale = 40, 
        λ = 5, α = 0, EZ = SA[0, 1e-3, 0], Δ = 0.3, d = 0, τ = 1);
    ha = nanoribbonSA(p);
    hz = nanoribbonSZ(p);
    pha = parametric(ha, @onsite!((o; dμ = 0, dEZy = 0) -> o - dμ * σ0τz +  dEZy*σyτ0))
    phz = parametric(hz, @onsite!((o; dμ = 0, dEZy = 0) -> o - dμ * σ0τz +  dEZy*σyτ0))
    mata = similarmatrix(pha, Quantica.flatten)
    matz = similarmatrix(phz, Quantica.flatten)
    dμs = range(0., 2.5, 10)
    boundarySA = @showprogress [optimize(fA(dμ, mata, pha), 0.00001, 2).minimizer for dμ in dμs]
    boundarySZ = @showprogress [optimize(fZ(dμ, matb, phb), 0.00001, 2).minimizer for dμ in dμs]
    return dμs, boundarySA, boundarySZ
end

function gapA( mat, ph; kw...)
    bloch!(mat, ph(; kw...), 0)
    λ = maximum(real.(first(eigs(mat; sigma = 0.00001, nev = 4))))
    return λ
end

function gapZ(mat, ph; kw...)
    k0 = 2.5692321  # that's the Bloch phase where normal zigzag bands cross zero
    bloch!(mat, ph(; kw...), k0)
    λ = maximum(real.(first(eigs(mat; sigma = 0.00001, nev = 4))))
    return λ
end

fA(dμ, mat, ph) = x -> gapA(mat, ph; dEZy = x, dμ)
fZ(dμ, mat, ph) = x -> gapZ(mat, ph; dEZy = x, dμ)