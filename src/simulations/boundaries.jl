## Phase diagram (superconducting, armchair, optim)

prs = Params(Ln = 1440, Ls = 0, scale = 40, 
        λ = 5, α = 0, EZ = SA[0, 1e-3, 0], Δ = 0.3, d = 0, τ = 1);


using Optim, ProgressMeter
h = nanoribbonSA(prs);
ph = parametric(h, @onsite!((o; dμ = 0, dEZy = 0) -> o - dμ * σ0τz +  dEZy*σyτ0))
mat = similarmatrix(ph, Quantica.flatten)

function gap(; mat = mat, ph = ph, kw...)
    bloch!(mat, ph(; kw...), 0)
    λ = maximum(real.(first(eigs(mat; sigma = 0.00001, nev = 4))))
    return λ
end

fu(dμ) = x -> gap(; dEZy = x, dμ)
##
dμs = range(0., 2.5, 100)
boundarySA = @showprogress [optimize(fu(dμ), 0.00001, 2).minimizer for dμ in dμs]

## Phase diagram (superconducting, zigzag, optim)
using Optim, ProgressMeter
h = nanoribbonS(prs, (1,0));
ph = parametric(h, @onsite!((o; dμ = 0, dEZy = 0) -> o - dμ * σ0τz +  dEZy*σyτ0))
mat = similarmatrix(ph, Quantica.flatten)

function gapZ(; mat = mat, ph = ph, kw...)
    bloch!(mat, ph(; kw...), 2.5692321)  # that's the Bloch phase where normal zigzag bands cross zero
    λ = maximum(real.(first(eigs(mat; sigma = 0.00001, nev = 4))))
    return λ
end

fl(dμ) = x -> gapZ(; dEZy = x, dμ)
##
dμs = range(0., 2.5, 100)
boundarySZ = @showprogress [optimize(fl(dμ), 0.00001, 2).minimizer for dμ in dμs]