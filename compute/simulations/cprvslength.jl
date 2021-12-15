using MBSinBLG, Quantica, Arpack, Parameters
using ElectronDisplay

function test(; α = 0, Ey = 4)
    p = Params(Ls = 0, scale = 40, λ = 5, α = α,   
        μN = 0.6, Δ = 1, d = 0, τ = 1, EZ = SA[0, Ey, 0]);
    lnlist = 50:25:500
    x, y = cprvsln(lnlist, p, nev = 8)
    return x, y
end

function cprvsln(list, p; nev = 16, kw...)
    elist = zeros(Float64, nev, length(list))
    for i in 1:length(list)
        ph = MBSinBLG.rectangle_squid(reconstruct(p, Ln = list[i]); kw...)
        s = spectrum(ph(ϕ = π), method = ArpackPackage(sigma = 1e-7, nev = nev))
        elist[:,i] = s.energies
    end
    return list, elist
end


x_b, y_b = test(α = 0, Ey = 4)
x_d, y_d = test(α = 1, Ey = 4)