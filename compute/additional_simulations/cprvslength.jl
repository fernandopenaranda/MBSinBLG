using MBSinBLG, Quantica, Arpack, Parameters
using ElectronDisplay

function test(; α = 0, Ey = 4, kw...)
    p = Params(Ls = 0, scale = 40, λ = 5, α = α,   
        μN = 0.6, Δ = 1, d = 0, τ = 1, EZ = SA[0, Ey, 0]);
    lnlist = 50:25:500
    x, y = cprvsln(lnlist, p, nev = 8; kw...)
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


x_b, y_b = test(α = 0, Ey = 4, W = 1440)
x_d, y_d = test(α = 1, Ey = 4, W = 1440)

x_b2, y_b2 = test(α = 0, Ey = 4, W = 644)
x_d2, y_d2 = test(α = 1, Ey = 4, W = 644)

function plottest(x, y,  ylims = (-0.1, 0.1))
    fig = Figure(resolution = (500, 300), font = "Times New Roman") 
    axa = Axis(fig[1,1])
   
    lines!(axa, x, y[1,:])
    for i in 2:size(y, 1)
        lines!(axa, x, y[i,:])
    end
    axa.xlabel = " Ln [nm]"
    axa.ylabel = "E [meV]"
    ylims!(axa, ylims)
    fig
end
