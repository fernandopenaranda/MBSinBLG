function runfig7()
    p = Params(Ln = 195, W = 2500, Ls = 20, Ws = 20, scale = 40, λ = 5,
    Δ = 1, d = 0, τ = 1);
    ϕlist = -2π:π/10:2π
    println("1")
    _, esa = spectrumvsphase(ϕlist, 
        reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6, α = 0), nev = 16)
    println("2")
    _, esb = spectrumvsphase(ϕlist, 
        reconstruct(p, EZ = SA[0, 4, 0], μN = 1.192, α = 0), nev = 16)
    println("3")
    _, esc = spectrumvsphase(ϕlist, 
        reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6, α = 2), nev = 16)
    println("4")
    ϕs, esd = spectrumvsphase(ϕlist, 
        reconstruct(p, EZ = SA[0, 4, 0], μN = 1.192, α = 2), nev = 16)
    savespectrumvsphase("fig7", p, ϕs, esa, esb, esc, esd)    
end

# Extra: checking W and Ln dependence

# function Wtest(Wlist, p, ϕ0, panel; nev = 4, kw...)
#     elist = zeros(Float64, nev, length(Wlist))
#     for i in 1:length(Wlist)
#         p = reconstruct(p, W = Wlist[i])
#         ph = rectangle_squid(p; kw...)
#         s = spectrum(ph(phi = ϕ0), method = ArpackPackage(sigma = 1e-7, nev = nev))
#         elist[:,i] = s.energies
    
#         println(s.energies, " W: ",Wlist[i])
#     end
#     return Wlist, elist
# end

# function testW(Ln0, ϕ0, panel)
#     Wlist = 2000:200:2500
#     p = Params(Ln = Ln0, Ls = 20, Ws = 20, scale = 40, λ = 5, Δ = 1, d = 0, τ = 1);
#     if panel == "a"
#         p = reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6, α = 0)
#     elseif panel == "b"
#         p =  reconstruct(p, EZ = SA[0, 4, 0], μN = 1.192, α = 0)
#     elseif panel == "c"
#         p = reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6, α = 2)
#     else
#         p =  reconstruct(p, EZ = SA[0, 4, 0], μN = 1.192, α = 2)
#     end

#     Wtest(Wlist, p, ϕ0, panel)
# end
 
# function Lntest(Lnlist, p, ϕ0, panel; nev = 4, kw...)
#     elist = zeros(Float64, nev, length(Lnlist))
#     for i in 1:length(Lnlist)
#         p = reconstruct(p, Ln = Lnlist[i])
#         ph = rectangle_squid(p; kw...)
#         s = spectrum(ph(phi = ϕ0), method = ArpackPackage(sigma = 1e-7, nev = nev))
#         elist[:,i] = s.energies
    
#         println(s.energies, " Ln: ",Lnlist[i])
#     end
#     return Lnlist, elist
# end

# function testLn(W0, ϕ0, panel, Lnlist )
     
#     p = Params(W = W0, Ls = 20, Ws = 20, scale = 40, λ = 5, Δ = 1, d = 0, τ = 1);
#     if panel == "a"
#         p = reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6, α = 0)
#     elseif panel == "b"
#         p =  reconstruct(p, EZ = SA[0, 4, 0], μN = 1.192, α = 0)
#     elseif panel == "c"
#         p = reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6, α = 2)
#     else
#         p =  reconstruct(p, EZ = SA[0, 4, 0], μN = 1.192, α = 2)
#     end

#     Lntest(Lnlist, p, ϕ0, panel)
# end
