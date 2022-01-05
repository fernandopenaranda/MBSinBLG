function runfig7()
    p = Params(Ln = 400, W = 2000, Ls = 0, scale = 40, λ = 5, α = 0,   
       Δ = 1, d = 0, τ = 1);
    ϕlist = -2π:π/6:2π
    println("1")
    _, esa = cpr(ϕlist, reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6, α = 0), nev = 16)
    println("2")
    _, esb = cpr(ϕlist, reconstruct(p, EZ = SA[0, 4, 0], μN = 1.192, α = 0), nev = 16)
    println("3")
    _, esc = cpr(ϕlist, reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6, α = 1), nev = 16)
    println("4")
    ϕs, esd = cpr(ϕlist, reconstruct(p, EZ = SA[0, 4, 0], μN = 1.192, α = 1), nev = 16)
    savecprs("fig7", p, ϕs, esa, esb, esc, esd)    
end
