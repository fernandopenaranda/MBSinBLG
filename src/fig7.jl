function runfig7()
    p = Params(Ln = 400, W = 2000, Ls = 40, Ws = 40, scale = 40, λ = 5,
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


    p = Params(Ln = 800, W = 2000, Ls = 20, Ws = 20, scale = 40, λ = 5,
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

    p = Params(Ln = 200, W = 2000, Ls = 20, Ws = 20, scale = 40, λ = 5,
    Δ = 1, d = 0, τ = 1);

    ϕlist = -2π:π/2:2π
    println("1")
    _, esa = cpr(ϕlist, reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6, α = 0), nev = 16)
    println("2")
    _, esb = cpr(ϕlist, reconstruct(p, EZ = SA[0, 4, 0], μN = 1.192, α = 0), nev = 16)
    println("3")
    _, esc = cpr(ϕlist, reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6, α = 1), nev = 16)
    println("4")
    ϕs, esd = cpr(ϕlist, reconstruct(p, EZ = SA[0, 4, 0], μN = 1.192, α = 1), nev = 16)
    savecprs("fig7", p, ϕs, esa, esb, esc, esd)    
    

    p = Params(Ln = 200, W = 3000, Ls = 20, Ws = 20, scale = 40, λ = 5,
    Δ = 1, d = 0, τ = 1);

    ϕlist = -2π:π/2:2π
    println("1")
    _, esa = cpr(ϕlist, reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6, α = 0), nev = 16)
    println("2")
    _, esb = cpr(ϕlist, reconstruct(p, EZ = SA[0, 4, 0], μN = 1.192, α = 0), nev = 16)
    println("3")
    _, esc = cpr(ϕlist, reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6, α = 1), nev = 16)
    println("4")
    ϕs, esd = cpr(ϕlist, reconstruct(p, EZ = SA[0, 4, 0], μN = 1.192, α = 1), nev = 16)
    savecprs("fig7", p, ϕs, esa, esb, esc, esd)    

    p = Params(Ln = 350, W = 2500, Ls = 20, Ws = 20, scale = 40, λ = 5,
    Δ = 1, d = 0, τ = 1);

    ϕlist = -2π:π/2:2π
    println("1")
    _, esa = cpr(ϕlist, reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6, α = 0), nev = 16)
    println("2")
    _, esb = cpr(ϕlist, reconstruct(p, EZ = SA[0, 4, 0], μN = 1.192, α = 0), nev = 16)
    println("3")
    _, esc = cpr(ϕlist, reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6, α = 1), nev = 16)
    println("4")
    ϕs, esd = cpr(ϕlist, reconstruct(p, EZ = SA[0, 4, 0], μN = 1.192, α = 1), nev = 16)
    savecprs("fig7", p, ϕs, esa, esb, esc, esd)    

    p = Params(Ln = 350, W = 3000, Ls = 20, Ws = 20, scale = 40, λ = 5,
    Δ = 1, d = 0, τ = 1);

    ϕlist = -2π:π/2:2π
    println("1")
    _, esa = cpr(ϕlist, reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6, α = 0), nev = 16)
    println("2")
    _, esb = cpr(ϕlist, reconstruct(p, EZ = SA[0, 4, 0], μN = 1.192, α = 0), nev = 16)
    println("3")
    _, esc = cpr(ϕlist, reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6, α = 2), nev = 16)
    println("4")
    ϕs, esd = cpr(ϕlist, reconstruct(p, EZ = SA[0, 4, 0], μN = 1.192, α = 2), nev = 16)
    savecprs("fig7", p, ϕs, esa, esb, esc, esd)    
end


p = Params(Ln = 400, W = 2000, Ls = 20, Ws = 20,scale = 40, λ = 5, α = 0,   
Δ = 1, d = 0, τ = 1);
ϕlist = -2π:π/2:2π
println("1")
_, esa = cpr(ϕlist, reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6, α = 0), nev = 16)
println("2")
_, esb = cpr(ϕlist, reconstruct(p, EZ = SA[0, 4, 0], μN = 1.192, α = 0), nev = 16)
println("3")
_, esc = cpr(ϕlist, reconstruct(p, EZ = SA[0, 2, 0], μN = 0.6, α = 1), nev = 16)
println("4")
ϕs, esd = cpr(ϕlist, reconstruct(p, EZ = SA[0, 4, 0], μN = 1.192, α = 1), nev = 16)
savecprs("fig7", p, ϕs, esa, esb, esc, esd)    