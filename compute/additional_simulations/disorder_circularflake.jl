
function buildham(LN0; disorder = false)
    p = Params(Ln = LN0, Ls = 0, scale = 40, U = 1e-7, λ = 20, α = 0, EZ = SA[0, 0, 0], 
       Δ = 0., μN = 0.); #presets
    h = MBSinBLG.circular_hamiltonian(p,LN0/3,LN0/2, disorder = disorder) # builds a doughnut-shaped system
    return h
end

function psicircular(sp, itstart, itend)
    psi = zeros(Int(size(h,1)*2));
    for i in itstart:itend
        psi .+=  abs.(sp.states[:,i]).^2;
    end
    return psi
end

#RUN

## Compute
h = buildham(1152, disorder = false)
sp1  = spectrum(h, method = ArpackPackage(sigma = 0.00001im, nev = 100))
sp2 = dosKPM(h, order = 5000, resolution = 3, ket = randomkets(5, maporbitals = true) )

h = buildham(1152, disorder = true)
sp1d  = spectrum(h, method = ArpackPackage(sigma = 0.00001im, nev = 100))
sp2d = dosKPM(h, order = 5000, resolution = 3, ket = randomkets(5, maporbitals = true) )



## Plot spin currents and DOS
# without Anderson disorder
MBSinBLG.ldosonlatticeup(psicircular(sp1,1, 100), h)
MBSinBLG.ldosonlatticedown(psicircular(sp1,1, 100), h)
plot(sp2[1], sp2[2], xlims = (-15,15))

#with anderson disorder
MBSinBLG.ldosonlatticeup(psicircular(sp1d,1, 100), h)
MBSinBLG.ldosonlatticedown(psicircular(sp1d,1, 100), h)
plot(sp2d[1], sp2d[2], xlims = (-15,15))