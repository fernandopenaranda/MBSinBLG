using Parameters, StaticArrays, CSV, DataFrames

function runfig2()
    p = Params(Ln = 1440, Ls = 0, scale = 40, U = 1e-7, λ = 10, α = 2, EZ = SA[0, 0, 0], 
        Δ = 0., μN = 0.);
    p = reconstruct(p, EZ = SA[0, 0.001, 0])
    p_EZ = reconstruct(p, EZ = SA[0, 2.8, 0])
    mat = fig2bands([p, p, p_EZ, p_EZ], 
       ["armchair", "zigzag", "armchair", "zigzag"])
    writecsv(mat, [p, p, p_EZ, p_EZ], "fig2", ["1", "2", "3", "4"])
end

fig2bands(listp::Array, edgelist) =
    [fig2bands(listp[i], edgelist[i]) for i in 1:length(listp)]

function fig2bands(p, which)
    kpoints = 301
    numbands = 56
    if which == "armchair"
        kmin, kmax = (-0.15, 0.15)
        axis = (0, 1)
        h = nanoribbonN(p, axis);
        b = bandstructure(h, 
            cuboid((kmin, kmax); subticks = kpoints), 
            mapping = x -> 2pi*x, 
            method = ArpackPackage(nev=numbands, maxiter = 300, sigma=-0.00001));
    else
        kmin, kmax = (0.2, 0.8)
        axis = (1, 0)
        h = nanoribbonN(p, axis);
        b = bandstructure(h, cuboid((kmin, kmax), subticks = kpoints),  
            mapping = x -> 2pi*x, 
            method = ArpackPackage(nev=numbands, maxiter = 300, sigma=-0.00001))

    end
    karray = collect(range(kmin, kmax, length = kpoints))
    ϵs = zeros(ComplexF64, numbands, kpoints)
    for i in 1:length(karray)
        ϵs[:,i], _ = b.diag((karray[i],))
    end
    return karray, ϵs, b
end

