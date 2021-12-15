using Quantica, ArnoldiMethod, Arpack

include("model.jl")
include("plots.jl")
include("random_interface.jl")

using CSV, Dates, DataFrames
BLAS.set_num_threads(15)

p = Params(Ln = 2000, Ls = 0, scale = 40, λ = 6, α = 0, EZ = SA[0, .5, 0],  
    μN = 1.5, Δ = 0.3, d = 0, τ = 1);

presets_fig = Fig3_presets(0.0, 0.0, 0, 4, false)
data = ldosonlattice_averaged_sc(p, presets_fig)
savepsi("fig3", p, data[1])

presets_fig = Fig3_presets(0.0, 0.0, 0, 4, true)
data = ldosonlattice_averaged_sc(p, presets_fig)
savepsi("fig3", p, data[1])

presets_fig = Fig3_presets(0.0, 1*π/180, 0, 4, true)
data = ldosonlattice_averaged_sc(p, presets_fig)
savepsi("fig3", p, data[1])

presets_fig = Fig3_presets(0.0, 1.5*π/180, 0, 4, true)
data = ldosonlattice_averaged_sc(p, presets_fig)
savepsi("fig3", p, data[1])


presets_fig = Fig3_presets(0.001, 1*π/180, 5, 4, true)
data = ldosonlattice_averaged_sc(p, presets_fig)
savepsi("fig3", p, data[1])

presets_fig = Fig3_presets(0.001, 1.5*π/180, 5, 4, true)
data = ldosonlattice_averaged_sc(p, presets_fig)
savepsi("fig3", p, data[1])

presets_fig = Fig3_presets(0.005, 1*π/180, 5, 4, true)
data = ldosonlattice_averaged_sc(p, presets_fig)
savepsi("fig3", p, data[1])

presets_fig = Fig3_presets(0.005, 1.5*π/180, 5, 4, true)
data = ldosonlattice_averaged_sc(p, presets_fig)
savepsi("fig3", p, data[1])

p = Params(Ln = 3000, Ls = 0, scale = 40, λ = 6, α = 0, EZ = SA[0, .5, 0],  
    μN = 1.5, Δ = 0.3, d = 0, τ = 1);

presets_fig = Fig3_presets(0.0, 0.0, 0, 4, false)
data = ldosonlattice_averaged_sc(p, presets_fig)
savepsi("fig3", p, data[1])

presets_fig = Fig3_presets(0.0, 0.0, 0, 4, true)
data = ldosonlattice_averaged_sc(p, presets_fig)
savepsi("fig3", p, data[1])

presets_fig = Fig3_presets(0.0, 1*π/180, 0, 4, true)
data = ldosonlattice_averaged_sc(p, presets_fig)
savepsi("fig3", p, data[1])

presets_fig = Fig3_presets(0.0, 1.5*π/180, 0, 4, true)
data = ldosonlattice_averaged_sc(p, presets_fig)
savepsi("fig3", p, data[1])


presets_fig = Fig3_presets(0.001, 1*π/180, 5, 4, true)
data = ldosonlattice_averaged_sc(p, presets_fig)
savepsi("fig3", p, data[1])

presets_fig = Fig3_presets(0.001, 1.5*π/180, 5, 4, true)
data = ldosonlattice_averaged_sc(p, presets_fig)
savepsi("fig3", p, data[1])

presets_fig = Fig3_presets(0.005, 1*π/180, 5, 4, true)
data = ldosonlattice_averaged_sc(p, presets_fig)
savepsi("fig3", p, data[1])

presets_fig = Fig3_presets(0.005, 1.5*π/180, 5, 4, true)
data = ldosonlattice_averaged_sc(p, presets_fig)
savepsi("fig3", p, data[1])


p = Params(Ln = 2000, Ls = 0, scale = 40, λ = 6, α = 0, EZ = SA[0, 1, 0],  
    μN = 1.5, Δ = 0.5, d = 0, τ = 1);

presets_fig = Fig3_presets(0.0, 0.0, 0, 4, false)
data = ldosonlattice_averaged_sc(p, presets_fig)
savepsi("fig3", p, data[1])

presets_fig = Fig3_presets(0.0, 0.0, 0, 4, true)
data = ldosonlattice_averaged_sc(p, presets_fig)
savepsi("fig3", p, data[1])

# presets_fig = Fig3_presets(0.0, 1*π/180, 0, 4, true)
# data = ldosonlattice_averaged_sc(p, presets_fig)
# savepsi("fig3", p, data[1])

# presets_fig = Fig3_presets(0.0, 1.5*π/180, 0, 4, true)
# data = ldosonlattice_averaged_sc(p, presets_fig)
# savepsi("fig3", p, data[1])


# presets_fig = Fig3_presets(0.001, 1*π/180, 5, 4, true)
# data = ldosonlattice_averaged_sc(p, presets_fig)
# savepsi("fig3", p, data[1])

# presets_fig = Fig3_presets(0.001, 1.5*π/180, 5, 4, true)
# data = ldosonlattice_averaged_sc(p, presets_fig)
# savepsi("fig3", p, data[1])

# presets_fig = Fig3_presets(0.005, 1*π/180, 5, 4, true)
# data = ldosonlattice_averaged_sc(p, presets_fig)
# savepsi("fig3", p, data[1])

# presets_fig = Fig3_presets(0.005, 1.5*π/180, 5, 4, true)
# data = ldosonlattice_averaged_sc(p, presets_fig)
# savepsi("fig3", p, data[1])


## save and load
#savespsi("fig3", p_3b, datafig_3b[1])
#psi = readfig3psi(datapath)
# angle = 0.0 
# h0 = rectangle_randombounds_sc(p_0, angle, 0, sidecontacts = true);
##


p_3a = Params(Ln = 1400, Ls = 0, scale = 40, λ = 4, α = 0, EZ = SA[0, 1, 0],  
        μN = 1.2, Δ = 1, d = 0, τ = 1);

presets_3a = Fig3_presets(0.0, 0.0, 0, 4)
datafig_3a = ldosonlattice_averaged_sc(p_3a, presets_3a)
fig_3a = ldosonlattice(datafig_3a[1], datafig_3a[2])

# p_3b = Params(Ln = 1400, Ls = 0, scale = 40, λ = 4, α = 0, EZ = SA[0, 4, 0],  
#             μN = 0.4, Δ = 1, d = 0, τ = 1);

# p_3b = Params(Ln = 600, Ls = 0, scale = 40, λ = 4, α = 0, EZ = SA[0, 4, 0],  
#             μN = 0.4, Δ = 1, d = 0, τ = 1);

p_3b = Params(Ln = 1000, Ls = 0, scale = 40, λ = 6, α = 0, EZ = SA[0, .5, 0],  
    μN = 1.5, Δ = 0.3, d = 0, τ = 1);

presets_3b = Fig3_presets(0.0, 0.0, 0, 4)
@time datafig_3b = ldosonlattice_averaged_sc(p_3b, presets_3b)
fig_3b = ldosonlattice(datafig_3b[1], datafig_3b[2])

p_3c = Params(Ln = 2000, Ls = 0, scale = 40, λ = 6, α = 0, EZ = SA[0, .5, 0],  
    μN = 1.5, Δ = 0.3, d = 0, τ = 1);

presets_3c = Fig3_presets(0.00, 1*π/180, 0, 4)

presets_3d = Fig3_presets(0.005, 1*π/180, 5, 4)
@time datafig_3c = ldosonlattice_averaged_sc(p_3c, presets_3c)
fig_3b = ldosonlattice(datafig_3b[1], datafig_3b[2])

# @everywhere BLAS.set_num_threads(10) 

p_3b = Params(Ln = 1400, Ls = 0, scale = 40, λ = 4, α = 0, EZ = SA[0, 4, 0],  
            μN = 0.4, Δ = 1, d = 0, τ = 1);

#############
# LDOS panels e-f
#############


##
indexa = 63202
indexb = 20599

indexb = 33755
##

ldosa = ldos(p_3a, indexa, presets_3a.η, presets_3a.angle)
ldosb = ldos(p_3b, indexb, presets_3b.η, presets_3b.angle)

# ldosa = ldos(palarge, "a", indexab, 0, order = 30000) 
# ldosb = ldos(pblarge, "b", indexab, 0, order = 30000)
# ldosc = ldos(pcdlarge, "c", indexc, 0, order = 30000)
# ldosd = ldos(pcdlarge, "d", indexd, pi/180, order = 30000)

adata = DataFrame(x = ldosa[1], y = ldosa[2])
bdata = DataFrame(x = ldosb[1], y = ldosb[2])
# cdata = DataFrame(x = ldosc[1], y = ldosc[2])
# ddata = DataFrame(x = ldosd[1], y = ldosd[2]) 

CSV.write("data/fig3_ldosa", adata; delim = '\t')
CSV.write("data/fig3_ldosb", bdata; delim = '\t')
# CSV.write("data/fig3_ldosc", cdata; delim = '\t')
# CSV.write("data/fig3_ldosd", ddata; delim = '\t')

## Plot panelrightfig3(q[1], dataa.y, datab.y, datac.y, datad.y)
#panelrightfig3(q[1], adata.y, bdata.y, cdata.y, ddata.y)


############################
# Additional checks (Majorana oscillations, disorder, rotation,
#BDI to class D transition....)
############################



# n majoranas!?

####### sweep in muN
pbtest2 = Params(Ln = 400, Ls = 0, scale = 40, λ = 6, α = 2, EZ = SA[0, 6, 0], 
    μN = 0.85, Δ = 2.5, d = 0, τ = 1); # splitted Majos
sp = spectrumsweep(pbtest2, 0.3:0.05:0.85)
spectrumvsmu(sp[1], sp[2])

pbtest2 = Params(Ln = 800, Ls = 0, scale = 40, λ = 6, α = 2, EZ = SA[0, 6, 0], 
    μN = 0.85, Δ = 2.5, d = 0, τ = 1); # splitted Majos
sp = spectrumsweep(pbtest2, 0.3:0.05:0.85)
spectrumvsmu(sp[1], sp[2])

pbtest2 = Params(Ln = 800, Ls = 0, scale = 40, λ = 6, α = 0, EZ = SA[0, 5.5, 0],  
μN = 0.85, Δ = 2.7, d = 0, τ = 1); #unsplitted Majos
sp = spectrumsweep(pbtest2, 0.3:0.05:0.85)

spectrumvsmu(sp[1], sp[2])
pbtest2 = Params(Ln = 800, Ls = 0, scale = 40, λ = 6, α = 0, EZ = SA[0, 5.5, 0],  
μN = 0.85, Δ = 2.5, d = 0, τ = 1); #crossover Majos
sp = spectrumsweep(pbtest2, 0.3:0.05:0.85)
spectrumvsmu(sp[1], sp[2])

pbtest2 = Params(Ln = 800, Ls = 0, scale = 40, λ = 10, α = 0, EZ = SA[0, 5.5, 0],  
μN = 0.85, Δ = 2.7, d = 0, τ = 1); #unsplitted Majos larger lambda
sp = spectrumsweep(pbtest2, 0.7:0.05:2.6)
spectrumvsmu(sp[1], sp[2], (-0.4,0.4))
######### sweep in λ
pbtest2 = Params(Ln = 800, Ls = 0, scale = 40, λ = 10, α = 0, EZ = SA[0, 5.5, 0],  
μN = 0.4, Δ = 2.5, d = 0, τ = 1);
sp = spectrumsweepλ(pbtest2, 0:0.5:10)
spectrumvsλ(sp[1], sp[2], (-0.4,0.4))


pbtest2 = Params(Ln = 400, Ls = 0, scale = 40, λ = 10, α = 0, EZ = SA[0, 5, 0],  
μN = 0.4, Δ = 2.5, d = 0, τ = 1);
sp = spectrumsweepλ(pbtest2, 0:0.5:10)
spectrumvsλ(sp[1], sp[2], (-0.4,0.4))
######### sweep in EZ
pbtest2 = Params(Ln = 400, Ls = 0, scale = 40, λ = 6, α = 2, EZ = SA[0, 6, 0], 
    μN = 0.85, Δ = 2.5, d = 0, τ = 1);
sp = spectrumsweepb(pbtest2, 0:0.5:10)
spectrumvsb(sp[1], sp[2], (-0.4,0.4))

#########  sweep in Rashba
pbtest2 = Params(Ln = 800, Ls = 0, scale = 40, λ = 6, α = 0, EZ = SA[0, 6, 0], 
    μN = 0.4, Δ = 2.5, d = 0, τ = 1); # splitted Majos
sp = spectrumsweepα(pbtest2, 0:0.25:4)
spectrumvsα(sp[1], sp[2], (-0.4,0.4))

pbtest2 = Params(Ln = 800, Ls = 0, scale = 40, λ = 6, α = 0, EZ = SA[0, 6, 0], 
    μN = 0.6, Δ = 2.5, d = 0, τ = 1); # splitted Majos
sp = spectrumsweepα(pbte
st2, 0:0.25:4)
spectrumvsα(sp[1], sp[2], (-0.4,0.4))
############## reducing lambda

# spectrum random model
# h = rectangle_randombounds_sc(pbtest2, 0, 0.01, sidecontacts = true);
# sp = spectrum(h, method = ArpackPackage(nev = 16, sigma = 1e-7))
# hist(real(sp.energies), bins = 100)



########## testing rotation of the 8MBS phase
pbtest2 = Params(Ln = 800, Ls = 0, scale = 40, λ = 6, α = 0, EZ = SA[0, 5.5, 0],  
    μN = 0.4, Δ = 2.7, d = 0, τ = 1);
sp = splittingvsrotation(pbtest2, 0:.5:5)
spectrumvsangle(sp[1], sp[2], (-0.4,0.4))

pbtest2 = Params(Ln = 400, Ls = 0, scale = 40, λ = 6, α = 0, EZ = SA[0, 5.5, 0],  
    μN = 0.4, Δ = 2.7, d = 0, τ = 1);
ldosonlattice_averaged_sc(pbtest2, 0.0, 1*π/180, nummodes = 8)


#### tests

p_3b = Params(Ln = 2000, Ls = 0, scale = 40, λ = 6, α = 0, EZ = SA[0, 0.01, 0],  
μN = 1.5, Δ = 1, d = 0, τ = 1);

sp = spectrumsweepb(p_3b, 0:0.5:5)
spectrumvsb(sp[1], sp[2], (-0.4,0.4))


p_3b = Params(Ln = 2000, Ls = 0, scale = 40, λ = 6, α = 0, EZ = SA[0, .5, 0],  
    μN = 1.5, Δ = 0.3, d = 0, τ = 1);

sp = spectrumsweepb(p_3b, 0:0.5:5)
spectrumvsb(sp[1], sp[2], (-0.4,0.4))


p_3b = Params(Ln = 2000, Ls = 0, scale = 40, λ = 6, α = 0, EZ = SA[0, 0.01, 0],  
μN = 1.5, Δ = 0.3, d = 0, τ = 1);

sp = spectrumsweepb(p_3b, 0:0.5:5)
spectrumvsb(sp[1], sp[2], (-0.4,0.4))

