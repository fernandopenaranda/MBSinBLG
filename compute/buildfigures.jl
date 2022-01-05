using MBSinBLG

########
#fig2
########
pathfig2 = ""
fig2plot()

########
#fig3
########
pathfig3 = ""
fig3plot()

########
#fig4
########
pathfig4a = "MBSinBLG/data/fig4/panela"
pathfig4b = "MBSinBLG/data/fig4/panelb"
fig4plot(pathfig4a, pathfig4b)
pathfig4c = "MBSinBLG/data/fig4/panela"
pathfig4d = "MBSinBLG/data/fig4/panelb"
psic = readfig4psi(pathfig4c)
psid = readfig4psi(pathfig4d)
h0 = rectangle_randombounds_sc(p, 0, 0)
figc = ldosonlattice(psic, h0, "reds")
figd = ldosonlattice(psic, h0, "reds")

########
#fig5
########
pathfig5a = "MBSinBLG/data/fig5/panela"
pathfig5b = "MBSinBLG/data/fig5/panelb"
pathfig5c = "MBSinBLG/data/fig5/panelc"
pathfig5d = "MBSinBLG/data/fig5/paneld"
fig5plot(pathfig5a, pathfig5b, pathfig5c, pathfig5d)

########
#fig6
########
pathfig6a = ""
pathfig6b1 = ""
pathfig6b2 = ""
pathfig6c = ""
pathfig6d1 = ""
pathfig6d2 = ""

psia = readfig4psi(pathfig6a)
psib1 = readfig4psi(pathfig6b1)
psib2 = readfig4psi(pathfig6b2)
psic = readfig4psi(pathfig6c)
psid1 = readfig4psi(pathfig6d1)
psid2 = readfig4psi(pathfig6d2)

h0 = rectangle_randombounds_sc(p, 2π/180, 0);
# a calculation of h0 with the same dimensions as in fig6 is required in order to 
# plot the wfs on top of the lattice.
figa = ldosonlattice(psia, h0, "reds")
figb1 = ldosonlattice(psib1, h0, "reds")
figb2 = ldosonlattice(psib2, h0, "blues")
figc = ldosonlattice(psic, h0, "reds")
figd1 = ldosonlattice(psid1, h0, "reds")
figd2 = ldosonlattice(psid2, h0, "blues")

pathfig6e = ""
pathfig6f = ""
fig6plot(pathfig6e, pathfig6f)

########
#fig7
########
pathfig7 = "MBSinBLG/data/fig7"
fig7plot(pathfig7)
