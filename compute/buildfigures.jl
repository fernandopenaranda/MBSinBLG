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
pathfig6 = ""
fig6plot()

########
#fig7
########
pathfig7 = "MBSinBLG/data/fig7"
fig7plot(pathfig7)
