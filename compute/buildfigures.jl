using CairoMakie, VegaLite, Colors, LaTeXStrings, MakieTeX
using MBSinBLG

# DATA PATHS
pathfig2 = "MBSinBLG/data/fig2"

pathfig3 = "MBSinBLG/data/fig3"

pathfig4a = "MBSinBLG/data/fig4/panela"
pathfig4b = "MBSinBLG/data/fig4/panelb"
pathfig4c = "MBSinBLG/data/fig4/panela"
pathfig4d = "MBSinBLG/data/fig4/panelb"

pathfig5a = "MBSinBLG/data/fig5/panela"
pathfig5b = "MBSinBLG/data/fig5/panelb"
pathfig5c = "MBSinBLG/data/fig5/panelc"
pathfig5d = "MBSinBLG/data/fig5/paneld"

pathfig6a1 = "MBSinBLG/data/fig6/panela_red"
pathfig6a2 = "MBSinBLG/data/fig6/panela_blue"
pathfig6b  = "MBSinBLG/data/fig6/panelb"
pathfig6c1 = "MBSinBLG/data/fig6/panelc_red"
pathfig6c2 = "MBSinBLG/data/fig6/panelc_blue"
pathfig6d  = "MBSinBLG/data/fig6/paneld"
pathfig6e  = "MBSinBLG/data/fig6/panele"
pathfig6f  = "MBSinBLG/data/fig6/panelf"

pathfig7 = "MBSinBLG/data/fig7"

# fig2
fig2plot(pathfig2)

# fig3
fig3plot(pathfig3)

# fig4
figab = fig4plot(pathfig4a, pathfig4b) #fig a and b
figc = ldosonlattice(pathfig4c, "reds", Ln = 1000)
figd = ldosonlattice(pathfig4d, "reds", Ln = 1000)

# fig5
fig5plot(pathfig5a, pathfig5b, pathfig5c, pathfig5d)

# fig6
figa1 = ldosonlattice(pathfig6a1, "reds", angle = 2)
figa2 = ldosonlattice(pathfig6a2, "blues", angle = 2)
figb  = ldosonlattice(pathfig6b, "reds", angle = 2)
figc1 = ldosonlattice(pathfig6c1, "reds", angle = 2)
figc2 = ldosonlattice(pathfig6c2, "blues", angle = 2) 
figd  = ldosonlattice(pathfig6d, "reds", angle = 2)
figef = fig6plot(pathfig6e, pathfig6f)

# fig7
fig7plot(pathfig7)
