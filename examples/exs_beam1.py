# -*- coding: utf-8 -*-
#
# example exs5
# ----------------------------------------------------------------
# PURPOSE
#    Analysis of a simply supported beam.
# ----------------------------------------------------------------

# REFERENCES
#     Göran Sandberg 94-03-08
#     Karl-Gunnar Olsson 95-09-28
#     Ola Dahlblom 2004-09-21
# ----------------------------------------------------------------

import numpy as np
import calfem.core as cfc
import calfem.vis_mpl as cfv

# ----- Topology -------------------------------------------------

edof = np.array([
    [1, 2, 3, 4],
    [3, 4, 5, 6],
])

# ----- Stiffness matrix K and load vector f ---------------------

K = np.mat(np.zeros((6, 6)))
f = np.mat(np.zeros((6, 1)))
f[2] = -10000.

# ----- Element stiffness matrices  ------------------------------

E = 2.1e11
I = 2510e-8
ep = np.array([E, I])
ex1 = np.array([0., 3.])
ex2 = np.array([3., 9.])
Ke1 = cfc.beam1e(ex1, ep)
Ke2 = cfc.beam1e(ex2, ep)

# ----- Assemble Ke into K ---------------------------------------

K = cfc.assem(edof[0,:], K, Ke1)
K = cfc.assem(edof[1,:], K, Ke2)

# ----- Solve the system of equations and compute support forces -

bc = np.array([1, 5])
(a, r) = cfc.solveq(K, f, bc)

# ----- Section forces -------------------------------------------

ed = cfc.extract_ed(edof, a)

es1, ed1, ec1 = cfc.beam1s(ex1, ep, ed[0, :], nep=10)
es2, ed2, ec2 = cfc.beam1s(ex2, ep, ed[1, :], nep=10)

es = np.vstack(([0,0],es1,es2,[0,0]))
eD = np.vstack((0,ed1,ed2,0))
ec = np.vstack((ex1[0],ec1+ex1[0],ec2+ex2[0],ex2[1]))

# ----- Draw deformed beam ----------------------------------------
cfv.figure(1,fig_size=(6,2.5))
cfv.plt.plot([0, 9],[0,0],'b',linewidth=0.5)
cfv.plt.plot(ec,eD,'b')
cfv.plt.axis([-1,10,-0.03, 0.01])
cfv.title("Displacements")

# ----- Draw shear force diagram ----------------------------------
cfv.figure(2,fig_size=(6,2.5))
cfv.plt.plot([0, 9],[0,0],'b',linewidth=0.5)
cfv.plt.plot(ec,es[:,0],'b')
cfv.plt.axis([-1,10,-8000, 5000])
cfv.title("Shear force")
cfv.gca().invert_yaxis()

# ----- Draw bending moment diagram -------------------------------
cfv.figure(3,fig_size=(6,2.5))
cfv.plt.plot([0, 9],[0,0],'b',linewidth=0.5)
cfv.plt.plot(ec,es[:,1],'b')
cfv.plt.axis([-1,10,-5000, 25000])
cfv.title("Bending moment")
cfv.gca().invert_yaxis()
cfv.showAndWait()

#------------------------ end -----------------------------------