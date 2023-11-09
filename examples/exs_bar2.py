# -*- coding: utf-8 -*-
#
# example exs3
# ----------------------------------------------------------------
# PURPOSE
#    Analysis of a plane truss.
# ----------------------------------------------------------------

# REFERENCES
#     Ola Dahlblom 2004-09-07
#     Jonas Lindemann 2009-01-25
#     Ola Dahlblom 2023-02-02
# ----------------------------------------------------------------

import numpy as np
import calfem.core as cfc
import calfem.utils as cfu
import calfem.vis_mpl as cfv

cfu.disp_h1("Analysis of a plane truss.")

# ----- Topology matrix Edof -------------------------------------

edof = np.array([
    [1, 2, 5, 6], 
    [5, 6, 7, 8], 
    [3, 4, 5, 6]
])

# ----- Stiffness matrix K and load vector f ---------------------

K = np.array(np.zeros((8, 8)))
f = np.array(np.zeros((8, 1)))

# ----- Element properties ---------------------------------------

E = 2.0e11
A1 = 6.0e-4
A2 = 3.0e-4
A3 = 10.0e-4
ep1 = [E, A1]
ep2 = [E, A2]
ep3 = [E, A3]

# ----- Element coordinates --------------------------------------

ex1 = np.array([0.0, 1.6])
ex2 = np.array([1.6, 1.6])
ex3 = np.array([0.0, 1.6])

ey1 = np.array([0.0, 0.0])
ey2 = np.array([0.0, 1.2])
ey3 = np.array([1.2, 0.0])

# ----- Element stiffness matrices  ------------------------------

Ke1 = cfc.bar2e(ex1, ey1, ep1)
Ke2 = cfc.bar2e(ex2, ey2, ep2)
Ke3 = cfc.bar2e(ex3, ey3, ep3)

# ----- Assemble Ke into K ---------------------------------------

K = cfc.assem(edof[0, :], K, Ke1)
K = cfc.assem(edof[1, :], K, Ke2)
K = cfc.assem(edof[2, :], K, Ke3)

cfu.disp_h2("Stiffness matrix K:")
cfu.disp_array(K)

# ----- Solve the system of equations ----------------------------

bc = np.array([1, 2, 3, 4, 7, 8])
f[5] = -80e3
a, r = cfc.solveq(K, f, bc)

cfu.disp_h2("Displacements a:")
cfu.disp_array(a)

cfu.disp_h2("Reaction forces r:")
cfu.disp_array(r)

# ----- Element forces -------------------------------------------

ed1 = cfc.extract_ed(edof[0, :], a)
N1 = cfc.bar2s(ex1, ey1, ep1, ed1)
ed2 = cfc.extract_ed(edof[1, :], a)
N2 = cfc.bar2s(ex2, ey2, ep2, ed2)
ed3 = cfc.extract_ed(edof[2, :], a)
N3 = cfc.bar2s(ex3, ey3, ep3, ed3)

cfu.disp_h2("Element forces r:")

print("N1 = ")
print(N1)
print("N2 = ")
print(N2)
print("N3 = ")
print(N3)

# ----- Draw deformed frame ---------------------------------------

plotpar = [2, 1, 0]
sfac = cfv.scalfact2(ex3, ey3, ed1, 0.1)
print("sfac=")
print(sfac)

cfv.figure(1)
cfv.eldraw2(ex1, ey1, plotpar)
cfv.eldraw2(ex2, ey2, plotpar)
cfv.eldraw2(ex3, ey3, plotpar)

plotpar = [1, 2, 1]
cfv.eldisp2(ex1, ey1, ed1, plotpar, sfac)
cfv.eldisp2(ex2, ey2, ed2, plotpar, sfac)
cfv.eldisp2(ex3, ey3, ed3, plotpar, sfac)
cfv.axis([-0.4, 2.0, -0.4, 1.4])
plotpar1 = 2
cfv.scalgraph2(sfac, [1e-3, 0, -0.3], plotpar1)
cfv.title("Displacements")

# ----- Draw normal force diagram --------------------------------

plotpar = [2, 1]
sfac = cfv.scalfact2(ex1, ey1, N2[:, 0], 0.1)
cfv.figure(2)
cfv.secforce2(ex1, ey1, N1[:, 0], plotpar, sfac)
cfv.secforce2(ex2, ey2, N2[:, 0], plotpar, sfac)
cfv.secforce2(ex3, ey3, N3[:, 0], plotpar, sfac)
cfv.axis([-0.4, 2.0, -0.4, 1.4])
cfv.scalgraph2(sfac, [5e4, 0, -0.3], plotpar1)
cfv.title("Normal force")

cfv.show_and_wait()
