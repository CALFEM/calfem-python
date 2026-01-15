# -*- coding: utf-8 -*-
#
# example exs_bar2
# ----------------------------------------------------------------
# PURPOSE
#    Analysis of a plane truss.
# ----------------------------------------------------------------

# REFERENCES
#     Ola Dahlblom 2004-09-07
#     Jonas Lindemann 2009-01-25
#     Ola Dahlblom 2023-02-02
#     Henrik Danielsson 2026-01-08
# ----------------------------------------------------------------

# ----- Import necessary modules ---------------------------------

import numpy as np
import calfem.core as cfc
import calfem.utils as cfu
import calfem.vis_mpl as cfv

# ----- Topology matrix Edof -------------------------------------

edof = np.array([
    [1, 2, 5, 6], 
    [5, 6, 7, 8], 
    [3, 4, 5, 6]
])

# ----- Stiffness matrix K and load vector f ---------------------

K = np.array(np.zeros((8, 8)))
f = np.array(np.zeros((8, 1)))

f[5] = -80e3                       

# ----- Element properties ---------------------------------------

E = 2.0e11          # (N/m^2)
A1 = 6.0e-4         # (m^2)
A2 = 3.0e-4         # (m^2)
A3 = 10.0e-4        # (m^2)    
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

cfu.disp_h2("Ke1")
cfu.disp_array(Ke1, tablefmt='plain')
cfu.disp_h2("Ke2")
cfu.disp_array(Ke2, tablefmt='plain')
cfu.disp_h2("Ke3")
cfu.disp_array(Ke3, tablefmt='plain')

# ----- Assemble Ke into K ---------------------------------------

K = cfc.assem(edof[0, :], K, Ke1)
K = cfc.assem(edof[1, :], K, Ke2)
K = cfc.assem(edof[2, :], K, Ke3)

cfu.disp_h2("Stiffness matrix K (N/m):")
cfu.disp_array(K, tablefmt='plain')

# ----- Boundary conditions and nodal loads ----------------------

bc_dofs = np.array([1, 2, 3, 4, 7, 8])  
bc_vals = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

# ----- Solve the system of equations ----------------------------

a, r = cfc.solveq(K, f, bc_dofs, bc_vals)

cfu.disp_h2("Displacements a (m):")
cfu.disp_array(a)

cfu.disp_h2("Reaction forces r (N):")
cfu.disp_array(r)

# ----- Element forces -------------------------------------------

ed1 = cfc.extract_ed(edof[0, :], a)
ed2 = cfc.extract_ed(edof[1, :], a)
ed3 = cfc.extract_ed(edof[2, :], a)

es1 = cfc.bar2s(ex1, ey1, ep1, ed1)
es2 = cfc.bar2s(ex2, ey2, ep2, ed2)
es3 = cfc.bar2s(ex3, ey3, ep3, ed3)

cfu.disp_h2("Element normal forces (N):")

cfu.disp("N1")
cfu.disp_array(es1)
cfu.disp("N2")
cfu.disp_array(es2)
cfu.disp("N3")
cfu.disp_array(es3)

# ----- Draw deformed frame ---------------------------------------

plotpar1 = [2, 1, 0]
sfac = cfv.scalfact2(ex3, ey3, ed1, 0.1)
print("sfac =", sfac)

cfv.figure(1)
cfv.eldraw2(ex1, ey1, plotpar1)
cfv.eldraw2(ex2, ey2, plotpar1)
cfv.eldraw2(ex3, ey3, plotpar1)

plotpar2 = [1, 2, 1]
cfv.eldisp2(ex1, ey1, ed1, plotpar2, sfac)
cfv.eldisp2(ex2, ey2, ed2, plotpar2, sfac)
cfv.eldisp2(ex3, ey3, ed3, plotpar2, sfac)
cfv.axis([-0.4, 2.0, -0.4, 1.4])
plotpar3 = 2
cfv.scalgraph2(sfac, [1e-3, 0, -0.3], plotpar3)
cfv.title("Displacements")

# ----- Draw normal force diagram --------------------------------

plotpar = [2, 1]
sfac = cfv.scalfact2(ex1, ey1, es2[:, 0], 0.1)
cfv.figure(2)
cfv.secforce2(ex1, ey1, es1[:, 0], plotpar, sfac)
cfv.secforce2(ex2, ey2, es2[:, 0], plotpar, sfac)
cfv.secforce2(ex3, ey3, es3[:, 0], plotpar, sfac)
cfv.axis([-0.4, 2.0, -0.4, 1.4])
cfv.scalgraph2(sfac, [5e4, 0, -0.3], plotpar3)
cfv.title("Normal force")

cfv.show_and_wait()