# -*- coding: utf-8 -*-
#
# example exs_beam2
# ----------------------------------------------------------------
# PURPOSE
#    Analysis of a plane frame.
# ----------------------------------------------------------------

# REFERENCES
#     Göran Sandberg 94-03-08
#     Karl-Gunnar Olsson 95-09-28
#     Anders Olsson 99-03-01
#     Ola Dahlblom 2004-09-14
#     Ola Dahlblom 2019-12-16
#     Ola Dahlblom 2023-02-03
# ----------------------------------------------------------------

import numpy as np
import calfem.core as cfc
import calfem.utils as cfu
import calfem.vis_mpl as cfv

# ----- Topology -------------------------------------------------

edof = np.array([
    [4, 5, 6, 1, 2, 3], 
    [7, 8, 9, 10, 11, 12], 
    [4, 5, 6, 7, 8, 9]
])

# ----- Stiffness matrix K and load vector f ---------------------

K = np.array(np.zeros((12, 12)))
f = np.array(np.zeros((12, 1)))
f[3] = 2.0e3

# ----- Element stiffness and element load matrices  -------------

E = 200.0e9
A1 = 2.0e-3
A2 = 6.0e-3
I1 = 1.6e-5
I2 = 5.4e-5

ep1 = np.array([E, A1, I1])
ep3 = np.array([E, A2, I2])
ex1 = np.array([0, 0])
ex2 = np.array([6, 6])
ex3 = np.array([0, 6])
ey1 = np.array([4, 0])
ey2 = np.array([4, 0])
ey3 = np.array([4, 4])
eq1 = np.array([0, 0])
eq2 = np.array([0, 0])
eq3 = np.array([0, -10e3])

Ke1 = cfc.beam2e(ex1, ey1, ep1)
Ke2 = cfc.beam2e(ex2, ey2, ep1)
Ke3, fe3 = cfc.beam2e(ex3, ey3, ep3, eq3)

# ----- Assemble Ke into K ---------------------------------------

K = cfc.assem(edof[0, :], K, Ke1)
K = cfc.assem(edof[1, :], K, Ke2)
K, f = cfc.assem(edof[2, :], K, Ke3, f, fe3)

# ----- Solve the system of equations and compute reactions ------

bc = np.array([1, 2, 3, 10, 11])
a, r = cfc.solveq(K, f, bc)

cfu.disp_array(a, ["a"])
cfu.disp_array(r, ["r"])

# ----- Section forces -------------------------------------------

ed = cfc.extract_ed(edof, a)

es1, edi1, ec1 = cfc.beam2s(ex1, ey1, ep1, ed[0, :], eq1, nep=21)
es2, edi2, ec2 = cfc.beam2s(ex2, ey2, ep1, ed[1, :], eq2, nep=21)
es3, edi3, ec3 = cfc.beam2s(ex3, ey3, ep3, ed[2, :], eq3, nep=21)

cfu.disp_h2("es1")
cfu.disp_array(es1, ["N", "Vy", "Mz"])
cfu.disp_h2("edi1")
cfu.disp_array(edi1, ["u1", "v1"])
cfu.disp_h2("es2")
cfu.disp_array(es2, ["N", "Vy", "Mz"])
cfu.disp_h2("edi2")
cfu.disp_array(edi2, ["u1", "v1"])
cfu.disp_h2("es3")
cfu.disp_array(es3, ["N", "Vy", "Mz"])
cfu.disp_h2("edi3")
cfu.disp_array(edi3, ["u1", "v1"])

# ----- Draw deformed frame ---------------------------------------

plotpar = [2, 1, 0]
sfac = cfv.scalfact2(ex3, ey3, edi3, 0.1)
print("sfac=")
print(sfac)

cfv.figure(1)
cfv.eldraw2(ex1, ey1, plotpar)
cfv.eldraw2(ex2, ey2, plotpar)
cfv.eldraw2(ex3, ey3, plotpar)

plotpar = [1, 2, 1]
cfv.dispbeam2(ex1, ey1, edi1, plotpar, sfac)
cfv.dispbeam2(ex2, ey2, edi2, plotpar, sfac)
cfv.dispbeam2(ex3, ey3, edi3, plotpar, sfac)
cfv.axis([-1.5, 7.5, -0.5, 5.5])
plotpar1 = 2
cfv.scalgraph2(sfac, [1e-2, 0.5, 0], plotpar1)
cfv.title("Displacements")

# ----- Draw normal force diagram --------------------------------

plotpar = [2, 1]
sfac = cfv.scalfact2(ex1, ey1, es1[:, 0], 0.2)
cfv.figure(2)
cfv.secforce2(ex1, ey1, es1[:, 0], plotpar, sfac)
cfv.secforce2(ex2, ey2, es2[:, 0], plotpar, sfac)
cfv.secforce2(ex3, ey3, es3[:, 0], plotpar, sfac)
cfv.axis([-1.5, 7.5, -0.5, 5.5])
cfv.scalgraph2(sfac, [3e4, 1.5, 0], plotpar1)
cfv.title("Normal force")

# ----- Draw shear force diagram ---------------------------------

plotpar = [2, 1]
sfac = cfv.scalfact2(ex3, ey3, es3[:, 1], 0.2)
cfv.figure(3)
cfv.secforce2(ex1, ey1, es1[:, 1], plotpar, sfac)
cfv.secforce2(ex2, ey2, es2[:, 1], plotpar, sfac)
cfv.secforce2(ex3, ey3, es3[:, 1], plotpar, sfac)
cfv.axis([-1.5, 7.5, -0.5, 5.5])
cfv.scalgraph2(sfac, [3e4, 0.5, 0], plotpar1)
cfv.title("Shear force")

# ----- Draw moment diagram --------------------------------------

plotpar = [2, 1]
sfac = cfv.scalfact2(ex3, ey3, es3[:, 2], 0.2)
print("sfac=")
print(sfac)

cfv.figure(4)
cfv.secforce2(ex1, ey1, es1[:, 2], plotpar, sfac)
cfv.secforce2(ex2, ey2, es2[:, 2], plotpar, sfac)
cfv.secforce2(ex3, ey3, es3[:, 2], plotpar, sfac)
cfv.axis([-1.5, 7.5, -0.5, 5.5])
cfv.scalgraph2(sfac, [3e4, 0.5, 0], plotpar1)
cfv.title("Moment")

cfv.show_and_wait()
