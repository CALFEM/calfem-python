# -*- coding: utf-8 -*-
#
# example exs_beam1
# ----------------------------------------------------------------
# PURPOSE
#    Analysis of a simply supported beam.
# ----------------------------------------------------------------

# REFERENCES
#     Ola Dahlblom 2015-11-13
#     Ola Dahlblom 2019-12-11
#     Ola Dahlblom 2022-07-11
# ----------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import calfem.core as cfc
import calfem.utils as cfu
import calfem.vis_mpl as cfv

# ----- Topology -------------------------------------------------

edof = np.array([
     [1, 2, 3, 4], 
     [3, 4, 5, 6]
])

# ----- Stiffness matrix K and load vector f ---------------------

K = np.array(np.zeros((6, 6)))
f = np.array(np.zeros((6, 1)))
f[2] = -10e3

# ----- Element stiffness and element load matrices  -------------

E = 210e9
I = 2510e-8

ep = np.array([E, I])
ex1 = np.array([0, 3])
ex2 = np.array([3, 9])
eq1 = np.array([0])
eq2 = np.array([0])

Ke1 = cfc.beam1e(ex1, ep)
Ke2 = cfc.beam1e(ex2, ep)

# ----- Assemble Ke into K ---------------------------------------

cfc.assem(edof[0, :], K, Ke1)
cfc.assem(edof[1, :], K, Ke2)

# ----- Solve the system of equations and compute reactions ------

bc = np.array([1, 5])
a, r = cfc.solveq(K, f, bc)

cfu.disp_array(a, ["a"])
cfu.disp_array(r, ["r"])

# ----- Section forces -------------------------------------------

ed = cfc.extract_ed(edof, a)

es1, ed1, ec1 = cfc.beam1s(ex1, ep, ed[0, :], eq1, nep=4)
es2, ed2, ec2 = cfc.beam1s(ex2, ep, ed[1, :], eq2, nep=7)

cfu.disp_h2("es1")
cfu.disp_array(es1, ["V1", "M1"])
cfu.disp_h2("ed1")
cfu.disp_array(ed1, ["v1"])
cfu.disp_h2("ec1")
cfu.disp_array(ec1, ["x1"])
cfu.disp_h2("es2")
cfu.disp_array(es2, ["V2", "M2"])
cfu.disp_h2("ed2")
cfu.disp_array(ed2, ["v2"])
cfu.disp_h2("ec2")
cfu.disp_array(ec2, ["x2"])

# ----- Draw deformed beam ---------------------------------------

cfv.figure(1)
plt.plot([0, 9], [0, 0], color=(0.8, 0.8, 0.8))
plt.plot(
    np.concatenate(([0], ec1[:, 0], 3 + ec2[:, 0], [9]), 0),
    np.concatenate(([0], ed1[:, 0], ed2[:, 0], [0]), 0),
    color=(0.0, 0.0, 0.0),
)
cfv.title("displacements")

# ----- Draw shear force diagram----------------------------------

cfv.figure(2)
plt.plot([0, 9], [0, 0], color=(0.8, 0.8, 0.8))
plt.plot(
    np.concatenate(([0], ec1[:, 0], 3 + ec2[:, 0], [9]), 0),
    -np.concatenate(([0], es1[:, 0], es2[:, 0], [0]), 0) / 1000,
    color=(0.0, 0.0, 0.0),
)
cfv.title("shear force")

# ----- Draw moment diagram----------------------------------

cfv.figure(3)
plt.plot([0, 9], [0, 0], color=(0.8, 0.8, 0.8))
plt.plot(
    np.concatenate(([0], ec1[:, 0], 3 + ec2[:, 0], [9]), 0),
    -np.concatenate(([0], es1[:, 1], es2[:, 1], [0]), 0) / 1000,
    color=(0.0, 0.0, 0.0),
)
cfv.title("bending moment")

cfv.show_and_wait()
