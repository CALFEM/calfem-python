# -*- coding: utf-8 -*-
#
# example exs_beam2
# ----------------------------------------------------------------
# PURPOSE
#    Analysis of a combined beam and bar structure.
# ----------------------------------------------------------------

# REFERENCES
#     Ola Dahlblom 2015-11-16
#     Ola Dahlblom 2019-12-19
#     Ola Dahlblom 2023-02-02
#     Copyright (c)  Division of Structural Mechanics and
#                    Division of Solid Mechanics.
#                    Lund University
# ----------------------------------------------------------------

import numpy as np
import calfem.core as cfc
import calfem.utils as cfu
import calfem.vis_mpl as cfv

#np.set_printoptions(precision=3, suppress=True)

# ----- Topology -------------------------------------------------

edof1 = np.array([
    [1,  2,  3,  4,  5,  6],
    [4,  5,  6,  7,  8,  9],
    [7,  8,  9, 10, 11, 12]    
])

edof2 = np.array([
    [13, 14,  4,  5],
    [13, 14,  7,  8]    
])

# ----- Stiffness matrix K and load vector f ---------------------

K = np.array(np.zeros((14, 14)))
f = np.array(np.zeros((14, 1)))

# ----- Element stiffness and element load matrices  -------------

E = 200.e9
A1 = 4.e-3
I1 = 5.4e-5
A2 = 1.e-3

ep1 = np.array([E, A1, I1])
ep4 = np.array([E, A2])

eq1 = np.array([0, 0])
eq2 = np.array([0, -10e+3])

ex1 = np.array([0, 2])
ex2 = np.array([2, 4])
ex3 = np.array([4, 6])
ex4 = np.array([0, 2])
ex5 = np.array([0, 4])
ey1 = np.array([2, 2])
ey2 = np.array([2, 2])
ey3 = np.array([2, 2])
ey4 = np.array([0, 2])
ey5 = np.array([0, 2])

Ke1 = cfc.beam2e(ex1, ey1, ep1)
Ke2, fe2 = cfc.beam2e(ex2, ey2, ep1, eq2)
Ke3, fe3 = cfc.beam2e(ex3, ey3, ep1, eq2)
Ke4 = cfc.bar2e(ex4, ey4, ep4)
Ke5 = cfc.bar2e(ex5, ey5, ep4)

# ----- Assemble Ke into K ---------------------------------------

K = cfc.assem(edof1[0, :], K, Ke1)
K, f = cfc.assem(edof1[1, :], K, Ke2, f, fe2)
K, f = cfc.assem(edof1[2, :], K, Ke3, f, fe3)
K = cfc.assem(edof2[0, :], K, Ke4)
K = cfc.assem(edof2[1, :], K, Ke5)

# ----- Solve the system of equations and compute reactions ------

bc = np.array([1, 2, 3, 13, 14])
a, r = cfc.solveq(K, f, bc)

cfu.disp_h2("Displacements a:")
cfu.disp_array(a)

cfu.disp_h2("Reaction forces r:")
cfu.disp_array(r)

# ----- Section forces -------------------------------------------

ed1 = cfc.extract_ed(edof1, a)
ed2 = cfc.extract_ed(edof2, a)

es1, _, _ = cfc.beam2s(ex1, ey1, ep1, ed1[0, :], eq1, nep=11)
es2, _, _ = cfc.beam2s(ex2, ey2, ep1, ed1[1, :], eq2, nep=11)
es3, _, _ = cfc.beam2s(ex3, ey3, ep1, ed1[2, :], eq2, nep=11)
es4 = cfc.bar2s(ex4, ey4, ep4, ed2[0, :])
es5 = cfc.bar2s(ex5, ey5, ep4, ed2[1, :])

cfu.disp_h2("es1 = ")
cfu.disp_array(es1, headers=["N", "Q", "M"])
cfu.disp_h2("es2 = ")
cfu.disp_array(es2, headers=["N", "Q", "M"])
cfu.disp_h2("es3 = ")
cfu.disp_array(es3, headers=["N", "Q", "M"])
cfu.disp_h2("es4 = ")
cfu.disp_array(es4, headers=["N"])
cfu.disp_h2("es5 = ")
cfu.disp_array(es5, headers=["N"])

