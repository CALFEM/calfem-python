# -*- coding: utf-8 -*-
#
# example exs_bar2_la
# ----------------------------------------------------------------
# PURPOSE
#    Analysis of a plane truss using loops.
# ----------------------------------------------------------------

# REFERENCES
#     P-E Austrell 1994-03-08
#     K-G Olsson 1995-09-28
#     O Dahlblom 2004-08-31
#     J Lindemann 2009-01-25
#     O Dahlblom 2023-02-02
# ----------------------------------------------------------------

import numpy as np
import calfem.core as cfc
import calfem.utils as cfu

cfu.disp_h1("Analysis of a plane truss using loops")

# ----- Topology matrix Edof -------------------------------------

edof = np.array([
    [1, 2, 5, 6],
    [3, 4, 7, 8],
    [5, 6, 9, 10],
    [7, 8, 11, 12],
    [7, 8, 5, 6],
    [11, 12, 9, 10],
    [3, 4, 5, 6],
    [7, 8, 9, 10],
    [1, 2, 7, 8],
    [5, 6, 11, 12],
])

# ----- Stiffness matrix K and load vector f ---------------------

K = np.zeros([12, 12])
f = np.zeros([12, 1])
f[10] = 0.5e6 * np.sin(np.pi / 6)
f[11] = -0.5e6 * np.cos(np.pi / 6)

# ----- Element properties ---------------------------------------

A = 25.0e-4
E = 2.1e11
ep = [E, A]

# ----- Element coordinates --------------------------------------

ex = np.array([
    [0.0, 2.0],
    [0.0, 2.0],
    [2.0, 4.0],
    [2.0, 4.0],
    [2.0, 2.0],
    [4.0, 4.0],
    [0.0, 2.0],
    [2.0, 4.0],
    [0.0, 2.0],
    [2.0, 4.0],
])

ey = np.array([
    [2.0, 2.0],
    [0.0, 0.0],
    [2.0, 2.0],
    [0.0, 0.0],
    [0.0, 2.0],
    [0.0, 2.0],
    [0.0, 2.0],
    [0.0, 2.0],
    [2.0, 0.0],
    [2.0, 0.0],
])

# ----- Create element stiffness matrices Ke and assemble into K -

for elx, ely, eltopo in zip(ex, ey, edof):
    Ke = cfc.bar2e(elx, ely, ep)
    cfc.assem(eltopo, K, Ke)

cfu.disp_h2("Stiffness matrix K:")
cfu.disp_array(K)

# ----- Solve the system of equations ----------------------------

bc = np.array([1, 2, 3, 4])
a, r = cfc.solveq(K, f, bc)

cfu.disp_h2("Displacements a:")
cfu.disp_array(a)

cfu.disp_h2("Reaction forces r:")
cfu.disp_array(r)

# ----- Element forces -------------------------------------------

ed = cfc.extract_ed(edof, a)
N = np.zeros([edof.shape[0]])

cfu.disp_h2("Element forces:")

i = 0
for elx, ely, eld in zip(ex, ey, ed):
    es = cfc.bar2s(elx, ely, ep, eld)
    N[i] = es[0][0]
    print("N%d = %g" % (i + 1, N[i]))
    i += 1
