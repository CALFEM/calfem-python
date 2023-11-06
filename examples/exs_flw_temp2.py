# -*- coding: utf-8 -*-
#
# example exs2
# ----------------------------------------------------------------
# PURPOSE
#    Analysis of one dimensional heat flow.
# ----------------------------------------------------------------

# REFERENCES
#     P-E Austrell 1994-03-08
#     K-G Olsson 1995-09-28
#     O Dahlblom 2004-09-07
#     J Lindemann 2009-01-25
# ----------------------------------------------------------------

import numpy as np
import calfem.core as cfc
import calfem.utils as cfu

# ----- Topology matrix Edof -------------------------------------

edof = np.array([
    [1, 2],
    [2, 3],
    [3, 4],
    [4, 5],
    [5, 6]
])

# ----- Stiffness matrix K and load vector f ---------------------

K = np.mat(np.zeros((6, 6)))
f = np.mat(np.zeros((6, 1)))
f[3] = 10.0

# ----- Element properties ---------------------------------------

ep1 = 25.0
ep2 = 24.3
ep3 = 0.4
ep4 = 17.0
ep5 = 7.7

# ----- Element stiffness matrices  ------------------------------

Ke1 = cfc.spring1e(ep1)
Ke2 = cfc.spring1e(ep2)
Ke3 = cfc.spring1e(ep3)
Ke4 = cfc.spring1e(ep4)
Ke5 = cfc.spring1e(ep5)

# ---- Assemble Ke into K ---------------------------------------

cfc.assem(edof[0, :], K, Ke1)
cfc.assem(edof[1, :], K, Ke2)
cfc.assem(edof[2, :], K, Ke3)
cfc.assem(edof[3, :], K, Ke4)
cfc.assem(edof[4, :], K, Ke5)

cfu.disp_h2("Stiffness matrix K:")
cfu.disp_array(K)

# ----- Solve the system of equations ----------------------------

bc = np.array([1, 6])
bcVal = np.array([-17.0, 20.0])
a, r = cfc.solveq(K, f, bc, bcVal)

cfu.disp_h2("Temperatures a:")
cfu.disp_array(a)

cfu.disp_h2("Reaction flows r:")
cfu.disp_array(r)

# ----- Element flows -------------------------------------------

ed1 = cfc.extract_ed(edof[0, :], a)
ed2 = cfc.extract_ed(edof[1, :], a)
ed3 = cfc.extract_ed(edof[2, :], a)
ed4 = cfc.extract_ed(edof[3, :], a)
ed5 = cfc.extract_ed(edof[4, :], a)

q1 = cfc.spring1s(ep1, ed1)
q2 = cfc.spring1s(ep2, ed2)
q3 = cfc.spring1s(ep3, ed3)
q4 = cfc.spring1s(ep4, ed4)
q5 = cfc.spring1s(ep5, ed5)

cfu.disp_h2("Element flows r:")

print("q1 = "+str(q1))
print("q2 = "+str(q2))
print("q3 = "+str(q3))
print("q4 = "+str(q4))
print("q5 = "+str(q5))
