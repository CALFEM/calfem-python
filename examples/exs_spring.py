# -*- coding: utf-8 -*-
#
# example exs1
# ----------------------------------------------------------------
# PURPOSE
#     Linear elastic spring analysis. Introduction to the basic
#     steps in the finite element method.
# ----------------------------------------------------------------

# REFERENCES
#     P-E Austrell 1994-03-08
#     K-G Olsson 1995-09-28
#     O Dahlblom 2004-09-06
#     J Lindemann 2009-01-25
# ----------------------------------------------------------------

# ----- import necesarry mooules

import numpy as np
import calfem.core as cfc
import calfem.utils as cfu

# ----- Topology matrix Edof

edof = np.array([
    [1, 2],      # element 1 between node 1 and 2
    [2, 3],      # element 2 between node 2 and 3
    [2, 3]       # element 3 between node 2 and 3
])

# ----- Stiffness matrix K and load vector f

K = np.zeros((3, 3))
f = np.zeros((3, 1))

# ----- Element stiffness matrices

k = 1500.
ep1 = k
ep2 = 2.*k
Ke1 = cfc.spring1e(ep1)
Ke2 = cfc.spring1e(ep2)

# ----- Assemble Ke into K

cfc.assem(edof[0, :], K, Ke2)
cfc.assem(edof[1, :], K, Ke1)
cfc.assem(edof[2, :], K, Ke2)

cfu.disp_h2("Stiffness matrix K:")
cfu.disp_array(K)

# f[1] corresponds to edof 2

f[1] = 100.0

# ----- Solve the system of equations

bc = np.array([1, 3])
a, r = cfc.solveq(K, f, bc)

cfu.disp_h2("Displacements a:")
cfu.disp_array(a)

cfu.disp_h2("Reaction forces r:")
cfu.disp_array(r)

# ----- Caculate element forces

ed1 = cfc.extract_ed(edof[0, :], a)
ed2 = cfc.extract_ed(edof[1, :], a)
ed3 = cfc.extract_ed(edof[2, :], a)

es1 = cfc.spring1s(ep2, ed1)
es2 = cfc.spring1s(ep1, ed2)
es3 = cfc.spring1s(ep2, ed3)

cfu.disp_h2("Element forces N:")
print("N1 = "+str(es1))
print("N2 = "+str(es2))
print("N3 = "+str(es3))
