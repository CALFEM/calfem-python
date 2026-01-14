# -*- coding: utf-8 -*-
#
# example exs_spring
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
#     H Danielsson 2026-01-08
# ----------------------------------------------------------------

# ----- Import necessary modules ---------------------------------

import numpy as np
import calfem.core as cfc
import calfem.utils as cfu

# ----- Topology matrix Edof -------------------------------------

edof = np.array([
    [1, 2],                 # element 1 between nodes 1 and 2
    [2, 3],                 # element 2 between nodes 2 and 3
    [2, 3]                  # element 3 between nodes 2 and 3
])

# ----- Stiffness matrix K and load vector f ---------------------

K = np.zeros((3, 3))        # empty global stiffness matrix
f = np.zeros((3, 1))        # empty global load vector

f[1] = 100.0                # (N), f[1] corresponds to dof 2


# ----- Element stiffness matrices Ke ----------------------------

k = 1500.                   # (N/m)
ep1 = 2.*k                  # spring stiffness, element 1
ep2 = k                     # spring stiffness, element 2
ep3 = 2.*k                  # spring stiffness, element 3
Ke1 = cfc.spring1e(ep1)     # element stiffness matrix, element 1
Ke2 = cfc.spring1e(ep2)     # element stiffness matrix, element 2
Ke3 = cfc.spring1e(ep3)     # element stiffness matrix, element 3

cfu.disp_h2("Ke1")
cfu.disp_array(Ke1, tablefmt='plain')
cfu.disp_h2("Ke2")
cfu.disp_array(Ke2, tablefmt='plain')
cfu.disp_h2("Ke3")
cfu.disp_array(Ke3, tablefmt='plain')

# ----- Assemble Ke into K ---------------------------------------

cfc.assem(edof[0, :], K, Ke1)   # assemble element stiffness matrix 1
cfc.assem(edof[1, :], K, Ke2)   # assemble element stiffness matrix 2
cfc.assem(edof[2, :], K, Ke3)   # assemble element stiffness matrix 3

cfu.disp_h2("Stiffness matrix K (N/m):")
cfu.disp_array(K, tablefmt='plain')

# ----- Boundary conditions and nodal loads ----------------------

bc_dofs = np.array([1, 3])       # dofs with prescribed displacment (0)
bc_vals = np.array([0.0, 0.0])   # prescribed displacments

# ----- Solve the system of equations ----------------------------

a, r = cfc.solveq(K, f, bc_dofs, bc_vals)

cfu.disp_h2("Displacements a (m):")
cfu.disp_array(a, tablefmt='plain')

cfu.disp_h2("Reaction forces r (N):")
cfu.disp_array(r, tablefmt='plain')

# ----- Caculate element forces ----------------------------------

ed1 = cfc.extract_ed(edof[0, :], a) # Nodal displacements, element 1
ed2 = cfc.extract_ed(edof[1, :], a) # Nodal displacements, element 2
ed3 = cfc.extract_ed(edof[2, :], a) # Nodal displacements, element 3

es1 = cfc.spring1s(ep1, ed1)        # Element force, element 1
es2 = cfc.spring1s(ep2, ed2)        # Element force, element 2
es3 = cfc.spring1s(ep3, ed3)        # Element force, element 3

cfu.disp_h2("Element forces (N):")
print("N1 = "+str(es1))
print("N2 = "+str(es2))
print("N3 = "+str(es3))