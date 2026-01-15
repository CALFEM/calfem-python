"""
Test script for exs_beam1 example.
Analysis of a simply supported beam.
"""

import numpy as np
import calfem.core as cfc
import calfem.utils as cfu
import matplotlib.pyplot as plt

np.set_printoptions(precision=4, linewidth=120)

# Element topology matrix (DOFs only, no element numbers)
edof = np.array([
    [1, 2, 3, 4],  # Element 1: DOFs 1,2,3,4 (node 1-2)
    [3, 4, 5, 6]   # Element 2: DOFs 3,4,5,6 (node 2-3)
])

cfu.disp_h1("Element topology matrix (edof):")
cfu.disp_array(edof, fmt=".0f", headers=["Element", "DOF1", "DOF2", "DOF3", "DOF4"])

# Initialize global system
K = np.zeros((6, 6))
f = np.zeros((6, 1))
f[2] = -10000  # Load at DOF 3 (index 2 in 0-based indexing)

cfu.disp("Global stiffness matrix initialized.")
cfu.disp_h1("Load vector:")
cfu.disp_array(f, headers=["f"])

# Material and geometric properties
E = 210e9      # Young's modulus [Pa]
I = 2510e-8    # Moment of inertia [m⁴]
ep = [E, I]    # Element properties

# Element coordinates [m]
ex1 = np.array([0, 3])  # Element 1: from x=0 to x=3
ex2 = np.array([3, 9])  # Element 2: from x=3 to x=9

# Compute element stiffness matrices
Ke1 = cfc.beam1e(ex1, ep)
cfu.disp_h1("Element 1 stiffness matrix:")
cfu.disp_array(Ke1)

Ke2 = cfc.beam1e(ex2, ep)
cfu.disp_h1("Element 2 stiffness matrix:")
cfu.disp_array(Ke2)

# Assemble global stiffness matrix
K = cfc.assem(edof[0], K, Ke1)
K = cfc.assem(edof[1], K, Ke2)

cfu.disp("")
cfu.disp("Global stiffness matrix assembled successfully")

# Boundary conditions (simply supported beam)
bc = np.array([
    [1, 0],  # DOF 1 = 0 (vertical displacement at left support)
    [5, 0]   # DOF 5 = 0 (vertical displacement at right support)
])

bc_dof = np.array([1, 5])
bc_value = np.array([0.0, 0.0])

# Solve the system
a, r = cfc.solveq(K, f, bc_dof, bc_value)
cfu.disp_h1("Displacements:")
cfu.disp_array(a)
cfu.disp_h1("Reaction forces [N]:")
cfu.disp_array(r)

# Extract element displacements
Ed = cfc.extract_ed(edof, a)

# Compute section forces and internal displacements
# For beam1s, we don't need eq parameter for point loads
es1, edi1, eci1 = cfc.beam1s(ex1, ep, Ed[0], nep=6)
es2, edi2, eci2 = cfc.beam1s(ex2, ep, Ed[1], nep=11)

cfu.disp_h1("Element 1 section forces [N, Nm]:")
cfu.disp_array(es1[:5])  # Show first 5 points
cfu.disp_h1("Element 2 section forces [N, Nm]:")
cfu.disp_array(es2[:5])  # Show first 5 points

