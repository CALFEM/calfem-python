# -*- coding: utf-8 -*-

"""Example 06

Solves a plane stress 2D problem using a structured mesh.
Shows how to draw von Mises effective stress as an element value with 
drawElementValues(). Shows use of GmshMesher attribute 'nodesOnCurve' 
(dictionary that says which nodes are on a given geometry curve)
"""

import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_mpl as cfv
import calfem.utils as cfu
import calfem.core as cfc
import numpy as np

from math import sqrt
import sys

cfu.enableLogging()

# ---- Define problem variables ---------------------------------------------

t = 0.2
v = 0.35
E = 2.1e9
ptype = 1
ep = [ptype, t]
D = cfc.hooke(ptype, E, v)

# ---- Define geometry ------------------------------------------------------

cfu.info("Creating geometry...")

g = cfg.geometry()

# Just a shorthand. We use this to make the circle arcs.

s2 = 1 / sqrt(2)

points = [
    [0, 3],
    [3, 3],
    [4 - s2, 3 - s2],
    [4, 2],  # 0-3
    [4 + s2, 3 - s2],
    [5, 3],
    [8, 3],
    [0, 0],  # 9-13
    [3, 0],
    [4 - s2, s2],
    [4, 1],
    [4 + s2, s2],  # 14-18
    [5, 0],
    [8, 0],
    [4, 3],
    [4, 0]
]

for xp, yp in points:
    g.point([xp * 0.1, yp * 0.1])




lines = [
    [0, 1],  # 0-1
    [5, 6],
    [6, 13],
    [13, 12],  # 4-7
    [8, 7],
    [7, 0],
]  

for s in lines:
    g.line(s)

# Points in circle arcs are [start, center, end]

circle_arcs = [
    [1, 14, 2],
    [2, 14, 3],
    [3, 14, 4],
    [4, 14, 5],
    [12, 15, 11],
    [11, 15, 10],
    [10, 15, 9],
    [9, 15, 8]
] 

for c in circle_arcs:
    g.circle(c)

# Create surfaces from the lines and arcs

g.surface([0, 6, 7, 8, 9, 1, 2, 3, 10, 11, 12, 13, 4, 5])

# Create markers

right_side = 20
left_side = 30

g.line_marker(2, right_side)  # Right side
g.line_marker(5, left_side)  # Right side

# ---- Create mesh ----------------------------------------------------------

cfu.info("Meshing geometry...")

# Create mesh

mesh = cfm.GmshMesh(geometry=g, return_boundary_elements=True)
mesh.el_type = 3
mesh.dofs_per_node = 2
mesh.el_size_factor = 0.01

coords, edof, dofs, bdofs, element_markers, b_elements = mesh.create()

# ---- Solve problem --------------------------------------------------------

cfu.info("Assembling system matrix...")

nDofs = np.size(dofs)
ex, ey = cfc.coordxtr(edof, coords, dofs)
K = np.zeros([nDofs, nDofs])

for eltopo, elx, ely in zip(edof, ex, ey):
    Ke = cfc.planqe(elx, ely, ep, D)
    cfc.assem(eltopo, K, Ke)

cfu.info("Solving equation system...")

f = np.zeros([nDofs, 1])

bc = np.array([], "i")
bc_val = np.array([], "f")

bc, bc_val = cfu.apply_bc(bdofs, bc, bc_val, left_side, 0.0, 0)

#cfu.apply_force(bdofs, f, right_side, 10e5, 1)
# boundaryElements, coords, dofs, F, marker, q
cfu.apply_traction_linear_element(
    b_elements, coords, dofs, f, right_side, [1e5, 0.0]
)
   

a, r = cfc.solveq(K, f, bc, bc_val)

cfu.info("Computing element forces...")

ed = cfc.extract_eldisp(edof, a)
von_mises = np.zeros((edof.shape[0]))

stress_table = np.zeros((edof.shape[0], 3))

# For each element:

for i in range(edof.shape[0]):
    # Determine element stresses and strains in the element.

    es, et = cfc.planqs(ex[i, :], ey[i, :], ep, D, ed[i, :])

    # Calc and append effective stress to list
    
    von_mises[i] = sqrt(pow(es[0], 2) - es[0] * es[1] + pow(es[1], 2) + 3 * es[2])

    stress_table[i, :] = es

# ---- Tabulate results -----------------------------------------------------

cfu.disp_array(stress_table, ["sigx", "sigy", "tauxy"])

# ---- Visualise results ----------------------------------------------------

cfu.info("Visualising...")

cfv.figure()
cfv.draw_geometry(
    g,
    draw_points=True,
    label_curves=True,
    label_points=True,
    title="Example 6a - Geometry",
)

cfv.figure()
cfv.draw_mesh(
    coords,
    edof,
    dofs_per_node=mesh.dofs_per_node,
    el_type=mesh.el_type,
    title="Example 6b - Meshing",
)

cfv.figure()
cfv.draw_element_values(
    von_mises,
    coords,
    edof,
    mesh.dofs_per_node,
    mesh.el_type,
    None,
    draw_elements=False,
    draw_undisplaced_mesh=False,
    title="Example 6c - Effective stress",
)

cfv.figure()
cfv.draw_displacements(
    a,
    coords,
    edof,
    mesh.dofs_per_node,
    mesh.el_type,
    draw_undisplaced_mesh=True,
    title="Example 6d - Displacements",
)

# Make use of attribute 'nodesOnCurve' in GmshMesher to draw some arrows on
# the right hand side of the mesh: --> currently not working

# rightSideNodes = set()

# 4 and 5 are the IDs of the curves where we applied the forces.

# for curveID in [4, 5]:
#     # Get the nodes, without duplicates.
#     rightSideNodes = rightSideNodes.union(set(mesh.nodesOnCurve[curveID]))

# for i in rightSideNodes:
#     # Position of the node with displacements.
#     x = coords[i, 0] + a[i * 2, 0]
#     y = coords[i, 1] + a[i * 2 + 1, 0]

#     # A poor man's force indicator. Could also use vv.plot()
#     cfv.add_text(r"$\rightarrow$", (x, y), color="g")

# Enter main loop
cfv.show_and_wait()

print("Done.")
