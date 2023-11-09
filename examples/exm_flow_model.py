# -*- coding: utf-8 -*-
#
# example exm0_mpl.py
# ----------------------------------------------------------------
# PURPOSE
#    Setup a finite element flow model using the mesh functions
#    in CALFEM.
# ----------------------------------------------------------------
#
# REFERENCES
#     J Lindemann  2021-12-29
# ----------------------------------------------------------------

# ----- Import needed modules ------------------------------------

import numpy as np
import calfem.core as cfc
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_mpl as cfv
import calfem.utils as cfu

# ----- Problem parameters ---------------------------------------

w = 100.0
h = 10.0
t = 1.0
d = h / 2

D = np.identity(2, "float")
ep = [1.0, 1]

# ----- Create geometry object -----------------------------------

g = cfg.Geometry()

g.point([0, 0])  # point 1
g.point([w, 0])  # point 2
g.point([w, h])  # point 3
g.point([w - w / 2 + t / 2, h])  # point 4
g.point([w - w / 2 + t / 2, h - d])  # point 5
g.point([w - w / 2 - t / 2, h - d])  # point 6
g.point([w - w / 2 - t / 2, h])  # point 7
g.point([0, h])  # point 8

# ----- Create lines between points ------------------------------

left_side = 80
right_side = 90

g.spline([0, 1])
g.spline([1, 2])
g.spline([2, 3], marker=left_side)  # marker just to name
g.spline([3, 4])
g.spline([4, 5])
g.spline([5, 6])
g.spline([6, 7], marker=right_side)
g.spline([7, 0])

# ----- Make surface area ----------------------------------------

g.surface([0, 1, 2, 3, 4, 5, 6, 7])

# ----- Mesh generation ------------------------------------------

el_type = 3  # quadrature element
dofs_per_node = 1  # 1 dof

# ----- Set mesh paramters ---------------------------------------

mesh = cfm.GmshMesh(g)
mesh.el_size_factor = 1.0
mesh.el_type = el_type
mesh.dofs_per_node = dofs_per_node

# ----- Create mesh ----------------------------------------------

coords, edof, dofs, bdofs, elementmarkers = mesh.create()

# ----- Assemble elements ----------------------------------------

nDofs = np.size(dofs)
ex, ey = cfc.coordxtr(edof, coords, dofs)
K = np.zeros([nDofs, nDofs])

for eltopo, elx, ely in zip(edof, ex, ey):
    Ke = cfc.flw2i4e(elx, ely, ep, D)
    cfc.assem(eltopo, K, Ke)

# ----- Force vector ---------------------------------------------

f = np.zeros([nDofs, 1])

# ----- Boundary conditions --------------------------------------

bc = np.array([], int)
bcVal = np.array([], int)

bc, bcVal = cfu.applybc(bdofs, bc, bcVal, left_side, 0.0)
bc, bcVal = cfu.applybc(bdofs, bc, bcVal, right_side, 10.0)

# ----- Solve equation system ------------------------------------

a, r = cfc.solveq(K, f, bc, bcVal)
ed = cfc.extractEldisp(edof, a)

# ----- Calculating element forces -------------------------------

maxFlow = []  # empty list to store flow

for i in range(edof.shape[0]):
    es, et, eci = cfc.flw2i4s(ex[i, :], ey[i, :], ep, D, ed[i, :])
    maxFlow.append(np.sqrt(pow(es[0, 0], 2) + pow(es[0, 1], 2)))

# ----- Visualize results ----------------------------------------

cfv.figure()
cfv.draw_geometry(g, title="Geometry")

cfv.figure()
cfv.draw_element_values(
    maxFlow, coords, edof, dofs_per_node, el_type, None, title="Max flows"
)

cfv.figure()
cfv.draw_nodal_values(a, coords, edof, dofs_per_node=dofs_per_node, el_type=el_type)

cfv.showAndWait()
