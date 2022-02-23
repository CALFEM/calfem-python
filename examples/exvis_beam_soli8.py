# -*- coding: utf-8 -*-
"""
Example: Vedo visualisation

Example of using vedo to visualise results from a structured 3D problem
"""

import sys

import numpy as np

from vedo import *
import vtk

import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis as cfv
import calfem.core as cfc
import calfem.utils as cfu
import calfem.vis_vedo_utils as cfvu

# ---- Define geometry ------------------------------------------------------

print("Defining geometry...")

g = cfg.geometry()

# Add Points

l = 5.0
h = 0.5
w = 0.3

n_el_x = 5
n_el_y = 5
n_el_z = 50

marker_fixed_left = 45
marker_fixed_right = 46
marker_top = 47
right_support = 48

g.point([0, 0, 0], 0)
g.point([0.0, 0.0, w/2.0], 1)
g.point([0, 0, w], 2)
g.point([l, 0, w], 3)
g.point([l, 0, 0], 4, marker=11)  # Set some markers no reason.
g.point([0, h, 0], 5, marker=11)  # (markers can be given to points as well
# as curves and surfaces)
g.point([0, h, w], 6, marker=11)
g.point([l, h, w], 7)
g.point([l, h, 0], 8)

# Add splines

g.spline([0, 1, 2], 0, marker=33, el_on_curve=n_el_x)
g.spline([2, 3], 1, marker=23, el_on_curve=n_el_z)
g.spline([3, 4], 2, marker=right_support, el_on_curve=n_el_x)
g.spline([4, 0], 3, el_on_curve=n_el_z)
g.spline([0, 5], 4, el_on_curve=n_el_y)
g.spline([2, 6], 5, el_on_curve=n_el_y)
g.spline([3, 7], 6, el_on_curve=n_el_y)
g.spline([4, 8], 7, el_on_curve=n_el_y)
g.spline([5, 6], 8, el_on_curve=n_el_x)
g.spline([6, 7], 9, el_on_curve=n_el_z)
g.spline([7, 8], 10, el_on_curve=n_el_x)
g.spline([8, 5], 11, el_on_curve=n_el_z)

# Add surfaces

g.structuredSurface([0, 1, 2, 3], 0)
g.structuredSurface([8, 9, 10, 11], 1, marker=marker_top)
g.structuredSurface([0, 4, 8, 5], 2, marker=marker_fixed_left)
g.structuredSurface([1, 5, 9, 6], 3)
g.structuredSurface([2, 6, 10, 7], 4, marker=marker_fixed_right)
g.structuredSurface([3, 4, 11, 7], 5)

g.structuredVolume([0, 1, 2, 3, 4, 5], 0, marker=90)

# ---- Create mesh ----------------------------------------------------------

# Element type 5 is hexahedron. (See user manual for more element types)

el_type = 5
dofs_per_node = 3
el_nodes = 8

# Create mesh

print("Generating mesh...")
coords, edof, dofs, bdofs, elementmarkers = cfm.mesh(
    g, el_type, 1, dofs_per_node)

print("Extracting element coordinates...")
ex, ey, ez = cfc.coord_extract(edof, coords, dofs)

t = 0.2
v = 0.35
E = 2e9
ep = [3]
D = cfc.hooke(4, E, v)

# ---- Solve problem --------------------------------------------------------

print("Assembling stiffness matrix...")

n_dofs = np.size(dofs)
K = np.zeros((n_dofs, n_dofs))

for eltopo, elx, ely, elz, el_marker in zip(edof, ex, ey, ez, elementmarkers):
    Ke = cfc.soli8e(elx, ely, elz, [2], D)
    cfc.assem(eltopo, K, Ke)

bc = np.array([], 'i')
bcVal = np.array([], 'i')

bc, bcVal = cfu.apply_bc_3d(bdofs, bc, bcVal, marker_fixed_left, 0.0)
#bc, bcVal = cfu.apply_bc_3d(bdofs, bc, bcVal, marker_fixed_right, 0.0)
#bc, bcVal = cfu.apply_bc_3d(bdofs, bc, bcVal, right_support, 0.0, dimension=2)
#bc, bcVal = cfu.apply_bc_3d(bdofs, bc, bcVal, right_support, 0.0, dimension=3)

f = np.zeros([n_dofs, 1])
cfu.apply_force_total_3d(bdofs, f, marker_top, value=-10e5, dimension=2)

print("Solving equation system")
a, r = cfc.solveq(K, f, bc, bcVal)

print("Extracting element displacements")
ed = cfc.extract_eldisp(edof, a)

# ---- Extract element forces -----------------------------------------------

print("Calculating element forces")

# We can visualise stresses per element or per node
#
# sigv - von mises stresses per element
# sigv_hex - von mises stresses per element node
# sigv_el - all von mieses stresses att integration points
#
# cfvu.von_mises(...) calculates von mises stresses at
# integration points
#
# cfvu.sig_to_hex(...) extracts von mises stresses at nodes
# corresponding to element node

n_ip = ep[0]*ep[0]*ep[0]

sigv = np.zeros(ed.shape[0])
sigv_hex = np.zeros((ed.shape[0], el_nodes))
sigv_el = np.zeros((ed.shape[0], n_ip))

i = 0

for elx, ely, elz, eld in zip(ex, ey, ez, ed):
    et, es, eci = cfc.soli8s(elx, ely, elz, ep, D, eld)
    sigv[i] = np.mean(cfvu.von_mises_3d(es))
    sigv_hex[i] = cfvu.sigv_to_hex(cfvu.von_mises_3d(es))
    sigv_el[i] = cfvu.von_mises_3d(es)
    i += 1

print("Max mises stress = ", np.max(sigv_el))
print("Min mises stress = ", np.min(sigv_el))

# ---- Visualisation --------------------------------------------------------

print("Visualising results...")

nodes, topo, node_dofs, node_displ, node_scalars = cfvu.convert_to_node_topo(
    edof, ex, ey, ez, ed, sigv_hex, dofs_per_node=3)

npoint = nodes.shape[0]
nel = topo.shape[0]
nnd = topo.shape[1]

if nnd == 4:
    ct = vtk.VTK_TETRA
elif nnd == 8:
    ct = vtk.VTK_HEXAHEDRON
else:
    print("Topology not supported.")
    sys.exit(1)

celltypes = [ct] * nel

max_deflection = np.max(node_displ)
print("Max deflection = ", max_deflection)

new_coords = nodes + node_displ*0.01*l/max_deflection

ug_cell = UGrid([nodes, topo, celltypes])
ug_cell.points(new_coords)
ug_cell.celldata["cell_stress"] = sigv

ug_point = UGrid([nodes, topo, celltypes])
ug_point.points(new_coords)
ug_point.pointdata["node_stress"] = node_scalars

msh_cell = ug_cell.tomesh()
msh_cell.cmap("jet", "cell_stress", on="cells")
msh_cell.addScalarBar()

msh_point = ug_point.tomesh()
msh_point.cmap("jet", "node_stress", on="points")
msh_point.addScalarBar()

plt = Plotter(shape=(1, 2))
plt.show([msh_cell], at=0)
plt.show([msh_point], at=1, interactive=True)
