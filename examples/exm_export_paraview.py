# -*- coding: utf-8 -*-

"""Export results to ParaView

Solves a plane stress problem, a plate with a hole in tension, using
4-node isoparametric elements (plani4e/plani4s) and writes the mesh,
displacements and stresses to a VTU file with calfem.io.write_mesh().

Open the file in ParaView (https://www.paraview.org) to visualise the
results, for example:

* Color by "von_mises" (element values) or "displacement" (nodal values).
* Filters > Warp By Vector with "displacement" to show the deformed shape.
* Filters > Cell Data to Point Data to get smooth stress contours.

Requires meshio: pip install calfem-python[io]
"""

from pathlib import Path

import numpy as np

import calfem.core as cfc
import calfem.geometry as cfg
import calfem.io as cfio
import calfem.mesh as cfm
import calfem.utils as cfu
import calfem.vis_mpl as cfv

cfu.enableLogging()

# ---- Define problem variables ---------------------------------------------

L = 0.4             # Plate length [m]
H = 0.2             # Plate height [m]
r = 0.04            # Hole radius [m]
t = 0.01            # Thickness [m]
E = 210e9           # Young's modulus [Pa]
v = 0.3             # Poisson's ratio
F = 4e5             # Total tensile force on the right edge [N]

ptype = 1           # Plane stress
ir = 2              # 2x2 Gauss integration
ep = [ptype, t, ir]
D = cfc.hooke(ptype, E, v)

left_marker = 10
right_marker = 20

# ---- Define geometry ------------------------------------------------------

cfu.info("Creating geometry...")

g = cfg.Geometry()

g.point([0.0, 0.0])         # 0
g.point([L, 0.0])           # 1
g.point([L, H])             # 2
g.point([0.0, H])           # 3

g.point([L/2, H/2])         # 4 - hole center
g.point([L/2 + r, H/2])     # 5
g.point([L/2, H/2 + r])     # 6
g.point([L/2 - r, H/2])     # 7
g.point([L/2, H/2 - r])     # 8

g.spline([0, 1])                        # 0
g.spline([1, 2], marker=right_marker)   # 1
g.spline([2, 3])                        # 2
g.spline([3, 0], marker=left_marker)    # 3

g.circle([5, 4, 6])                     # 4
g.circle([6, 4, 7])                     # 5
g.circle([7, 4, 8])                     # 6
g.circle([8, 4, 5])                     # 7

g.surface([0, 1, 2, 3], holes=[[4, 5, 6, 7]])

# ---- Create mesh ----------------------------------------------------------

cfu.info("Meshing geometry...")

mesh = cfm.GmshMesh(g)
mesh.el_type = 3            # 4-node quadrilaterals
mesh.dofs_per_node = 2
mesh.el_size_factor = 0.01

coords, edof, dofs, bdofs, elementmarkers = mesh.create()

# ---- Solve problem --------------------------------------------------------

cfu.info("Assembling system matrix...")

n_dofs = np.size(dofs)
n_el = edof.shape[0]
ex, ey = cfc.coordxtr(edof, coords, dofs)

K = np.zeros([n_dofs, n_dofs])

for eltopo, elx, ely in zip(edof, ex, ey):
    Ke, _ = cfc.plani4e(elx, ely, ep, D)
    cfc.assem(eltopo, K, Ke)

cfu.info("Solving equation system...")

f = np.zeros([n_dofs, 1])

bc = np.array([], "i")
bc_val = np.array([], "f")
bc, bc_val = cfu.apply_bc(bdofs, bc, bc_val, left_marker, 0.0, 0)

# Distribute the total force equally over the nodes on the right edge
n_right_nodes = len(bdofs[right_marker]) // mesh.dofs_per_node
cfu.apply_force(bdofs, f, right_marker, F/n_right_nodes, 1)

a, r = cfc.solveq(K, f, bc, bc_val)

# ---- Element stresses -----------------------------------------------------

cfu.info("Computing element stresses...")

ed = cfc.extract_eldisp(edof, a)

stress = np.zeros((n_el, 3))
von_mises = np.zeros(n_el)

for i in range(n_el):
    # Stresses in the Gauss points, averaged over the element
    es, et = cfc.plani4s(ex[i, :], ey[i, :], ep, D, ed[i, :])
    sigx, sigy, tauxy = np.mean(es, axis=0)

    stress[i, :] = [sigx, sigy, tauxy]
    von_mises[i] = np.sqrt(sigx**2 - sigx*sigy + sigy**2 + 3*tauxy**2)

cfu.info(f"Max effective stress: {von_mises.max()/1e6:.1f} MPa")

# ---- Export to ParaView ---------------------------------------------------

filename = Path(__file__).with_suffix(".vtu")

cfu.info(f"Writing {filename}...")

cfio.write_mesh(
    filename,
    coords,
    edof,
    mesh.dofs_per_node,
    mesh.el_type,
    point_data={
        "displacement": a,      # 2 dofs per node, written as 3D vectors
    },
    cell_data={
        "von_mises": von_mises,
        "sigx": stress[:, 0],
        "sigy": stress[:, 1],
        "tauxy": stress[:, 2],
    },
)

# The mesh can also be read back, e.g. to continue in another script

coords2, edof2, dofs2, el_type2 = cfio.read_mesh(filename, dofs_per_node=2)
cfu.info(f"Read back {edof2.shape[0]} elements of type {el_type2}.")

# ---- Visualise results ----------------------------------------------------

cfv.figure()
cfv.draw_element_values(
    von_mises,
    coords,
    edof,
    mesh.dofs_per_node,
    mesh.el_type,
    draw_elements=False,
    title="Effective stress",
)
cfv.colorbar()

cfv.show_and_wait()

print("Done.")
