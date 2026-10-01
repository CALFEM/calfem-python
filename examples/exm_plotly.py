# -*- coding: utf-8 -*-

"""Interactive visualisation with plotly

Shows the drawing functions in calfem.vis_plotly using two problems on the
same geometry, a plate with a hole:

1. Plane stress, the plate in tension (4-node elements, plani4e/plani4s):
   geometry, mesh, effective stress and displacements.
2. Heat flow, a temperature difference between the short sides (3-node
   elements, flw2te/flw2ts): temperature contours and heat flux.

vis_plotly has the same functions as vis_mpl, so the plots can be drawn
with matplotlib instead by changing the import to calfem.vis_mpl.

Hover the plots to see element numbers and values, zoom by dragging and
reset the view by double clicking. The figures are written to HTML files
that open in the web browser, the location is printed. In Jupyter they are
shown inline (see exm_plotly.ipynb).

Requires plotly: pip install calfem-python[plotly]
"""

from pathlib import Path

import numpy as np

import calfem.core as cfc
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.utils as cfu
import calfem.vis_plotly as cfv

cfu.enableLogging()

# ---- Geometry -------------------------------------------------------------

L = 0.4             # Plate length [m]
H = 0.2             # Plate height [m]
r = 0.04            # Hole radius [m]

left_marker = 10
right_marker = 20

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

# ===========================================================================
# 1. Plane stress
# ===========================================================================

cfu.info("Solving plane stress problem...")

t = 0.01            # Thickness [m]
E = 210e9           # Young's modulus [Pa]
v = 0.3             # Poisson's ratio
F = 4e5             # Total tensile force on the right edge [N]

ep = [1, t, 2]      # Plane stress, thickness, 2x2 Gauss points
D = cfc.hooke(1, E, v)

mesh = cfm.GmshMesh(g)
mesh.el_type = 3            # 4-node quadrilaterals
mesh.dofs_per_node = 2
mesh.el_size_factor = 0.015

coords, edof, dofs, bdofs, elementmarkers = mesh.create()

n_dofs = np.size(dofs)
ex, ey = cfc.coordxtr(edof, coords, dofs)

K = np.zeros([n_dofs, n_dofs])
for eltopo, elx, ely in zip(edof, ex, ey):
    Ke, _ = cfc.plani4e(elx, ely, ep, D)
    cfc.assem(eltopo, K, Ke)

f = np.zeros([n_dofs, 1])
bc = np.array([], "i")
bc_val = np.array([], "f")
bc, bc_val = cfu.apply_bc(bdofs, bc, bc_val, left_marker, 0.0, 0)

n_right_nodes = len(bdofs[right_marker]) // mesh.dofs_per_node
cfu.apply_force(bdofs, f, right_marker, F/n_right_nodes, 1)

a, r = cfc.solveq(K, f, bc, bc_val)

ed = cfc.extract_eldisp(edof, a)
von_mises = np.zeros(edof.shape[0])
for i, (elx, ely, eld) in enumerate(zip(ex, ey, ed)):
    es, et = cfc.plani4s(elx, ely, ep, D, eld)
    sigx, sigy, tauxy = np.mean(es, axis=0)
    von_mises[i] = np.sqrt(sigx**2 - sigx*sigy + sigy**2 + 3*tauxy**2)

# ---- Visualise ------------------------------------------------------------

cfv.figure()
cfv.draw_geometry(g, title="Geometry - point and curve numbers [markers]")

cfv.figure()
cfv.draw_mesh(coords, edof, mesh.dofs_per_node, mesh.el_type,
              filled=True, title="Mesh")

cfv.figure()
cfv.draw_element_values(
    von_mises/1e6,
    coords,
    edof,
    mesh.dofs_per_node,
    mesh.el_type,
    displacements=a,
    magnfac=50.0,
    draw_elements=False,
    title="Effective stress (displacements x 50)",
)
cfv.colorbar("MPa")

cfv.figure()
cfv.draw_displacements(
    a,
    coords,
    edof,
    mesh.dofs_per_node,
    mesh.el_type,
    draw_undisplaced_mesh=True,
    title="Displacements (automatic scale)",
)

# ===========================================================================
# 2. Heat flow
# ===========================================================================

cfu.info("Solving heat flow problem...")

k = 50.0            # Thermal conductivity [W/(m K)]
D_heat = k*np.eye(2)
ep_heat = [t]

mesh = cfm.GmshMesh(g)
mesh.el_type = 2            # 3-node triangles
mesh.dofs_per_node = 1
mesh.el_size_factor = 0.02

coords, edof, dofs, bdofs, elementmarkers = mesh.create()

n_dofs = np.size(dofs)
ex, ey = cfc.coordxtr(edof, coords, dofs)

K = np.zeros([n_dofs, n_dofs])
for eltopo, elx, ely in zip(edof, ex, ey):
    Ke = cfc.flw2te(elx, ely, ep_heat, D_heat)
    cfc.assem(eltopo, K, Ke)

f = np.zeros([n_dofs, 1])
bc = np.array([], "i")
bc_val = np.array([], "f")
bc, bc_val = cfu.apply_bc(bdofs, bc, bc_val, left_marker, 100.0, 0)
bc, bc_val = cfu.apply_bc(bdofs, bc, bc_val, right_marker, 20.0, 0)

T, Q = cfc.solveq(K, f, bc, bc_val)

ed = cfc.extract_eldisp(edof, T)
flux = np.zeros((edof.shape[0], 2))
for i, (elx, ely, eld) in enumerate(zip(ex, ey, ed)):
    es, et = cfc.flw2ts(elx, ely, D_heat, eld)
    flux[i, :] = np.ravel(es)

# ---- Visualise ------------------------------------------------------------

cfv.figure()
cfv.draw_nodal_values_contourf(
    T,
    coords,
    edof,
    levels=16,
    dofs_per_node=mesh.dofs_per_node,
    el_type=mesh.el_type,
    draw_elements=True,
    title="Temperature",
    colorscale="RdBu_r",
)
cfv.colorbar("°C")

cfv.figure()
cfv.draw_nodal_values_contour(T, coords, edof, levels=16,
                              title="Isotherms", colorscale="RdBu_r")

cfv.figure()
cfv.draw_element_flux(
    flux,
    coords,
    edof,
    mesh.dofs_per_node,
    mesh.el_type,
    color_by_magnitude=True,
    draw_elements=True,
    title="Heat flux",
)
cfv.colorbar("W/m²")

# Save the last figure as an interactive web page. include_plotlyjs="cdn"
# loads plotly from the web instead of embedding it, giving a much smaller
# file.

flux_filename = Path(__file__).parent / "exm_plotly_flux.html"
cfv.save_figure(flux_filename, include_plotlyjs="cdn")
cfu.info(f"Heat flux figure saved to {flux_filename}")

cfv.show()

print("Done.")
