# -*- coding: utf-8 -*-

'''
Example 02
Creating geometry from B-Splines and circle arcs.
Also shows how to set ID numbers for geometry entities and how to specify element density. 
'''

import matplotlib.pyplot as plt
import matplotlib.collections
import numpy as np

import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_mpl as cfv

# ---- Define geometry ------------------------------------------------------

g = cfg.Geometry()

# Add points:
#  In this example we set the IDs manually.

g.point([-2,  0], ID=0)
# el_size determines the size of the elements near this point.
g.point([0,  1], ID=1, el_size=5)
# el_size is 1 by default. Larger number means less dense mesh.
g.point([1,  0], 2, el_size=5)
g.point([0, -2], 3)  # Size means the length of the sides of the elements.
g.point([0,  0], 4, el_size=5)
g.point([.5, .2], 5)
g.point([-.5, .5], 6)
g.point([-.7, -.5], 7)

# Add curves:

# The 3 points that define the circle arc are [start, center, end].
# The arc must be smaller than Pi.

g.circle([1, 4, 2], 2)

# BSplines are similar to Splines, but do not necessarily pass through the
# control points.

g.bspline([5, 6, 7, 5], 5)
g.bspline([1, 0, 3, 2], 4)

# Add surface:

g.surface([4, 2], [[5]])

# Markers do not have to be set when the curve is created. It can be done afterwards.
#  Set marker=80 for curves 2 and 4:

for curveID in [2, 4]:
    g.curveMarker(curveID, 80)

# ---- Generate mesh --------------------------------------------------------

mesh = cfm.GmshMesh(g)

# Element type 2 is triangle. (3 is quad. See user manual for more element types)

mesh.el_type = 3

# Degrees of freedom per node.

mesh.dofs_per_node = 2
mesh.el_size_factor = 0.05
# mesh.gmsh_exec_path = "D:\\vsmn20-software\\gmsh\gmsh.exe"

coords, edof, dofs, bdofs, elementmarkers = mesh.create()

# ---- Visualise mesh -------------------------------------------------------

# Hold left mouse button to pan.
# Hold right mouse button to zoom.

# Draw the geometry.

cfv.draw_geometry(
    g,
    label_curves=True,
    title="Example 2 - Geometry"
)

# New figure window

cfv.figure()

# Draws the mesh.

cfv.draw_mesh(
    coords=coords,
    edof=edof,
    dofs_per_node=mesh.dofs_per_node,
    el_type=mesh.el_type,
    filled=True,
    title="Example 2 - Mesh"
)

# Enter main loop

cfv.show_and_wait()
