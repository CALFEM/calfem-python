# -*- coding: utf-8 -*-

"""Example 05

This example shows how to make an unstructured 3D mesh (tetrahedron elements, which calfem cant actually use).
It also demonstrates how to do subplots and create two axes that are viewed from the same camera.
"""

import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_vtk as cfv

import vtk

# ---- Define geometry ------------------------------------------------------

g = cfg.geometry()

g.point([0, 0, 0], 0)
g.point([1, 0, 0], 1)
g.point([0, 1, 0], 2)
g.point([0, 1, 1], 3, el_size=0.1)
g.point([0.5, -0.3, 0], 4)
g.point([-0.3, 0.5, 0], 5)
g.point([0.75, 0.75, 0], 6)

g.spline([0, 4, 1])
g.spline([1, 6, 2])
g.spline([2, 5, 0])
g.spline([0, 3])
g.spline([3, 2])
g.spline([3, 1])

g.ruled_surface([0, 1, 2])
g.ruled_surface([0, 5, 3])
g.ruled_surface([1, 5, 4])
g.ruled_surface([2, 3, 4])

g.volume([0, 1, 2, 3])

# ---- Create mesh ----------------------------------------------------------

# Element type 4 is tetrahedron. (See user manual for more element types).

el_type = 4

# Degrees of freedom per node.

dofs_per_node = 1

# Create mesh

coords, edof, dofs, bdofs, elementmarkers = cfm.mesh(g, el_type, 0.3, dofs_per_node)

# coords, edof, dofs, bdofs, _ = cfm.mesh(
#        g, el_type, 0.3, dofs_per_node, gmsh_exec_path="D:\\vsmn20-software\\gmsh\gmsh.exe")


# ---- Visualise mesh -------------------------------------------------------

cfv.draw_mesh(coords, edof, el_type)
