# -*- coding: utf-8 -*-

'''Example 05

This example shows how to make an unstructured 3D mesh (tetrahedron elements, which calfem cant actually use).
It also demonstrates how to do subplots and create two axes that are viewed from the same camera.
''' 

import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis as cfv

# ---- Define geometry ------------------------------------------------------

g = cfg.geometry()

g.point([0, 0, 0],      0)
g.point([1, 0, 0],      1)
g.point([0, 1, 0],      2)
g.point([0, 1, 1],      3, el_size=0.1)
g.point([0.5, -0.3, 0], 4)
g.point([-0.3, 0.5, 0], 5)
g.point([0.75, 0.75, 0],6)

g.spline([0,4,1])
g.spline([1,6,2])
g.spline([2,5,0])
g.spline([0,3])
g.spline([3,2])
g.spline([3,1])

g.ruledSurface([0,1,2])
g.ruledSurface([0,5,3])
g.ruledSurface([1,5,4])
g.ruledSurface([2,3,4])

g.volume([0,1,2,3])

# ---- Create mesh ----------------------------------------------------------

# Element type 4 is tetrahedron. (See user manual for more element types).

el_type = 4 

# Degrees of freedom per node.

dofs_per_node = 1 

# Create mesh

coords, edof, dofs, bdofs, elementmarkers = cfm.mesh(g, el_type, 0.3, dofs_per_node)

#coords, edof, dofs, bdofs, _ = cfm.mesh(
#        g, el_type, 0.3, dofs_per_node, gmsh_exec_path="D:\\vsmn20-software\\gmsh\gmsh.exe")


# ---- Visualise mesh -------------------------------------------------------

# Create two axes that are viewed from the same camera:

cfv.figure()
a1 = cfv.subplot(121)
a2 = cfv.subplot(122)
cam = cfv.camera3d()
a1.camera = a2.camera = cam

# Draw geometry and mesh

cfv.draw_geometry(g, axes=a1)
cfv.draw_mesh(coords=coords, edof=edof, dofs_per_node=dofs_per_node, el_type=el_type, filled=False, axes=a2)

# Enter main loop

cfv.show_and_wait()