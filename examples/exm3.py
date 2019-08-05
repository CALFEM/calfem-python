# -*- coding: utf-8 -*-

'''Example 03

Shows structured meshing in 2D.
'''

import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis as cfv

# ---- Define geometry ------------------------------------------------------

g = cfg.Geometry()

# Add Points:

g.point([0,0])
g.point([1.2, 0])
g.point([1, 1.3])
g.point([0, 1])
g.point([2, 0.5])

# Add Splines:
# The first four curves are structured curves, i.e the number of nodes along 
# the curves is pre-determined. Parameter el_on_curve states how many elements 
# are placed along the curve. Parameters el_distrib_type and el_distrib_val are 
# optional parameters that specify how elements are distributed.
#   "bump" means elements are bunched up at the ends or the middle of the curve.
#       In this case el_distrib_val is smaller than 1, so elements crowd at the edges.
#   "progression" means each element along the curve is larger/smaller than the previous one.
#       A larger el_distrib_val makes the elements larger at the end of the curves. 

g.spline([0,1], el_on_curve=10, el_distrib_type="bump", el_distrib_val=0.2)
g.spline([1,2], el_on_curve=20, el_distrib_type="progression", el_distrib_val=1.1)
g.spline([2,3], el_on_curve=10, el_distrib_type="bump", el_distrib_val=0.2)
g.spline([0,3], el_on_curve=20, el_distrib_type="progression", el_distrib_val=1.1) #Change order of points to reverse progression distribution
g.spline([2, 4, 1])

# Add Surfaces:
#  A structured surface must contain 4 curves that have the parameter 'el_on_curve' 
#  defined. The number of elements on two opposite curves must be the same 
#  (In this case, curves 0 & 2 and 1 & 3).

g.structuredSurface([0,1,2,3]) 
g.surface([4,1])

# ---- Create mesh ----------------------------------------------------------

mesh = cfm.GmshMesh(g)

# Element type 3 is quad. (2 is triangle. See user manual for more element types)

mesh.el_type = 3 

# Degrees of freedom per node.

mesh.dofs_per_node = 1 

# mesh.gmsh_exec_path = "D:\\vsmn20-software\\gmsh\gmsh.exe"

coords, edof, dofs, bdofs, elementmarkers = mesh.create()

# ---- Visualise mesh -------------------------------------------------------

# Draw geometry

cfv.draw_geometry(g)

# Draw mesh

cfv.figure()
cfv.draw_mesh(coords, edof, dofs_per_node=mesh.dofs_per_node, el_type=mesh.el_type, filled=True)

# Enter main loop

cfv.show_and_wait()