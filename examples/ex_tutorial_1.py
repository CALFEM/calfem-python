# -*- coding: utf-8 -*-
"""
Created on Sat Mar  3 22:08:29 2018

@author: Jonas Lindemann
"""

import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_mpl as cfv

# ----- Define geometry

g = cfg.Geometry()

g.point([0.0, 0.0]) # point 0
g.point([5.0, 0.0], marker=20) # point 1
g.point([2.5, 4.0]) # point 2

g.spline([0, 1]) # line 0
g.spline([1, 2]) # line 1
g.spline([2, 0], marker=10) # line 2

g.surface([0, 1, 2])

# ----- Create mesh

mesh = cfm.GmshMesh(g)

mesh.elType = 2 # Degrees of freedom per node.
mesh.dofsPerNode = 1 # Factor that changes element sizes.
mesh.elSizeFactor = 0.15

coords, edof, dofs, bdofs, elementmarkers = mesh.create()

print(bdofs)

cfv.draw_geometry(g)

cfv.figure() 

# ----- Draw the mesh.

cfv.draw_mesh(
    coords=coords, 
    edof=edof, 
    dofs_per_node=mesh.dofsPerNode, 
    el_type=mesh.elType, 
    filled=True, 
    title="Example 01"
    ) 

cfv.showAndWait()