# -*- coding: utf-8 -*-

'''Example 10

The use case from the user manual. 
The example does not contain anything that is not covered in the previous examples.
'''

import calfem.core as cfc
import calfem.vis as cfv
import calfem.utils as cfu
import calfem.shapes as cfs
import calfem.solver as cfslv

import numpy as np

cfu.enableLogging()
  
# ---- General parameters ---------------------------------------------------

# Define marker constants instead of using numbers in the code

cfu.info("Creating rectangle")

rect = cfs.Rectangle(5.0, 1.0, element_type=3, dofs_per_node=1, max_area=0.08)

rect.t = 1
rect.ep = [rect.t, 1]

rect.D = np.diag([1.7, 1.7])

cfu.info("Creating mesh...")

mesh = cfs.ShapeMesh(rect)    

# ---- Solve problem --------------------------------------------------------

solver = cfslv.Flow2DSolver(mesh)

solver.addBC(rect.left_id, 0.0)
solver.addBC(rect.right_id, 120.0)
#solver.addForceTotal(rect.topId, -10e5, dimension=2)

results = solver.execute()       
       
# ---- Visualise results ----------------------------------------------------

cfu.info("Drawing results...")

cfv.figure() 
cfv.drawGeometry(rect.geometry(), title="Geometry")

cfv.figure() 
cfv.drawMesh(mesh.coords, mesh.edof, rect.dofs_per_node, rect.element_type, 
             filled=True, title="Mesh") #Draws the mesh.

cfv.figure() 
cfv.drawNodalValues(results.a, mesh.coords, mesh.edof, rect.dofs_per_node, rect.element_type)

cfv.showAndWait()