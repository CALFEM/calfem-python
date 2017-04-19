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

rect = cfs.Rectangle(5.0, 1.0, elementType=3, dofsPerNode=1, maxArea=0.08)

rect.t = 1
rect.ep = [rect.t, 1]

rect.D = np.diag([1.7, 1.7])

cfu.info("Creating mesh...")

mesh = cfs.ShapeMesh(rect)    

# ---- Solve problem --------------------------------------------------------

solver = cfslv.Flow2DSolver(mesh)

solver.addBC(rect.leftId, 0.0)
solver.addBC(rect.rightId, 120.0)
#solver.addForceTotal(rect.topId, -10e5, dimension=2)

results = solver.execute()       
       
# ---- Visualise results ----------------------------------------------------

cfu.info("Drawing results...")

cfv.figure() 
cfv.drawGeometry(rect.geometry(), title="Geometry")

cfv.figure() 
cfv.drawMesh(mesh.coords, mesh.edof, rect.dofsPerNode, rect.elementType, 
             filled=True, title="Mesh") #Draws the mesh.

cfv.figure() 
cfv.drawNodalValues(results.a, mesh.coords, mesh.edof, rect.dofsPerNode, rect.elementType)

print(results.a)

#cfv.figure()
#cfv.drawDisplacements(results.a, mesh.coords, mesh.edof, rect.dofsPerNode, rect.elementType, 
#                      doDrawUndisplacedMesh=False, title="Displacements", 
#                      magnfac=1)

#cfv.figure()
#cfv.drawElementValues(results.elForces, mesh.coords, mesh.edof, rect.dofsPerNode, rect.elementType, results.a, 
#                      doDrawMesh=True, doDrawUndisplacedMesh=False, 
#                      title="Effective Stress", magnfac=1)
                      
#cfv.colorBar().SetLabel("Effective stress")

cfu.info("Done drawing...")

cfv.showAndWait()