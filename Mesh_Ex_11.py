'''Example 10

The use case from the user manual. 
The example does not contain anything that is not covered in the previous examples.
'''

import calfem.core as cfc
import calfem.vis as cfv
import calfem.utils as cfu
import calfem.shapes as cfs
import calfem.solver as cfslv
  
# ---- General parameters ---------------------------------------------------

cfu.enableLogging()

# Define marker constants instead of using numbers in the code

print("Creating rectangle")

rect = cfs.Rectangle(5.0, 1.0, elementType=3, dofsPerNode=2, maxArea=0.08)
rect.t = 0.2
rect.v = 0.35
rect.E = 2e9
rect.ptype = 1
rect.ep = [rect.ptype, rect.t]
rect.D = cfc.hooke(rect.ptype, rect.E, rect.v)

print("Creating mesh...")

mesh = cfs.ShapeMesh(rect)    

# ---- Solve problem --------------------------------------------------------

solver = cfslv.Plan2DSolver(mesh)

solver.addBC(rect.leftId, 0.0)
solver.addForceTotal(rect.topId, -10e5, dimension=2)

results = solver.execute()       
       
# ---- Visualise results ----------------------------------------------------

print("Drawing results...")

cfv.figure() 
cfv.drawGeometry(rect.geometry(), title="Geometry")

cfv.figure() 
cfv.drawMesh(mesh.coords, mesh.edof, rect.dofsPerNode, rect.elementType, 
             filled=True, title="Mesh") #Draws the mesh.

cfv.figure()
cfv.drawDisplacements(results.a, mesh.coords, mesh.edof, rect.dofsPerNode, rect.elementType, 
                      doDrawUndisplacedMesh=False, title="Displacements", 
                      magnfac=1)

cfv.figure()
cfv.drawElementValues(results.elForces, mesh.coords, mesh.edof, rect.dofsPerNode, rect.elementType, results.a, 
                      doDrawMesh=True, doDrawUndisplacedMesh=False, 
                      title="Effective Stress", magnfac=1)
                      
#cfv.colorBar().SetLabel("Effective stress")

print("Done drawing...")

cfv.showAndWait()