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

class RectSolver(cfslv.Solver):
        
    def onInit(self):
        self.vonMises = []

    def onCreateKe(self, elx, ely, elementType):
        Ke = None
        if self.mesh.shape.elementType == 2:
            Ke = cfc.plante(elx, ely, self.mesh.shape.ep, self.mesh.shape.D)
        else:
            Ke = cfc.planqe(elx, ely, self.mesh.shape.ep, self.mesh.shape.D)
            
        return Ke
    
    def onApplyBCs(self, mesh, bc, bcVal):        
        bc, bcVal = cfu.applybc(mesh.bdofs, bc, bcVal, self.mesh.shape.leftId, 0.0)
        
    def onApplyLoads(self, mesh, f):
        cfu.applyforcetotal(mesh.bdofs, f, self.mesh.shape.topId, value = -10e5, dimension=2)
        
    def onCalcElForce(self, ex, ey, ed):
        if self.rect.elementType == 2: 
            es, et = cfc.plants(ex, ey, self.mesh.shape.ep, self.mesh.shape.D, ed)
            self.vonMises.append( np.math.sqrt( pow(es[0,0],2) - es[0,0]*es[0,1] + pow(es[0,1],2) + 3*pow(es[0,2],2) ) )
        else:
            es, et = cfc.planqs(ex, ey, self.rect.ep, self.rect.D, ed)
            self.vonMises.append( np.math.sqrt( pow(es[0],2) - es[0]*es[1] + pow(es[1],2) + 3*pow(es[2],2) ) )
        


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

solver = RectSolver(mesh)
       
# ---- Visualise results ----------------------------------------------------

print("Drawing results...")

cfv.figure() 
cfv.drawGeometry(rect.geometry(), title="Geometry")

cfv.figure() 
cfv.drawMesh(mesh.coords, mesh.edof, rect.dofsPerNode, rect.elementType, 
             filled=True, title="Mesh") #Draws the mesh.

cfv.figure()
cfv.drawDisplacements(solver.a, mesh.coords, mesh.edof, rect.dofsPerNode, rect.elementType, 
                      doDrawUndisplacedMesh=False, title="Displacements", 
                      magnfac=1)

#cfv.figure()
#cfv.drawElementValues(solver.vonMises, mesh.coords, mesh.edof, rect.dofsPerNode, rect.elementType, solver.a, 
#                      doDrawMesh=True, doDrawUndisplacedMesh=False, 
#                      title="Effective Stress", magnfac=1)
                      
#cfv.colorBar().SetLabel("Effective stress")

print("Done drawing...")

cfv.showAndWait()