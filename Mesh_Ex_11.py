'''Example 10

The use case from the user manual. 
The example does not contain anything that is not covered in the previous examples.
'''

import calfem.core as cfc
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis as cfv
import calfem.utils as cfu
import calfem.shapes as cfs

import numpy as np

from scipy.sparse import lil_matrix

# ---- General parameters ---------------------------------------------------

cfu.enableLogging()

t = 0.2
v = 0.35
E = 2e9
ptype = 1
ep = [ptype,t]
D = cfc.hooke(ptype, E, v)

# Define marker constants instead of using numbers in the code

print("Creating rectangle")
rect = cfs.Rectangle(5.0, 1.0, elementType=3, dofsPerNode=2, maxArea=0.08)

print("Creating mesh...")

mesh = cfs.ShapeMesh(rect)    

# ---- Solve problem --------------------------------------------------------

nDofs = np.size(mesh.dofs)
K = lil_matrix((nDofs,nDofs))
#ex, ey = cfc.coordxtr(mesh.edof, mesh.coords, mesh.dofs)

print("Assembling K... ("+str(nDofs)+")")

for eltopo, elx, ely in zip(mesh.edof, mesh.ex, mesh.ey):

    if rect.elementType == 2:
        Ke = cfc.plante(elx, ely, ep, D)
    else:
        Ke = cfc.planqe(elx, ely, ep, D)
        
        
    cfc.assem(eltopo, K, Ke)
    
print("Applying bc and loads...")

bc = np.array([],'i')
bcVal = np.array([],'i')

bc, bcVal = cfu.applybc(mesh.bdofs, bc, bcVal, rect.leftId, 0.0)

f = np.zeros([nDofs,1])

cfu.applyforcetotal(mesh.bdofs, f, rect.topId, value = -10e5, dimension=2)

print("Solving system...")

a,r = cfc.spsolveq(K, f, bc, bcVal)

print("Extracting ed...")

ed = cfc.extractEldisp(mesh.edof, a)
vonMises = []

# ---- Calculate elementr stresses and strains ------------------------------

print("Element forces... ")

for i in range(mesh.edof.shape[0]):
    
    # Handle triangle elements
        
    if rect.elementType == 2: 
        es, et = cfc.plants(mesh.ex[i,:], mesh.ey[i,:], 
                        ep, 
                        D, 
                        ed[i,:])
        
        vonMises.append( np.math.sqrt( pow(es[0,0],2) - es[0,0]*es[0,1] + pow(es[0,1],2) + 3*pow(es[0,2],2) ) )

    else:
        
        # Handle quad elements
        
        es, et = cfc.planqs(mesh.ex[i,:], mesh.ey[i,:], 
                        ep, 
                        D, 
                        ed[i,:])
        
        vonMises.append( np.math.sqrt( pow(es[0],2) - es[0]*es[1] + pow(es[1],2) + 3*pow(es[2],2) ) )
        
# ---- Visualise results ----------------------------------------------------

print("Drawing results...")

cfv.figure() 
cfv.drawGeometry(rect.geometry(), title="Geometry")

cfv.figure() 
cfv.drawMesh(mesh.coords, mesh.edof, rect.dofsPerNode, rect.elementType, 
             filled=True, title="Mesh") #Draws the mesh.

cfv.figure()
cfv.drawDisplacements(a, mesh.coords, mesh.edof, rect.dofsPerNode, rect.elementType, 
                      doDrawUndisplacedMesh=False, title="Displacements", 
                      magnfac=1)

cfv.figure()
cfv.drawElementValues(vonMises, mesh.coords, mesh.edof, rect.dofsPerNode, rect.elementType, a, 
                      doDrawMesh=True, doDrawUndisplacedMesh=False, 
                      title="Effective Stress", magnfac=1)
                      
cfv.colorBar().SetLabel("Effective stress")

print("Done drawing...")

cfv.showAndWait()