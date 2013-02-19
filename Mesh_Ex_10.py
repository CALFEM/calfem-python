'''Example 10
The use case from the user manual. 
The example does not contain anything that is not covered in the previous examples.
'''

import pycalfem_GeoData
import pycalfem_mesh
import pycalfem_vis as pcv
import visvis as vv
from pycalfem import *
from pycalfem_utils import *

#DEFINE GEOMETRY:
g = pycalfem_GeoData.GeoData() #Create a GeoData object that holds the geometry.

#Add points:
g.addPoint([0, 0])		#0
g.addPoint([1, 0])		#1
g.addPoint([1, 1])		#2
g.addPoint([0, 1])		#3
g.addPoint([0.2, 0.2])	#4
g.addPoint([0.8, 0.2])	#5
g.addPoint([0.8, 0.8])	#6
g.addPoint([0.2, 0.8])	#7

#Add curves:
g.addSpline([0, 1], marker=70)        #0
g.addSpline([2, 1])                   #1
g.addSpline([3, 2], marker=90)        #2
g.addSpline([0, 3])                   #3
g.addSpline([4, 5])                 #4
g.addSpline([5, 6])                 #5
g.addSpline([6, 7])                 #6
g.addSpline([7, 4])                 #7

#Add surfaces:
g.addSurface([0,1,2,3], holes=[[4,5,6,7]], marker = 55)
g.addSurface([4,5,6,7], marker = 66)

#MESHING:
elType = 3 #Element type 3 is quad.
dofsPerNode = 2 #Degrees of freedom per node.

mesher = pycalfem_mesh.GmshMesher(geoData = g,
							gmshExecPath = None, #Path to gmsh.exe. If None then the system PATH variable is queried. Relative and absolute paths work.
							elSizeFactor = 0.04, #Factor that changes element sizes.
							elType = elType, 
							dofsPerNode= dofsPerNode)

#Mesh the geometry:
# The first four return values are the same as those that trimesh2d() returns.
# value elementmarkers is a list of markers, and is used for finding the marker of a given element (index).
coords, edof, dofs, bdofs, elementmarkers = mesher.create()


#SOLVE:
t = 0.2
v = 0.35
E1 = 2e9
E2 = 0.2e9
ptype = 1
ep = [ptype,t]
D1 = hooke(ptype, E1, v)
D2 = hooke(ptype, E2, v)
Ddict = {55 : D1, 66 : D2} 

nDofs = size(dofs)
K = zeros([nDofs,nDofs])
ex, ey = coordxtr(edof, coords, dofs)

for eltopo, elx, ely, elMarker in zip(edof, ex, ey, elementmarkers):
    #Calc element stiffness matrix: Conductivity matrix D is taken 
    # from Ddict and depends on which region (which marker) the element is in.
    Ke = planqe(elx, ely, ep, Ddict[elMarker])
    assem(eltopo, K, Ke)

bc = array([],'i')
bcVal = array([],'i')
bc, bcVal = applybc(bdofs,bc,bcVal, 70, 0.0)

f = zeros([nDofs,1])
applyforce(bdofs, f, 90, value = -10e5, dimension=2)
a,r = solveq(K,f,bc,bcVal)

ed = extractEldisp(edof,a)
vonMises = []
for i in range(edof.shape[0]): #For each element:
    es, et = planqs(ex[i,:], ey[i,:], ep, Ddict[elementmarkers[i]], ed[i,:]) #Determine element stresses and strains in the element.
    vonMises.append( math.sqrt( pow(es[0],2) - es[0]*es[1] + pow(es[1],2) + 3*pow(es[2],2) ) )


#VISUALISATION:
pcv.drawGeometry(g, title="Geometry")
vv.figure() #New figure window
pcv.drawMesh(coords=coords, edof=edof, dofsPerNode=dofsPerNode, elType=elType, filled=True, title="Mesh") #Draws the mesh.
vv.figure()
pcv.drawDisplacements(a, coords, edof, dofsPerNode, elType, doDrawUndisplacedMesh=False, title="Displacements")
vv.figure()
pcv.drawElementValues(vonMises, coords, edof, dofsPerNode, elType, a, doDrawMesh=True, doDrawUndisplacedMesh=False, title="Effective Stress")
pcv.getColorbar().SetLabel("Effective stress")

# Enter main loop:
app = vv.use()
#app.Create()
app.Run()