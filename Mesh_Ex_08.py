'''Example 08
Shows how to load geometry from a gmsh geo file.
The geometry is the same as in example 06.
'''

from pycalfem import *
from pycalfem_utils import *
from pycalfem_mesh import *
import pycalfem_vis as pcv
import visvis as vv

# ---- Create Geometry
g = "examplegeo\ex8.geo" #Relative path to the geo file that contains the geometry.

elmType = 3 #3 Quads
dofsPerNode = 2

mesher = GmshMesher(geoData = g,
                    gmshExecPath = None, #"gmsh\gmsh.exe"
                    elmType = elmType, 
                    dofsPerNode= dofsPerNode)
coords, edof, dofs, bdofs, elementmarkers = mesher.create()

#The rest of the example is the same as example 06.
t = 0.2
v = 0.35
E = 2.1e9
ptype = 1
ep = [ptype,t]
D=hooke(ptype, E, v)

print "Assembling system matrix..."
nDofs = size(dofs)
ex, ey = coordxtr(edof, coords, dofs)
K = zeros([nDofs,nDofs])

for eltopo, elx, ely in zip(edof, ex, ey):
    Ke = planqe(elx, ely, ep, D)
    assem(eltopo, K, Ke)

print "Solving equation system..."
f = zeros([nDofs,1])
bc = array([],'i')
bcVal = array([],'i')
bc, bcVal = applybc(bdofs, bc, bcVal, 5, 0.0, 0)
applyforce(bdofs, f, 7, 10e5, 1)
a,r = solveq(K,f,bc,bcVal)

print "Drawing element mesh..."
pcv.drawDisplacements(a, coords, edof, dofsPerNode, elmType, doDrawUndisplacedMesh=True, title="Example 08 - Geometry from 06")

# Enter main loop:
app = vv.use()
app.Create()
app.Run()

print "Done."
