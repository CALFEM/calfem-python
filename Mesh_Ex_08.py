'''Example 08

Shows how to load geometry from a gmsh geo file.
The geometry is the same as in example 06.
'''

import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis as cfv
import calfem.utils as cfu
import visvis as vv

from calfem.core import *

# ---- Problem variables ----------------------------------------------------

t = 0.2
v = 0.35
E = 2.1e9
ptype = 1
ep = [ptype,t]
D = hooke(ptype, E, v)

# ---- Create Geometry ------------------------------------------------------

# Relative path to the geo file that contains the geometry.

g = "examplegeo\ex8.geo" 

elType = 3 #3 Quads
dofsPerNode = 2

# ---- Generate mesh --------------------------------------------------------

meshGen = cfm.GmshMeshGenerator(geometry = g,
                                gmshExecPath = None,
                                elType = elType, 
                                dofsPerNode = dofsPerNode)
                    
coords, edof, dofs, bdofs, elementmarkers = meshGen.create()

# ---- Solve problem --------------------------------------------------------

print("Assembling system matrix...")

nDofs = size(dofs)
ex, ey = coordxtr(edof, coords, dofs)
K = zeros([nDofs,nDofs])

for eltopo, elx, ely in zip(edof, ex, ey):
    Ke = planqe(elx, ely, ep, D)
    assem(eltopo, K, Ke)

print("Solving equation system...")

f = zeros([nDofs,1])
bc = array([],'i')
bcVal = array([],'i')
bc, bcVal = cfu.applybc(bdofs, bc, bcVal, 5, 0.0, 0)
cfu.applyforce(bdofs, f, 7, 10e5, 1)

a,r = solveq(K,f,bc,bcVal)

# ---- Visualise results ----------------------------------------------------

print("Visualising...")

cfv.drawDisplacements(a, coords, edof, dofsPerNode, elType, doDrawUndisplacedMesh=True, title="Example 08 - Geometry from 06")

# Enter main loop

app = vv.use()
app.Create()
app.Run()

print("Done.")
