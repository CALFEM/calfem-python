'''Example 06

Solves a plane stress 2D problem using a structured mesh.
Shows how to draw von Mises effective stress as an element value with 
drawElementValues(). Shows use of GmshMesher attribute 'nodesOnCurve' 
(dictionary that says which nodes are on a given geometry curve)
'''

import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis as cfv
import calfem.utils as cfu
import visvis as vv

from calfem.core import *
from math import sqrt

# ---- Define problem variables ---------------------------------------------

t = 0.2
v = 0.35
E = 2.1e9
ptype = 1
ep = [ptype,t]
D=hooke(ptype, E, v)

# ---- Define geometry ------------------------------------------------------

print("Creating geometry...")

g = cfg.Geometry()

# Just a shorthand. We use this to make the circle arcs.

s2 = 1/sqrt(2) 

points = [[0, 3], [2.5, 3], [3, 3], [4-s2, 3-s2], [4, 2],     #0-4
          [4+s2, 3-s2], [5, 3], [5.5, 3], [8,3], [0, 1.5],    #5-9
          [2.5, 1.5], [4, 1.5], [5.5, 1.5], [8, 1.5], [0, 0], #10-14
          [2.5, 0], [3, 0], [4-s2, s2], [4, 1], [4+s2, s2],   #15-19
          [5, 0], [5.5, 0], [8,0], [4,3], [4,0]]              #20-24
          
for xp, yp in points:
    g.addPoint([xp*0.1, yp*0.1])

splines = [[0,1], [1,2], [6,7], [7,8], [8,13],          #0-4
           [13,22], [22,21], [21,20], [16,15], [15,14], #5-9
           [14,9], [9,0], [9,10], [10,1], [10, 15],     #10-14
           [10,11], [11,4], [11,18], [11,12], [12,7],   #15-19
           [12,21], [12,13], [3,10], [5,12], [10,17],   #20-24
           [12,19]]                                     #25
           
for s in splines:
    g.addSpline(s, elOnCurve=5)
    
g.setCurveMarker(ID=4,  marker=7) #Assign marker 7 to the splines on the right.
g.setCurveMarker(ID=5,  marker=7) # We will apply a force on nodes with marker 7.
g.setCurveMarker(ID=10, marker=5) #Assign marker 5 to the splines on the left.
g.setCurveMarker(ID=11, marker=5) # The nodes with marker 5 will be locked in place.

# Points in circle arcs are [start, center, end]

circlearcs = [[2, 23, 3], [3, 23, 4], [4, 23, 5], [5, 23, 6],           #26-29
              [16, 24, 17], [17, 24, 18], [18, 24, 19], [19, 24, 20]]   #30-33
              
for c in circlearcs:
    g.addCircle(c, elOnCurve=5)

g.addStructuredSurface([11,12,13,0]) #0
g.addStructuredSurface([14, 12, 10, 9])
g.addStructuredSurface([8, 30, 24, 14])
g.addStructuredSurface([24, 31, 17, 15])
g.addStructuredSurface([15, 16, 27, 22]) #4
g.addStructuredSurface([22, 26, 1, 13]) 
g.addStructuredSurface([16, 18, 23, 28])
g.addStructuredSurface([19, 2, 29, 23])
g.addStructuredSurface([19, 21, 4, 3]) #8
g.addStructuredSurface([20, 6, 5, 21])
g.addStructuredSurface([25, 20, 7, 33])
g.addStructuredSurface([32, 17, 18, 25]) #11

# ---- Create mesh ----------------------------------------------------------

print("Meshing geometry...")

# 3 Quads

elType = 3 
dofsPerNode = 2

# Create mesh

meshGen = cfm.GmshMeshGenerator(geometry = g)
meshGen.elType = elType
meshGen.dofsPerNode = dofsPerNode

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

print("Computing element forces...")

ed = extractEldisp(edof,a)
vonMises = []

# For each element:

for i in range(edof.shape[0]): 

    # Determine element stresses and strains in the element.
    
    es, et = planqs(ex[i,:], ey[i,:], ep, D, ed[i,:]) 

    # Calc and append effective stress to list.    
    
    vonMises.append( math.sqrt( pow(es[0],2) - es[0]*es[1] + pow(es[1],2) + 3*es[2] ) ) 

    ## es: [sigx sigy tauxy]

# ---- Visualise results ----------------------------------------------------

print("Visualising...")

cfv.drawGeometry(g, drawPoints=False, labelCurves=True)

vv.figure()
cfv.drawElementValues(vonMises, coords, edof, dofsPerNode, elType, a, doDrawMesh=True, doDrawUndisplacedMesh=False, title="Example 06 effective stress")

vv.figure()
cfv.drawDisplacements(a, coords, edof, dofsPerNode, elType, doDrawUndisplacedMesh=True, title="Example 06")

# Make use of attribute 'nodesOnCurve' in GmshMesher to draw some arrows on 
# the right hand side of the mesh:

rightSideNodes = set()

# 4 and 5 are the IDs of the curves where we applied the forces.

for curveID in [4,5]: 

    # Get the nodes, without duplicates.

    rightSideNodes = rightSideNodes.union(set(meshGen.nodesOnCurve[curveID])) 
    
for i in rightSideNodes:

    # Position of the node with displacements.

    x = coords[i,0] + a[i*2, 0]     
    y = coords[i,1] + a[i*2+1, 0]
    
    # A poor man's force indicator. Could also use vv.plot()    
    
    cfv.addText("\rightarrow", (x, y), fontSize=20, color='g') 

# Enter main loop

app = vv.use()
app.Create()
app.Run()

print("Done.")
