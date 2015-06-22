'''Example 07

Meshing 8-node-isoparametric elements (second order incomplete quads).
Shows use of surfacemarkers/elementmarkers to apply different properties to 
elements in different regions. 
'''

import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis as cfv
import calfem.utils as cfu
import visvis as vv

from calfem.core import *

# ---- Problem constants ----------------------------------------------------

kx1 = 100
ky1 = 100
kx2 = 10
ky2 = 10
t = 1.0

# Gauss points or integration points 

n = 2 
ep = [t, n]

D1 = matrix([
    [kx1, 0.],
    [0., ky1]
])
D2 = matrix([
    [kx2, 0.],
    [0., ky2]
])

# markers 10 & 11 will be used to specify different regions with different 
# conductivity.

Ddict = {10 : D1, 11 : D2} 

# ---- Create Geometry ------------------------------------------------------

g = cfg.Geometry()

# Add Points:

points = [
    [0,0], 
    [0,100], 
    [0,150], 
    [100,0], 
    [150,0], 
    [100,-100], 
    [150,-100]
]

for p in points:
    g.addPoint(p)

# Add Splines:

g.addSpline([1,2], marker=2, elOnCurve=4)
g.addSpline([3,4], elOnCurve=4)
g.addCircle([1,0,3], elOnCurve = 10)
g.addCircle([2,0,4], elOnCurve = 10)
g.addSpline([3,5], elOnCurve = 6)
g.addSpline([5,6], marker=3, elOnCurve = 4)
g.addSpline([6,4], elOnCurve = 6)

# Add Surfaces:
#
# When we set markers for surfaces, and have 2D elements, we can find which 
# region an element is in via the list 'elementmarkers', which is returned by 
# GmshMesher.create()

g.addStructuredSurface([0,2,1,3], marker = 10)
g.addStructuredSurface([1,4,5,6], marker = 11)

# Element type 16 is 8-node-quad. (See gmsh manual for more element types)

elType = 16 

#Degrees of freedom per node.

dofsPerNode = 1 

# ---- Generate mesh --------------------------------------------------------

# gmshExecPath = Path to gmsh.exe. 
# If None then the system PATH variable is queried. 
# Relative and absolute paths work.

meshGen = cfm.GmshMeshGenerator(geometry = g,
                                elType = elType,
                                dofsPerNode = dofsPerNode) 

coords, edof, dofs, bdofs, elementmarkers = meshGen.create()

# ---- Solve problem --------------------------------------------------------

print("Assembling system matrix...")

nDofs = size(dofs)
ex, ey = coordxtr(edof, coords, dofs)

K = zeros([nDofs,nDofs])

for eltopo, elx, ely, elMarker in zip(edof, ex, ey, elementmarkers):

    # Calc element stiffness matrix: Conductivity matrix D is taken 
    # from Ddict and depends on which region (which marker) the element is in.

    Ke = flw2i8e(elx, ely, ep, Ddict[elMarker]) 
    assem(eltopo, K, Ke)

print("Solving equation system...")

f = zeros([nDofs,1])

bc = array([],'i')
bcVal = array([],'i')

bc, bcVal = cfu.applybc(bdofs,bc,bcVal,2,30.0)
bc, bcVal = cfu.applybc(bdofs,bc,bcVal,3,0.0)

a,r = solveq(K,f,bc,bcVal)

# ---- Compute element forces -----------------------------------------------

print("Computing element forces...")

ed = extractEldisp(edof,a)

for i in range(shape(ex)[0]):
    es, et, eci = flw2i8s(ex[i,:], ey[i,:], ep, Ddict[elementmarkers[i]], ed[i,:])

    # Do something with es, et, eci here.
   
# ---- Visualise results ----------------------------------------------------

print("Visualising...")

cfv.drawGeometry(g, title="Geometry")

vv.figure()

# 8-node quads are drawn as simple quads.

cfv.drawMesh(coords, edof, dofsPerNode, elType, filled=False)

vv.figure()
cfv.drawNodalValues(a, coords, edof, dofsPerNode, elType, title="Example 7")
cfv.getColorbar().SetLabel("Temperature")
cfv.addText("The bend has high conductivity", (125,125))
cfv.addText("This part has low conductivity", (160,-50))

# Enter main loop

app = vv.use()
app.Create()
app.Run()

print("Done.")















