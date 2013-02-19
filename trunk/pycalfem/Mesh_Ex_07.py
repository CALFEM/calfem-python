'''Example 07
Meshing 8-node-isoparametric elements (second order incomplete quads).
Shows use of surfacemarkers/elementmarkers to apply different properties to elements in different regions. 
'''
from pycalfem_GeoData import *
from pycalfem_mesh import *
from pycalfem_vis import *
import visvis as vv
from pycalfem import *
from pycalfem_utils import *

# ---- Problem constants
kx1 = 100
ky1 = 100
kx2 = 10
ky2 = 10
t = 1.0
n = 2 #Gauss points or integration points
ep = [t, n]

D1 = matrix([
    [kx1, 0.],
    [0., ky1]])
D2 = matrix([
    [kx2, 0.],
    [0., ky2]])
Ddict = {10 : D1, 11 : D2} #markers 10 & 11 will be used to specify different regions with different conductivity.

# ---- Create Geometry
g = GeoData()

#Add Points:
points = [[0,0], [0,100], [0,150], [100,0], [150,0], [100,-100], [150,-100]]
for p in points:
    g.addPoint(p)

#Add Splines:
g.addSpline([1,2], marker=2, elOnCurve=4)
g.addSpline([3,4], elOnCurve=4)
g.addCircle([1,0,3], elOnCurve = 10)
g.addCircle([2,0,4], elOnCurve = 10)
g.addSpline([3,5], elOnCurve = 6)
g.addSpline([5,6], marker=3, elOnCurve = 4)
g.addSpline([6,4], elOnCurve = 6)

#Add Surfaces:
# When we set markers for surfaces, and have 2D elements, we can find which region
# an element is in via the list 'elementmarkers', which is returned by GmshMesher.create()
g.addStructuredSurface([0,2,1,3], marker = 10)
g.addStructuredSurface([1,4,5,6], marker = 11)

elType = 16 #Element type 16 is 8-node-quad. (See gmsh manual for more element types)
dofsPerNode= 1 #Degrees of freedom per node.

mesher = GmshMesher(geoData = g,
                    gmshExecPath = None, #Path to gmsh.exe. If None then the system PATH variable is queried. Relative and absolute paths work.
                    elType = elType,
                    dofsPerNode= dofsPerNode)
coords, edof, dofs, bdofs, elementmarkers = mesher.create()

print "Assembling system matrix..."
nDofs = size(dofs)
ex, ey = coordxtr(edof, coords, dofs)

K = zeros([nDofs,nDofs])

for eltopo, elx, ely, elMarker in zip(edof, ex, ey, elementmarkers):
    #Calc element stiffness matrix: Conductivity matrix D is taken 
    # from Ddict and depends on which region (which marker) the element is in.
    Ke = flw2i8e(elx, ely, ep, Ddict[elMarker]) 
    assem(eltopo, K, Ke)


print "Solving equation system..."
f = zeros([nDofs,1])

bc = array([],'i')
bcVal = array([],'i')

bc, bcVal = applybc(bdofs,bc,bcVal,2,30.0)
bc, bcVal = applybc(bdofs,bc,bcVal,3,0.0)

a,r = solveq(K,f,bc,bcVal)


print "Computing element forces..."
ed = extractEldisp(edof,a)

for i in range(shape(ex)[0]):
    es, et, eci = flw2i8s(ex[i,:], ey[i,:], ep, Ddict[elementmarkers[i]], ed[i,:])
    #Do something with es, et, eci here.
   
print "Visualising..."
drawGeometry(g, title="Geometry")

vv.figure()
drawMesh(coords, edof, dofsPerNode, elType, filled=False)
#8-node quads are drawn as simple quads.

vv.figure()
drawNodalValues(a, coords, edof, dofsPerNode, elType, title="Example 7")
getColorbar().SetLabel("Temperature")
addText("The bend has high conductivity", (125,125))
addText("This part has low conductivity", (160,-50))

# Enter main loop:
app = vv.use()
app.Create()
app.Run()

print "Done."















