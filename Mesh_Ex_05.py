'''Example 05
This example shows how to make an unstructured 3D mesh (tetrahedron elements, which pycalfem cant actually use).
It also demonstrates how to do subplot and create two axes that are viewed from the same camera.
''' 

from pycalfem_GeoData import *
from pycalfem_mesh import *
import pycalfem_vis as pcv
import visvis as vv

g = GeoData()

g.addPoint([0, 0, 0],      0)
g.addPoint([1, 0, 0],      1)
g.addPoint([0, 1, 0],      2)
g.addPoint([0, 1, 1],      3, elSize=0.1)
g.addPoint([0.5, -0.3, 0], 4)
g.addPoint([-0.3, 0.5, 0], 5)
g.addPoint([0.75, 0.75, 0],6)

g.addSpline([0,4,1])
g.addSpline([1,6,2])
g.addSpline([2,5,0])
g.addSpline([0,3])
g.addSpline([3,2])
g.addSpline([3,1])

g.addRuledSurface([0,1,2])
g.addRuledSurface([0,5,3])
g.addRuledSurface([1,5,4])
g.addRuledSurface([2,3,4])

g.addVolume([0,1,2,3])

elmType = 4 #Element type 4 is tetrahedron. Pycalfem (See gmsh manual for more element types).
dofsPerNode= 1 #Degrees of freedom per node.

mesher = GmshMesher(geoData = g,
                    gmshExecPath = None,
                    elmType = elmType,
                    elmSizeFactor = 0.3,
                    dofsPerNode= dofsPerNode)

#Mesh the geometry:
coords, edof, dofs, bdofs, elementmarkers = mesher.create()

#VISUALISATION:

#Create two axes that are viewed from the same camera:
vv.figure()
a1 = vv.subplot(121)
a2 = vv.subplot(122)
cam = vv.cameras.ThreeDCamera()
a1.camera = a2.camera = cam
#Draw:
pcv.drawGeometry(g, axes=a1)#Draws the geometry.
pcv.drawMesh(coords=coords, edof=edof, dofsPerNode=dofsPerNode, elmType=elmType, filled=False, axes=a2) #Draws the mesh.

# Enter main loop:
app = vv.use()
app.Create()
app.Run()