'''Example 02
Creating geometry from B-Splines and circle arcs.
Also shows how to set ID numbers for geometry entities and how to specify element density. 
'''

from pycalfem_GeoData import *
from pycalfem_mesh import *
import pycalfem_vis as pcv
import visvis as vv

#DEFINE GEOMETRY:
g = GeoData()

#Add points:
#In this example we set the IDs manually.
g.addPoint([ -2,  0], ID=0)
g.addPoint([  0,  1], ID=1, elSize=5) #elSize determines the size of the elements near this point.
g.addPoint([  1,  0], 2, elSize=5)    # elSize is 1 by default. Larger number means less dense mesh.
g.addPoint([  0, -2], 3)
g.addPoint([  0,  0], 4, elSize=5)
g.addPoint([ .5, .2], 5)
g.addPoint([-.5, .5], 6)
g.addPoint([-.7,-.5], 7)
#Add curves:
g.addCircle([1, 4, 2], 2) #The 3 points that define the circle arc are [start, center, end]. The arc must be smaller than Pi.
g.addBSpline([5,6,7,5], 5) # BSplines are similar to Splines, but do not pass through all control points.
g.addBSpline([1,0,3,2], 4)
#Add surface:
g.addSurface([4,2], [[5]])

#Markers do not have to be set when the curve is created. It can be done afterwards.
# Set marker=80 for curves 2 and 4:
for curveID in [2, 4]:
    g.setCurveMarker(curveID, 80)

#MESHING:
elmType = 2 #Element type 2 is triangle. (3 is quad. See gmsh manual for more element types)
dofsPerNode= 2 #Degrees of freedom per node.

mesher = GmshMesher(geoData = g,
                    gmshExecPath = None, #Path to gmsh.exe. If None then the system PATH variable is queried. Relative and absolute paths work.
                    elmSizeFactor = 0.05, #Factor that changes element sizes.
                    elmType = elmType, 
                    dofsPerNode= dofsPerNode)

#Mesh the geometry:
# The first four return values are the same as those that trimesh2d() returns.
# value elementmarkers is a list of markers, and is used for finding the marker of a given element (index).
coords, edof, dofs, bdofs, elementmarkers = mesher.create()

#VISUALISATION:
#Hold left mouse button to pan.
#Hold right mouse button to zoom.
pcv.drawGeometry(g, labelCurves=True)#Draws the geometry.

vv.figure() #New figure window

pcv.drawMesh(coords=coords, edof=edof, dofsPerNode=dofsPerNode, elmType=elmType, filled=True, title="Example 02") #Draws the mesh.

# Enter main loop:
app = vv.use()
app.Create()
app.Run()