'''Example 03
Shows structured meshing in 2D.
'''

from pycalfem_GeoData import *
from pycalfem_mesh import *
import pycalfem_vis as pcv
import visvis as vv

#DEFINE GEOMETRY:
g = GeoData()

#Add Points:
g.addPoint([0,0])
g.addPoint([1.2, 0])
g.addPoint([1, 1.3])
g.addPoint([0, 1])
g.addPoint([2, 0.5])

#Add Splines:
#The first four curves are structured curves, i.e the number of nodes along the curves is pre-determined.
#Parameter elOnCurve states how many elements are placed along the curve.
#Parameters elDistribType and elDistribVal are optional parameters that specify how elements are distributed.
#   "bump" means elements are bunched up at the ends or the middle of the curve.
#       In this case elDistribVal is smaller than 1, so elements crowd at the edges.
#   "progression" means each element along the curve is larger/smaller than the previous one.
#       A larger elDistribVal makes the elements larger at the end of the curves. 
g.addSpline([0,1], elOnCurve=10, elDistribType="bump", elDistribVal=0.2)
g.addSpline([1,2], elOnCurve=20, elDistribType="progression", elDistribVal=1.1)
g.addSpline([2,3], elOnCurve=10, elDistribType="bump", elDistribVal=0.2)
g.addSpline([0,3], elOnCurve=20, elDistribType="progression", elDistribVal=1.1) #Change order of points to reverse progression distribution
g.addSpline([2, 4, 1])

#Add Surfaces:
#A structured surface must contain 4 curves that have the parameter 'elOnCurve' defined. 
#The number of elements on two opposite curves must be the same (In this case, curves 0 & 2 and 1 & 3).
g.addStructuredSurface([0,1,2,3]) 
g.addSurface([4,1])

#MESHING:
elType = 3 #Element type 3 is quad. (2 is triangle. See user manual for more element types)
dofsPerNode= 1 #Degrees of freedom per node.

mesher = GmshMesher(geoData = g,
                    gmshExecPath = None, #Path to gmsh.exe. If None then the system PATH variable is queried. Both relative and absolute paths work.
                    elSizeFactor = 0.05, #Factor that changes element sizes.
                    elType = elType, 
                    dofsPerNode= dofsPerNode)

#Mesh the geometry:
coords, edof, dofs, bdofs, elementmarkers = mesher.create()

#VISUALISATION:
pcv.drawGeometry(g)#Draws the geometry.
vv.figure() #New figure window
pcv.drawMesh(coords=coords, edof=edof, dofsPerNode=dofsPerNode, elType=elType, filled=True) #Draws the mesh.
# Enter main loop:
app = vv.use()
app.Create()
app.Run()
