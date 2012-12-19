'''Example 01
Shows how to create simple geometry from splines and ellipse arcs, and how to mesh a quad mesh in GmshMesher.
Also demonstrates drawGeometry(), drawMesh, and drawing texts and labels in a figure.
'''

from pycalfem_GeoData import *
from pycalfem_mesh import *
import pycalfem_vis as pcv
import visvis as vv

#DEFINE GEOMETRY:
g = GeoData() #Create a GeoData object that holds the geometry.

#Add points:
# The first parameter is the coordinates. These can be in 2D or 3D.
# The other parameters are not defined in this example. These parameters are
# ID, marker, and elSize.
# Since we do not specify an ID the points are automatically assigned IDs, starting from 0.
g.addPoint([0, 0])
g.addPoint([2, 0])
g.addPoint([2, 1])
g.addPoint([0, 1])
g.addPoint([0.5, 0.3])
g.addPoint([0.3, 0.7])
g.addPoint([0.7, 0.7])
g.addPoint([0.8, 0.5])
g.addPoint([1.7, 0.5])
g.addPoint([1.5, 0.5])
g.addPoint([1.7, 0.7])

#Add curves:
# There are four types of curves. In this example we create an ellipse arc and some splines.
# The first parameter is a list of point IDs that make the curve.
# Curves can have have IDs and markers. In this example the IDs are undefined so the curves are
# automatically assigned IDs. The markers can be used for identifying regions/boundaries in the
# model.
g.addEllipse([7,8,9,10], marker=50)   #0 An ellipse arc. Read the function doc for more information.
g.addSpline([0, 1], marker=80)        #1
g.addSpline([2, 1])                   #2
g.addSpline([3, 2])                   #3
g.addSpline([0, 3])                   #4
g.addSpline([7, 9], marker=50)        #5
g.addSpline([10, 9])                  #6
g.addSpline([4, 5, 6, 4])             #7 This is a closed spline. The start and end points are the same (4).

#Add a surface:
# Surfaces are defined by its curve boundaries.
# The first parameter is a list of curve IDs that specify the outer boundary of the surface.
# The second parameter is a list of lists of curve IDs that specify holes in the surface.
# In this example there are two holes.
# The boundaries and holes must be closed paths. We can see that [7] is closed because curve 7
# is a closed spline.
# addSurface creates a flat surface, so all curves must lie on the same plane.
g.addSurface([4,3,2,1], [[7], [5,6,0]])

#MESHING:
elmType = 3 #Element type 3 is quad. (2 is triangle. See gmsh manual for more element types)
dofsPerNode= 1 #Degrees of freedom per node.

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
pcv.drawGeometry(g)#Draws the geometry. Note that surfaces and volumes are not drawn at all by this function.

vv.figure() #New figure window

pcv.drawMesh(coords=coords, edof=edof, dofsPerNode=dofsPerNode, elmType=elmType, filled=True, title="Example 01") #Draws the mesh.

pcv.addText("This is a Text", pos=(1, -0.3), angle=45)  #Adds a text in world space

ourLabel = pcv.addLabel("This is a Label", pos=(100,200), angle=-45) #Adds a label in the screen space
ourLabel.text = "Label, changed." #We can change the attributes of labels and texts, such as color, text, and position.
ourLabel.textColor = 'r'  #Make it red. (1,0,0) would also have worked.
ourLabel.position = (20,30)

# Enter main loop:
app = vv.use()
app.Create()
app.Run()