'''Example 05

This example shows how to make an unstructured 3D mesh (tetrahedron elements, which calfem cant actually use).
It also demonstrates how to do subplots and create two axes that are viewed from the same camera.
''' 

import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis as cfv
import visvis as vv

# ---- Define geometry ------------------------------------------------------

g = cfg.Geometry()

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

# ---- Create mesh ----------------------------------------------------------

# Element type 4 is tetrahedron. (See user manual for more element types).

elType = 4 

# Degrees of freedom per node.

dofsPerNode= 1 

# Create mesh

coords, edof, dofs, bdofs, elementmarkers = cfm.createGmshMesh(geometry=g,
                                                               elType = elType,
                                                               elSizeFactor = 0.3,
                                                               dofsPerNode = dofsPerNode)

# ---- Visualise mesh -------------------------------------------------------

# Create two axes that are viewed from the same camera:

vv.figure()
a1 = vv.subplot(121)
a2 = vv.subplot(122)
cam = vv.cameras.ThreeDCamera()
a1.camera = a2.camera = cam

# Draw geometry and mesh

cfv.drawGeometry(g, axes=a1)
cfv.drawMesh(coords=coords, edof=edof, dofsPerNode=dofsPerNode, elType=elType, filled=False, axes=a2)

# Enter main loop

app = vv.use()
app.Create()
app.Run()