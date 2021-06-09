# import needed modules

import numpy as np
import calfem.core as cfc
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_mpl as cfv
import calfem.utils as cfu

# problem parameters
w = 100.0
h = 10.0
t = 1.0
d = h/2

# identity matrix
D = np.identity(2, 'float')
ep = [1.0, 1]


# create geometry object
g = cfg.Geometry()

g.point([0, 0])             # point 1
g.point([w, 0])             # point 2
g.point([w, h])             # point 3 
g.point([w-w/2+t/2, h])     # point 4
g.point([w-w/2+t/2, h-d])   # point 5
g.point([w-w/2-t/2, h-d])   # point 6
g.point([w-w/2-t/2, h])     # point 7
g.point([0, h])             # point 8

# create lines between points
g.spline([0, 1])
g.spline([1, 2])
g.spline([2, 3], marker=80)     # marker just to name
g.spline([3, 4])
g.spline([4, 5])
g.spline([5, 6])
g.spline([6, 7], marker=90)
g.spline([7, 0])

# make an surface area
g.surface([0, 1, 2, 3, 4, 5, 6, 7])

# plot geometry
cfv.draw_geometry(g)
cfv.showAndWait() 

# mesh generation
elType = 3          # quadrature element
dofsPerNode = 1     # 1 dof

# generate mesh paramters
meshGen = cfm.GmshMeshGenerator(g)
meshGen.elSizeFactor = 1.0              # factor that changes element sizes
meshGen.elType = elType
meshGen.dofsPerNode = dofsPerNode

# create mesh
coords, edof, dofs, bdofs, elementmarkers = meshGen.create()

# assembly
nDofs = np.size(dofs)
ex, ey = cfc.coordxtr(edof, coords, dofs)
K = np.zeros([nDofs, nDofs])

# enable loop over topology and element coordinates
for eltopo, elx, ely, in zip(edof, ex, ey):
    Ke = cfc.flw2i4e(elx, ely, ep, D)
    cfc.assem(eltopo, K, Ke)

# empty array to store loading force
f = np.zeros([nDofs, 1])    

# empty array to store boundary conditions
bc = np.array([], int)
bcVal = np.array([], int)

bc, bcVal = cfu.applybc(bdofs, bc, bcVal, 80, 0.0)
bc, bcVal = cfu.applybc(bdofs, bc, bcVal, 90, 10.0)

# solving equation
a, r = cfc.solveq(K, f, bc, bcVal)
ed = cfc.extractEldisp(edof, a)

maxFlow = []     # empty list to store flow

# calculating element force
for i in range(edof.shape[0]):
    es, et, eci = cfc.flw2i4s(ex[i, :], ey[i, :], ep, D, ed[i, :])
    maxFlow.append(np.sqrt(pow(es[0, 0], 2) + pow(es[0, 1], 2)))

# visualization
cfv.figure()
cfv.draw_geometry(g, title='Geometry')

cfv.figure()
cfv.draw_element_values(maxFlow, coords, edof, dofsPerNode, elType, None,
                        title='Max flows')

# cfv.colorBar().SetLabel("Flow")

cfv.figure()
cfv.draw_nodal_values(a, coords, edof, 
                      dofs_per_node=dofsPerNode, 
                      el_type=elType)

# cfv.colorBar().SetLabel("Node values")
cfv.showAndWait()
