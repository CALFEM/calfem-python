# -*- coding: utf-8 -*-

'''Example 07

Meshing 8-node-isoparametric elements (second order incomplete quads).
Shows use of surfacemarkers/elementmarkers to apply different properties to 
elements in different regions. 
'''

import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_mpl as cfv
import calfem.utils as cfu
import calfem.core as cfc

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.tri as tri

import numpy as np

# ---- Problem constants ----------------------------------------------------

kx1 = 100
ky1 = 100
kx2 = 10
ky2 = 10
t = 1.0

# Gauss points or integration points 

n = 2 
ep = [t, n]

D1 = np.matrix([
    [kx1, 0.],
    [0., ky1]
])
D2 = np.matrix([
    [kx2, 0.],
    [0., ky2]
])

# markers 10 & 11 will be used to specify different regions with different 
# conductivity.

Ddict = {10 : D1, 11 : D2} 

# ---- Create Geometry ------------------------------------------------------

g = cfg.geometry()

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
    g.point(p)

# Add Splines:

g.spline([1,2], marker=2, el_on_curve=4)
g.spline([3,4], el_on_curve=4)
g.circle([1,0,3], el_on_curve = 10)
g.circle([2,0,4], el_on_curve = 10)
g.spline([3,5], el_on_curve = 6)
g.spline([5,6], marker=3, el_on_curve = 4)
g.spline([6,4], el_on_curve = 6)

# Add Surfaces:
#
# When we set markers for surfaces, and have 2D elements, we can find which 
# region an element is in via the list 'elementmarkers', which is returned by 
# GmshMesher.create()

g.structuredSurface([0,2,1,3], marker = 10)
g.structuredSurface([1,4,5,6], marker = 11)

# Element type 16 is 8-node-quad. (See gmsh manual for more element types)

el_type = 16 

#Degrees of freedom per node.

dofs_per_node = 1 

# ---- Generate mesh --------------------------------------------------------

# gmshExecPath = Path to gmsh.exe. 
# If None then the system PATH variable is queried. 
# Relative and absolute paths work.

mesh = cfm.GmshMeshGenerator(g, el_type, dofs_per_node) 

coords, edof, dofs, bdofs, elementmarkers = mesh.create()

# ---- Solve problem --------------------------------------------------------

print("Assembling system matrix...")

n_dofs = np.size(dofs)
ex, ey = cfc.coordxtr(edof, coords, dofs)

K = np.zeros([n_dofs,n_dofs])

for eltopo, elx, ely, elMarker in zip(edof, ex, ey, elementmarkers):

    # Calc element stiffness matrix: Conductivity matrix D is taken 
    # from Ddict and depends on which region (which marker) the element is in.

    Ke = cfc.flw2i8e(elx, ely, ep, Ddict[elMarker]) 
    cfc.assem(eltopo, K, Ke)

print("Solving equation system...")

f = np.zeros([n_dofs,1])

bc = np.array([],'i')
bc_val = np.array([],'i')

bc, bc_val = cfu.applybc(bdofs,bc,bc_val,2,30.0)
bc, bc_val = cfu.applybc(bdofs,bc,bc_val,3,0.0)

a,r = cfc.solveq(K,f,bc,bc_val)

# ---- Compute element forces -----------------------------------------------

print("Computing element forces...")

ed = cfc.extractEldisp(edof,a)

for i in range(np.shape(ex)[0]):
    es, et, eci = cfc.flw2i8s(ex[i,:], ey[i,:], ep, Ddict[elementmarkers[i]], ed[i,:])

    # Do something with es, et, eci here.
   
# ---- Visualise results ----------------------------------------------------

print("Visualising...")

cfv.set_figure_dpi(100)

cfv.draw_geometry(g, title="Geometry")

# 8-node quads are drawn as simple quads.

cfv.figure()
cfv.draw_mesh(coords, edof, dofs_per_node, el_type, filled=False)

cfv.figure()
cfv.draw_nodal_values_shaded(a, coords, edof, title="Temperature", dofs_per_node=mesh.dofs_per_node, el_type=mesh.el_type, draw_elements=True)
cbar = cfv.colorbar(orientation="horizontal")
cbar.set_label("Temperature")


cfv.text("The bend has high conductivity", (125,125))
cfv.text("This part has low conductivity", (160,-50))

# Enter main loop

cfv.show_and_wait()

print("Done.")















