# -*- coding: utf-8 -*-
"""
Created on Sat Mar  3 22:08:29 2018

@author: Jonas Lindemann
"""

import calfem.core as cfc
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_mpl as cfv
import calfem.utils as cfu

import numpy as np
from math import *

# ----- Problem parameters

l = 5.0
h = 1.0
t = 0.2

v = 0.35
E = 2.1e9
ptype = 1
ep = [ptype,t]
D=cfc.hooke(ptype, E, v)

left_support = 10
right_support = 20
top_line = 30

# ----- Define geometry

g = cfg.Geometry()

g.point([0.0, 0.0], marker = left_support) # point 0
g.point([l, 0.0], marker = right_support) # point 1
g.point([l, h]) # point 2
g.point([0.0, h]) # point 2

g.spline([0, 1]) # line 0
g.spline([1, 2]) # line 1
g.spline([2, 3], marker = top_line) # line 2
g.spline([3, 0]) # line 2

g.surface([0, 1, 2, 3])

# ----- Create mesh

mesh = cfm.GmshMesh(g)

mesh.elType = 3 # Degrees of freedom per node.
mesh.dofsPerNode = 2 # Factor that changes element sizes.
mesh.elSizeFactor = 0.10

coords, edof, dofs, bdofs, elementmarkers = mesh.create()

# ----- Solve problem

nDofs = np.size(dofs)
ex, ey = cfc.coordxtr(edof, coords, dofs)

K = np.zeros([nDofs,nDofs])

for eltopo, elx, ely in zip(edof, ex, ey):
    Ke = cfc.planqe(elx, ely, ep, D)
    cfc.assem(eltopo, K, Ke)

bc = np.array([],'i')
bcVal = np.array([],'f')

bc, bcVal = cfu.applybc(bdofs, bc, bcVal, left_support, 0.0, 0)
bc, bcVal = cfu.applybc(bdofs, bc, bcVal, right_support, 0.0, 2)

f = np.zeros([nDofs,1])

cfu.applyforcetotal(bdofs, f, top_line, -10e5, 2)

a,r = cfc.solveq(K,f,bc,bcVal)

ed = cfc.extractEldisp(edof,a)
vonMises = []

for i in range(edof.shape[0]): 
    es, et = cfc.planqs(ex[i,:], ey[i,:], ep, D, ed[i,:]) 
    vonMises.append( sqrt( pow(es[0],2) - es[0]*es[1] + pow(es[1],2) + 3*es[2] ) ) 
    
# ----- Draw geometry

cfv.draw_geometry(g)

# ----- Draw the mesh.

cfv.figure() 
cfv.draw_mesh(
    coords=coords, 
    edof=edof, 
    dofs_per_node=mesh.dofsPerNode, 
    el_type=mesh.elType, 
    filled=True, 
    title="Example 01"
    ) 

# ----- Draw results

cfv.figure()
cfv.draw_element_values(vonMises, coords, edof, mesh.dofs_per_node, mesh.el_type, a, draw_elements=True, draw_undisplaced_mesh=False, title="Example 06 effective stress")

cfv.figure()
cfv.draw_displacements(a, coords, edof, mesh.dofs_per_node, mesh.el_type, draw_undisplaced_mesh=True, title="Example 06", magnfac=10.0)

cfv.showAndWait()