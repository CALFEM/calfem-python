#!/bin/env python

import numpy as np
import calfem.core as cfc
import calfem.utils as cfu
import calfem.vis as cfv
import calfem.mesh as cfm


# ---- Problem constants

kx = 50
ky = 50
t = 1.0
ep = [t]

D = np.matrix([
    [kx, 0.],
    [0., ky]
])

# ---- Problem geometry

vertices = np.array([
    [0.0, 0.0],
    [200., 0.0],
    [200., 70.0],
    [120.0, 70.0],
    [120.0, 20.0],
    [80.0, 20.0],
    [80.0, 70.0],
    [0.0, 70.0]
])

segments = np.array([
    [0,1,1],
    [1,2,1],
    [2,3,2],
    [3,4,1],
    [4,5,1],
    [5,6,1],
    [6,7,3],
    [7,0,1]
])

# ---- Create element mesh

print("Creating element mesh...")

coords, edof, dofs, bdofs = cfm.trimesh2d(vertices, segments, maxArea=20.0, dofs_per_node=1)

# ---- Assemble system matrix

print("Assemblig system matrix...")

nDofs = np.size(dofs)
ex, ey = cfc.coordxtr(edof, coords, dofs)

K = np.zeros([nDofs,nDofs])

for eltopo, elx, ely in zip(edof, ex, ey):
    Ke = cfc.flw2te(elx, ely, ep, D)
    cfc.assem(eltopo, K, Ke)

# ---- Solving equation system

print("Solving equation system...")

f = np.zeros([nDofs,1])

bc = np.array([],'i')
bcVal = np.array([],'i')

bc, bcVal = cfu.applybc(bdofs,bc,bcVal,2,30.0)
bc, bcVal = cfu.applybc(bdofs,bc,bcVal,3,0.0)

a, r = cfc.solveq(K,f,bc,bcVal)

# ---- Compute element forces

print("Computing element forces...")

ed = cfc.extractEldisp(edof,a)
qs, qt = cfc.flw2ts(ex, ey, D, ed)

# ---- Visualise results

print("Drawing element mesh...")

cfv.eliso2(ex,ey,ed)    
cfv.eldraw2(ex,ey)
cfv.showAndWait()

print("Done.")




