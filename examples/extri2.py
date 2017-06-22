#!/bin/env python

from calfem.core import *
from calfem.utils import *

def drawCustom(self, width, height):
    glPushMatrix()
    glBegin(GL_LINES)
    glColor(1.0, 0.0, 0.0, 1.0)
    glVertex(50,50,0)
    glVertex(100,100,0)
    glPopMatrix()

# ---- Problem constants

t = 0.1
v = 0.35
E = 2.1e9
ptype = 1
ep = [ptype,t]

D=hooke(ptype, E, v)

# ---- Problem geometry

l = 0.2
w = 0.05
h = 0.1

vertices = array([
    [0.0, h],
    [l, h],
    [l, 0.0],
    [l-w, 0.0],
    [l-w, h-w],
    [0.0, h-w]
])

segments = array([
    [0,1,1],
    [1,2,1],
    [2,3,3],
    [3,4,1],
    [4,5,1],
    [5,0,2]
])

# ---- Create element mesh

print("Creating element mesh...")

coords, edof, dofs, bdofs = trimesh2d(vertices, segments, maxArea=0.00005, dofsPerNode=2)

# ---- Assemble system matrix

print("Assemblig system matrix...")

nDofs = size(dofs)
ex, ey = coordxtr(edof, coords, dofs)

#eldraw2(ex,ey)

K = zeros([nDofs,nDofs])

for eltopo, elx, ely in zip(edof, ex, ey):
    Ke = plante(elx, ely, ep, D)
    assem(eltopo, K, Ke)

# ---- Solving equation system

print("Solving equation system...")

f = zeros([nDofs,1])

bc = array([],'i')
bcVal = array([],'i')

bc, bcVal = applybc(bdofs,bc,bcVal,2,0.0)
bc, bcVal = applybcnode(0, dofs, bc, bcVal, 0.0, 1)
bc, bcVal = applybcnode(0, dofs, bc, bcVal, 0.0, 2)
bc, bcVal = applybcnode(5, dofs, bc, bcVal, 0.0, 1)
bc, bcVal = applybcnode(5, dofs, bc, bcVal, 0.0, 2)

applyforce(bdofs,f,3,-10e3,2)
            
a,r = solveq(K,f,bc,bcVal)

# ---- Compute element forces

print("Computing element forces...")

ed = extractEldisp(edof,a)
es, et = plants(ex, ey, ep, D, ed)
ev = effmises(es, ptype)
esnv = stress2nodal(ev, edof)

# ---- Visualise results

print("Drawing element mesh...")

eldisp2(ex, ey, ed)
elval2(ex, ey, ev)
eliso2(ex, ey, esnv)

elementView = ElementView(None, -1, "")
elementView.ex = ex
elementView.ey = ey
elementView.ev = ev
elementView.showMesh = False
elementView.showElementValues = True
elementView.showNodalValues = False
elementView.drawCustom = drawCustom
elementView.Show()

waitDisplay()

print("Done.")



