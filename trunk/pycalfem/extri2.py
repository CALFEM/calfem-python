#!/bin/env python

from pycalfem import *
from pycalfem_utils import *
from pycalfem_classes import *

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

print "Creating element mesh..."

coords, edof, dofs, bdofs = trimesh2d(vertices, segments, maxArea=0.00005, dofsPerNode=2)

# ---- Assemble system matrix

print "Assemblig system matrix..."

nDofs = size(dofs)
ex, ey = coordxtr(edof, coords, dofs)

#eldraw2(ex,ey)

K = zeros([nDofs,nDofs])

for eltopo, elx, ely in zip(edof, ex, ey):
    Ke = plante(elx, ely, ep, D)
    assem(eltopo, K, Ke)

# ---- Solving equation system

print "Solving equation system..."

f = zeros([nDofs,1])

bc = array([],'i')
bcVal = array([],'i')

bc, bcVal = applybc(bdofs,bc,bcVal,2,0.0)
applyforce(bdofs,f,3,-10e3,2)
            
a,r = solveq(K,f,bc,bcVal)

# ---- Compute element forces

print "Computing element forces..."

ed = extract(edof,a)
es, et = plants(ex, ey, ep, D, ed)

# ---- Visualise results

print "Drawing element mesh..."

eldisp2(ex,ey,ed)

print "Done."




