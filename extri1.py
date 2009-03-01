#!/bin/env python

#from numpy import *
from pycalfem import *
from pycalfem_utils import *

# ---- Problem constants

kx = 50
ky = 50
t = 1.0
ep = [t]

D = matrix([
    [kx, 0.],
    [0., ky]
])

# ---- Problem geometry

vertices = array([
    [0.0, 70.0],
    [110.0, 70.0],
    [220.0, 70.0],
    [200.0, 40.0],
    [170.0, 10.0],
    [130.0, 0.0],
    [90.0, 0.0],
    [50.0, 10.0],
    [20.0, 40.0]
])

segments = array([
    [0,1,2],
    [1,2,3],
    [2,3,1],
    [3,4,1],
    [4,5,1],
    [5,6,1],
    [6,7,1],
    [7,8,1],
    [8,0,1]
])

# ---- Create element mesh

print "Creating element mesh..."

coords, edof, dofs, bdofs = trimesh2d(vertices, segments, maxArea=20.0, dofsPerNode=1)

# ---- Assemble system matrix

print "Assemblig system matrix..."

nDofs = size(dofs)
ex, ey = coordxtr(edof, coords, dofs)

K = zeros([nDofs,nDofs])

for etopo, eex, eey in zip(edof, ex, ey):
    Ke = flw2te(eex, eey, ep, D)
    assem(etopo, K, Ke)

# ---- Solving equation system

print "Solving equation system..."

f = zeros([nDofs,1])

bc = array([],'i')
bcVal = array([],'i')

bc, bcVal = applybc(bdofs,bc,bcVal,2,30.0)
bc, bcVal = applybc(bdofs,bc,bcVal,3,0.0)
            
a,r = solveq(K,f,bc,bcVal)

print a

print "Drawing element mesh..."
eldraw2(ex, ey)
show()




