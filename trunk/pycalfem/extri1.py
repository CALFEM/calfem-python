#!/bin/env python

from numpy import *
from pycalfem import *
from pycalfem_utils import *

vertices = array([
    [0.0, 0.0],
    [1.0, 0.0],
    [1.0, 2.0],
    [0.0, 1.0]
])

segments = array([
    [0,1,1],
    [1,2,2],
    [2,3,3],
    [3,0,4]
])


coords, edof, dofs, bcVert = trimesh2d(vertices, segments, maxArea=0.01, dofsPerNode=2)
ex, ey = coordxtr(edof, coords, dofs)

kx = 50
ky = 50
t = 1.0
ep = [t]

D = matrix([
    [kx, 0.],
    [0., ky]
])

Ke = flw2te(ex[0,:],ey[0,:],ep,D)
print Ke

#eldraw2(ex, ey)
#show()




