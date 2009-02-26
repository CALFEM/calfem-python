#!/bin/env python

from numpy import *
from pycalfem import *
from pycalfem_utils import *

vertices = array([
    [0.0, 0.0],
    [1.0, 0.0],
    [1.0, 1.0],
    [0.0, 1.0]
])

segments = array([
    [0,1,1],
    [1,2,1],
    [2,3,1],
    [3,0,1]
])

coords, bcVert, edof = triangle(vertices, segments, maxArea=0.005)
dofs = createdofs(coords,1)
ex, ey = coordxtr(edof, coords, dofs)

eldraw2(ex, ey)
show()

#print edof
