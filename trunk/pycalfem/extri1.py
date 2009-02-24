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

coord, bcVert, topo = triangle(vertices, segments, maxArea=0.01)

#print coord
#print bcVert
#print topo

nDof = 1
nElements = size(topo,0)

dof = arange(nElements*nDof).reshape(nElements,nDof)+1

ex, ey, = coordxtr(topo, coord, dof, 3)