import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_mpl as cfv
import calfem.utils as cfu
import calfem.core as cfc

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.tri as tri

import numpy as np

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

g.spline([1,2], marker=2, elOnCurve=4)
g.spline([3,4], elOnCurve=4)
g.circle([1,0,3], elOnCurve = 10)
g.circle([2,0,4], elOnCurve = 10)
g.spline([3,5], elOnCurve = 6)
g.spline([5,6], marker=3, elOnCurve = 4)
g.spline([6,4], elOnCurve = 6)

# Add Surfaces:
#
# When we set markers for surfaces, and have 2D elements, we can find which 
# region an element is in via the list 'elementmarkers', which is returned by 
# GmshMesher.create()

g.structuredSurface([0,2,1,3], marker = 10)
g.structuredSurface([1,4,5,6], marker = 11)

el_type = 16 
dofs_per_node = 1 

mesh = cfm.GmshMeshGenerator(g, el_type, dofs_per_node) 
coords, edof, dofs, bdofs, elementmarkers = mesh.create()


# ---- Visualise results ----------------------------------------------------

print("Visualising...")

mpl.rcParams['figure.dpi'] = 160

cfv.figure()

cfv.draw_geometry(g, title="Geometry")

cfv.point_in_geometry(g, [0.0, 0.0])

cfv.show_and_wait()

print("Done.")
