# -*- coding: iso-8859-15 -*-

"""
CALFEM Editor Example

Written by Karl Eriksson
"""

import calfem.vis as cfv
import calfem.editor as cfe
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_mpl as cfv
import calfem.utils as cfu
import calfem.core as cfc
import numpy as np

#polydrawing.run()
#print("")
#g = cfu.load_geometry("testGeom")
#cfe.run_and_load(g)

g = cfg.Geometry()

g.point([0.0, 0.0]) # point 0
g.point([150.0, 0.0]) # point 1
g.point([125, 140.0]) # point 2

g.spline([0, 1]) # line 0
g.spline([1, 2]) # line 1
g.spline([2, 0]) # line 2

g.surface([0, 1, 2])
g.setCurveMarker(0,20)
g.setPointMarker(0,"Hej")
g.setPointMarker(1,10)

print(g.points)
print(g.curves)
print(g.marker_dict)
#cfv.draw_geometry(g)
#cfv.show_and_wait_mpl()
#print(g)
g, dict = cfe.edit_geometry(g)
print(g, dict)

cfv.draw_geometry(g)
cfv.show_and_wait_mpl()
