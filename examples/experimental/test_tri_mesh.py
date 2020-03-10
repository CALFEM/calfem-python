# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.tri as tri

import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_mpl as cfv
import calfem.utils as cfu
import calfem.core as cfc

import numpy as np

coords = np.array([
        [0.000, 0.000],
        [1.000, 0.000],
        [2.000, 0.500],
        [0.000, 1.000],
        [1.000, 1.000],
        [1.750, 1.300],
        [1.000, 1.700]])
        
edof = np.array([
        [1, 2, 5],
        [5, 4, 1],
        [2, 3, 6],
        [6, 5, 2],
        [4, 5, 7],
        [5, 6, 7]])

values = [1, 2, 1, 2, 7, 4, 5]

plt.figure()
cfv.draw_nodal_values_contours(values, coords, edof)
plt.colorbar()

plt.figure()
cfv.draw_nodal_values_shaded(values, coords, edof)
plt.colorbar()

plt.show()