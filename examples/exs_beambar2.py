# -*- coding: utf-8 -*-
#
# example exs7
# ----------------------------------------------------------------
# PURPOSE
#    Set up a frame, consisting of both beams and bars, and
#    illustrate the calculations by use of graphics functions.
# ----------------------------------------------------------------

# REFERENCES
#     P-A Hansson  1994-01-20
#     K-G Olsson   1995-09-28
#     O Dahlblom   2004-10-07
#     J Lindemann  2021-12-29 (Python version)
# ----------------------------------------------------------------

import numpy as np
import calfem.core as cfc
import calfem.utils as cfu
import calfem.vis_mpl as cfv

np.set_printoptions(precision=3, suppress=True)

# ----- System matrices ------------------------------------------

K = np.zeros((18, 18))
f = np.zeros((18, 1))

f[12] = 1

coord = np.array([
    [0., 0.],
    [1., 0.],
    [0., 1.],
    [1., 1.],
    [0., 2.],
    [1., 2.]
])

dof = np.array([
    [1,  2,  3],
    [4,  5,  6],
    [7,  8,  9],
    [10, 11, 12],
    [13, 14, 15],
    [16, 17, 18]
])

# ----- Element properties, topology and coordinates -------------

ep1 = [1, 1, 1]
edof1 = np.array([
    [1,   2,   3,   7,   8,   9],
    [7,   8,   9,   13,  14,  15],
    [4,   5,   6,   10,  11,  12],
    [10,  11,  12,  16,  17,  18],
    [7,   8,   9,   10,  11,  12],
    [13,  14,  15,  16,  17,  18]
])

ex1, ey1 = cfc.coordxtr(edof1, coord, dof, 2)

ep2 = [1, 1]
edof2 = np.array([
    [1,   2,  10,  11],
    [7,   8,  16,  17],
    [7,   8,   4,   5],
    [13,  14,  10,  11]
])

ex2, ey2 = cfc.coordxtr(edof2, coord, dof, 2)

# ----- Draw the fe-mesh as a check of the model -----------------

cfv.figure(1)
cfv.eldraw2(ex1, ey1, [1, 3, 1])
cfv.eldraw2(ex2, ey2, [1, 2, 1])

# ----- Create and assemble element matrices ---------------------

for elx, ely, eltopo in zip(ex1, ey1, edof1):
    Ke = cfc.beam2e(elx, ely, ep1)
    K = cfc.assem(eltopo, K, Ke)

for elx, ely, eltopo in zip(ex2, ey2, edof2):
    Ke = cfc.bar2e(elx, ely, ep2)
    K = cfc.assem(eltopo, K, Ke)

# ----- Solve equation system ------------------------------------

bc = np.array([1, 2, 3, 4, 5, 6])
a, r = cfc.solveq(K, f, bc)

# ---- Extract element displacements and display the deformed mesh -

ed1 = cfc.extract_ed(edof1, a)
ed2 = cfc.extract_ed(edof2, a)

sfac = cfv.scalfact2(ex1, ey1, ed1, 0.1)

cfv.figure(2)
cfv.eldraw2(ex1, ey1)
cfv.eldraw2(ex2, ey2)
cfv.eldisp2(ex1, ey1, ed1, [2, 1, 1], sfac)
cfv.eldisp2(ex2, ey2, ed2, [2, 1, 1], sfac)

cfv.show_and_wait()


# -------------------------- end --------------------------------
