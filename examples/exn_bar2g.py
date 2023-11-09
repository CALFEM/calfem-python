# example exn_bar2g
# --------------------------------------------------------------------------
# PURPOSE
#    Analysis of a plane truss using second order theory.
# --------------------------------------------------------------------------

import numpy as np
import calfem.core as cfc
import calfem.utils as cfu

# ----- Topology -----

edof = np.array([[1, 2, 5, 6], [3, 4, 5, 6]])

# ----- Element properties and global coordinates -------------------------

E = 10e9
A1 = 4e-2
A2 = 1e-2
ep1 = np.array([E, A1])
ep2 = np.array([E, A2])

ex1 = np.array([0.0, 1.6])
ey1 = np.array([0.0, 0.0])
ex2 = np.array([0.0, 1.6])
ey2 = np.array([1.2, 0])

# ----- Initial values for the iteration ----------------------------------

eps = 1e-6  # Error norm
QX1 = 0.01
QX2 = 0
# Initial axial forces
QX01 = 1  # Axial force of the initial former iteration
n = 0  # Iteration counter

# ----- Iteration procedure -----------------------------------------------
while abs((QX1 - QX01) / QX01) > eps:
    n += 1

    K = np.zeros((6, 6))
    f = np.zeros((6, 1))
    f[4] = -10e6
    f[5] = -0.2e6

    Ke1 = cfc.bar2ge(ex1, ey1, ep1, QX1)
    Ke2 = cfc.bar2ge(ex2, ey2, ep2, QX2)
    K = cfc.assem(edof[0, :], K, Ke1)
    K = cfc.assem(edof[1, :], K, Ke2)
    bc = np.array([1, 2, 3, 4])
    a, r = cfc.solveq(K, f, bc)

    Ed = cfc.extract_ed(edof, a)

    QX01 = QX1
    es1, QX1 = cfc.bar2gs(ex1, ey1, ep1, Ed[0, :])
    es2, QX2 = cfc.bar2gs(ex2, ey2, ep2, Ed[1, :])

    if n > 20:
        print("The solution does not converge")
        break

# ----- Results -----------------------------------------------

cfu.disp_h1("Displacements:")
cfu.disp_array(a)
cfu.disp_h1("Normal forces:")
cfu.disp(QX1)
cfu.disp(QX2)
