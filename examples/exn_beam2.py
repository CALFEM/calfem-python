# example exn_beam2g
# ----------------------------------------------------------------
# PURPOSE
#    Geometrically nonlinear analysis of a plane frame.
# ----------------------------------------------------------------

import numpy as np
import calfem.core as cfc
import calfem.vis_mpl as cfv


# ----- Topology -------------------------------------------------

edof = np.array([
    [4, 5, 6, 1, 2, 3], 
    [7, 8, 9, 10, 11, 12], 
    [4, 5, 6, 7, 8, 9]
])

# ----- Element stiffness and element load matrices  -------------

E = 200e9
A1 = 2e-3
A2 = 6e-3
I1 = 1.6e-5
I2 = 5.4e-5

ep1 = np.array([E, A1, I1])
ep3 = np.array([E, A2, I2])

ex1 = np.array([0.0, 0.0])
ey1 = np.array([4.0, 0.0])
ex2 = np.array([6.0, 6.0])
ey2 = np.array([4.0, 0.0])
ex3 = np.array([0.0, 6.0])
ey3 = np.array([4.0, 4.0])

eq1 = np.array([0.0])
eq2 = np.array([0.0])
eq3 = np.array([-50e3])

# ----- Initial axial forces ----------------------------------------------

QX1 = 1e-4
QX2 = 0
QX3 = 0
QX01 = 1

# ----- Iteration for convergence -----------------------------------------

eps = 1e-6
n = 0

while abs((QX1 - QX01) / QX01) > eps:
    n += 1
    K = np.zeros([12, 12])
    f = np.zeros([12, 1])
    f[3, 0] = 10e3

    Ke1 = cfc.beam2ge(ex1, ey1, ep1, QX1)
    Ke2 = cfc.beam2ge(ex2, ey2, ep1, QX2)
    Ke3, fe3 = cfc.beam2ge(ex3, ey3, ep3, QX3, eq3)

    K = cfc.assem(edof[0, :], K, Ke1)
    K = cfc.assem(edof[1, :], K, Ke2)
    K, f = cfc.assem(edof[2, :], K, Ke3, f, fe3)

    bc = np.array([1, 2, 3, 10, 11])
    a, r = cfc.solveq(K, f, bc)

    ed = cfc.extract_ed(edof, a)

    QX01 = QX1
    es1 = cfc.beam2gs(ex1, ey1, ep1, ed[0, :], QX1, eq1)
    es2 = cfc.beam2gs(ex2, ey2, ep1, ed[1, :], QX2, eq2)
    es3 = cfc.beam2gs(ex3, ey3, ep3, ed[2, :], QX3, eq3)
    QX1 = es1[1]
    QX2 = es2[1]
    QX3 = es3[1]

    if n == 1:
        ed0 = ed

    if n > 20:
        print("The solution does not converge")
        break


# ----- Section forces ---------------------------------------

eq1 = np.array([0.0, eq1.item()])
eq2 = np.array([0.0, eq2.item()])
eq3 = np.array([0.0, eq3.item()])
es1, edi1, eci1 = cfc.beam2s(ex1, ey1, ep1, ed[0, :], eq1, 21)
es2, edi2, eci2 = cfc.beam2s(ex2, ey2, ep1, ed[1, :], eq2, 21)
es3, edi3, eci3 = cfc.beam2s(ex3, ey3, ep3, ed[2, :], eq3, 21)


# ----- Draw deformed frame ---------------------------------------

cfv.figure(1, fig_size=(6, 4))
plotpar = [3, 1, 0]
cfv.eldraw2(ex1, ey1, plotpar)
cfv.eldraw2(ex2, ey2, plotpar)
cfv.eldraw2(ex3, ey3, plotpar)
sfac = cfv.scalfact2(ex3, ey3, edi3, 0.1)
plotpar = [1, 2, 1]
cfv.eldisp2(ex1, ey1, ed[0, :], plotpar, sfac)
cfv.eldisp2(ex2, ey2, ed[1, :], plotpar, sfac)
cfv.eldisp2(ex3, ey3, ed[2, :], plotpar, sfac)
plotpar = [2, 4, 2]
cfv.eldisp2(ex1, ey1, ed0[0, :], plotpar, sfac)
cfv.eldisp2(ex2, ey2, ed0[1, :], plotpar, sfac)
cfv.eldisp2(ex3, ey3, ed0[2, :], plotpar, sfac)
cfv.plt.axis([-1.5, 7.5, -0.5, 5.5])
cfv.title("Displacements")


# ----- Draw normal force diagram --------------------------------
cfv.figure(2, fig_size=(6, 4))
plotpar = [2, 1]
sfac = cfv.scalfact2(ex1, ey1, es1[:, 0], 0.2)
cfv.secforce2(ex1, ey1, es1[:, 0], plotpar, sfac)
cfv.secforce2(ex2, ey2, es2[:, 0], plotpar, sfac)
cfv.secforce2(ex3, ey3, es3[:, 0], plotpar, sfac)
cfv.plt.axis([-1.5, 7.5, -0.5, 5.5])
cfv.title("Normal force")

# ----- Draw shear force diagram ---------------------------------
cfv.figure(3, fig_size=(6, 4))
plotpar = [2, 1]
sfac = cfv.scalfact2(ex3, ey3, es3[:, 1], 0.2)
cfv.secforce2(ex1, ey1, es1[:, 1], plotpar, sfac)
cfv.secforce2(ex2, ey2, es2[:, 1], plotpar, sfac)
cfv.secforce2(ex3, ey3, es3[:, 1], plotpar, sfac)
cfv.plt.axis([-1.5, 7.5, -0.5, 5.5])
cfv.title("Shear force")

# ----- Draw moment diagram --------------------------------------
cfv.figure(4, fig_size=(6, 4))
plotpar = [2, 1]
sfac = cfv.scalfact2(ex3, ey3, es3[:, 2], 0.2)
cfv.secforce2(ex1, ey1, es1[:, 2], plotpar, sfac)
cfv.secforce2(ex2, ey2, es2[:, 2], plotpar, sfac)
cfv.secforce2(ex3, ey3, es3[:, 2], plotpar, sfac)
cfv.plt.axis([-1.5, 7.5, -0.5, 5.5])
cfv.title("Bending moment")
cfv.show_and_wait()
# ------------------------ end -----------------------------------
