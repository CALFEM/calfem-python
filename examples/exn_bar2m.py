# example exn_bar2m
# --------------------------------------------------------------------------
# PURPOSE
#    Analysis of a plane truss considering material nonlinearity.
# --------------------------------------------------------------------------

import numpy as np
import calfem.core as cfc
import calfem.utils as cfu
import calfem.vis_mpl as cfv

edof = np.array([
    [1, 2, 5, 6], 
    [5, 6, 7, 8], 
    [3, 4, 5, 6]
])

bc = np.array([1, 2, 3, 4, 7, 8])

ex = np.array([
    [0.0, 1.6], 
    [1.6, 1.6], 
    [0.0, 1.6]
])

ey = np.array([[0.0, 0.0], [0.0, 1.2], [1.2, 0.0]])

E = np.array([200e9, 200e9, 200e9])
A = np.array([6.0e-4, 3.0e-4, 10.0e-4])
SY = 400e6
Ns = (SY * A).reshape(-1, 1)

dp = 4e3

incr = 100

a = np.zeros((8, 1))
r = np.zeros((8, 1))
es = np.zeros((3, 1))

plbar = 0
pl = np.array([0.0, 0.0])

# Forward Euler increamental solution

for i in range(incr):
    K = np.zeros((8, 8))
    df = np.zeros((8, 1))
    df[5] = -dp

    # Create and assemble element tangent stiffness matrix

    for j in range(3):
        ep = np.array([E[j], A[j]])
        Ke = cfc.bar2e(ex[j, :], ey[j, :], ep)
        K = cfc.assem(edof[j, :], K, Ke)

    # Stop iteration if determinant det(Kr) <= 0

    fdof = np.setdiff1d(np.arange(1, 9), bc) - 1
    Kr = K[np.ix_(fdof, fdof)]
    if np.linalg.det(Kr) <= 0:
        print("Determinant zero after increment ", i)
        break

    # Solve for the displacement increment and determine total displacements

    da, dr = cfc.solveq(K, df, bc)
    a += da
    r += dr

    # Determine normal forces in elements

    ded = cfc.extract_ed(edof, da)
    des = np.zeros((3, 1))

    for j in range(3):
        ep = np.array([E[j], A[j]])
        desj = cfc.bar2s(ex[j, :], ey[j, :], ep, ded[j, :])
        des[j, 0] = desj[0]

    es += des

    for j in range(3):
        if abs(es[j, 0]) >= Ns[j]:
            E[j] = 0

    # Determine if the stress in a bar has reached the yield stress

    newplbar = np.sum(abs(es) > Ns)

    if newplbar > plbar:
        plbar = newplbar
        print(
            plbar, "plastic elements for increment ", i + 1, " at load = ", (i + 1) * dp
        )

    # Save variables for curve plotting

    pl = np.vstack((pl, np.array([-a[5, 0], (i + 1) * dp])))

# Plot force-displacement relation

cfv.figure(1, fig_size=(7, 4))
cfv.plt.plot(pl[:, 0], pl[:, 1])
cfv.plt.xlabel("Displacement")
cfv.plt.ylabel("Force")
cfv.show_and_wait()
