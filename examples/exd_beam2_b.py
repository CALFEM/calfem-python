# example exd_beam2_b
# ----------------------------------------------------------------
# PURPOSE
#    Structural Dynamics, time integration, time dependent
#    boundary conditions.
#
#    Note: file exd_beam2_m.py must be in the same directory
# ----------------------------------------------------------------

from exd_beam2_m import *

# ----- Impact, center point, vertical beam ----------------------

dt = 0.002
T = 1

# ----- Boundary condition, initial condition --------------------

G = np.array([[0, 0], [0.1, 0.02], [0.2, -0.01], [0.3, 0], [T, 0]])

t, g = cfc.gfunc(G, dt)

bc = np.zeros((4, 1 + len(g)))
bc[0, :] = np.hstack((1, g))
bc[1, 0] = 2
bc[2, 0] = 3
bc[3, 0] = 14

a0 = np.zeros((15, 1))
da0 = np.zeros((15, 1))

# ----- Output parameters ----------------------------------------

times = np.arange(0.1, 1.1, 0.1)
dofs = np.array([1, 4, 11])

# ----- Time integration parameters ------------------------------

ip = np.array([dt, T, 0.25, 0.5])

# ----- Time integration -----------------------------------------

sol, dofhist = cfc.step2(K, [], M, [], a0, da0, bc, ip, times, dofs)

# ----- Plot time history for two DOFs ---------------------------

cfv.figure(1, fig_size=(7, 4))
cfv.plt.plot(t, dofhist["a"][0, :], "-")
cfv.plt.plot(t, dofhist["a"][1, :], "--")
cfv.plt.plot(t, dofhist["a"][2, :], "-.")
cfv.plt.xlim([0, 1])
cfv.plt.ylim([-0.02, 0.03])
cfv.plt.xlabel("time (sec)")
cfv.plt.ylabel("displacement (m)")
cfv.plt.title("Displacement(time) at the 1st, 4th and 11th degree of freedom")
cfv.text("solid line = bottom, vertical beam, x-direction", [0.2, 0.022])
cfv.text("dashed line = center, vertical beam, x-direction", [0.2, 0.017])
cfv.text("dashed-dotted line = center, horizontal beam, y-direction", [0.2, 0.012])
cfv.plt.grid()

# ----- Plot displacement for some time increments ----------------

cfv.figure(2, fig_size=(7, 5))
for i in range(5):
    Edb = cfc.extract_ed(edof, sol["a"][:, i])
    ext = ex + i * 3.5
    cfv.eldraw2(ext, ey, [2, 3, 1])
    cfv.eldisp2(ext, ey, Edb, [1, 2, 2], sfac=20)
    cfv.text(f"{times[i]:.1f}", [3.5 * i + 0.5, 1.5])
eyt = ey - 4
for i in range(5, 10):
    Edb = cfc.extract_ed(edof, sol["a"][:, i])
    ext = ex + (i - 5) * 3.5
    cfv.eldraw2(ext, eyt, [2, 3, 1])
    cfv.eldisp2(ext, eyt, Edb, [1, 2, 2], sfac=20)
    cfv.text(f"{times[i]:.1f}", [3.5 * (i - 5) + 0.5, -2.5])
cfv.title("Snapshots (sec), magnification = 20")
ax = cfv.gca()
ax.set_axis_off()
cfv.show_and_wait()

# ----- End -------------------------------------------------------
