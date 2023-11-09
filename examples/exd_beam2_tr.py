# example exd_beam2_tr
# ----------------------------------------------------------------
# PURPOSE
#    Structural Dynamics, time integration, reduced system.
#
#    Note: file exd_beam2_m.py must be in the same directory
# ----------------------------------------------------------------

from exd_beam2_m import *

# ----- Impact, center point, vertical beam ----------------------

dt = 0.002
T = 1
nev = 2

# ----- The load -------------------------------------------------

G = np.array([[0, 0], [0.15, 1], [0.25, 0], [T, 0]])

t, g = cfc.gfunc(G, dt)
f = np.zeros((15, len(g)))
f[3, :] = 1000 * g
fr = np.hstack((np.arange(1, nev + 1).reshape(-1, 1), Egv[:, :nev].T @ f))

# ----- Reduced system matrices ----------------------------------

kr = np.diag(np.diag(Egv[:, :nev].T @ K @ Egv[:, :nev]))
mr = np.diag(np.diag(Egv[:, :nev].T @ M @ Egv[:, :nev]))

# ----- Initial condition ----------------------------------------

ar0 = np.zeros((nev, 1))
dar0 = np.zeros((nev, 1))

# ----- Output parameters ----------------------------------------

times = np.arange(0.1, 1.1, 0.1)
dofsr = np.arange(1, nev + 1)
dofs = np.array([4, 11])

# ----- Time integration parameters ------------------------------

ip = np.array([dt, T, 0.25, 0.5])

# ----- Time integration -----------------------------------------
sol, dofhist = cfc.step2(kr, [], mr, fr, ar0, dar0, [], ip, times, dofsr)

# ----- Mapping back to original coordinate system ---------------

aR = Egv[:, :nev] @ sol["a"]
aRhist = Egv[dofs - 1, :nev] @ dofhist["a"]

# ----- Plot time history for two DOFs ---------------------------

cfv.figure(1, fig_size=(7, 4))
cfv.plt.plot(t, aRhist[0, :], "-")
cfv.plt.plot(t, aRhist[1, :], "--")
cfv.plt.xlim([0, 1])
cfv.plt.ylim([-0.01, 0.02])
cfv.plt.xlabel("time (sec)")
cfv.plt.ylabel("displacement (m)")
cfv.plt.title("Displacement(time) at the 4th and 11th degree of freedom")
cfv.text("solid line = impact point, x-direction", [0.3, 0.017])
cfv.text("dashed line = center, horizontal beam, y-direction", [0.3, 0.012])
cfv.text("TWO EIGENVECTORS ARE USED", [0.3, -0.007])
cfv.plt.grid()

# ----- Plot displacement for some time increments ----------------

cfv.figure(2, fig_size=(7, 5))

for i in range(5):
    Edb = cfc.extract_ed(edof, aR[:, i])
    ext = ex + i * 3.5
    cfv.eldraw2(ext, ey, [2, 3, 1])
    cfv.eldisp2(ext, ey, Edb, [1, 2, 2], sfac=25)
    cfv.text(f"{times[i]:.1f}", [3.5 * i + 0.5, 1.5])

eyt = ey - 4

for i in range(5, 10):
    Edb = cfc.extract_ed(edof, aR[:, i])
    ext = ex + (i - 5) * 3.5
    cfv.eldraw2(ext, eyt, [2, 3, 1])
    cfv.eldisp2(ext, eyt, Edb, [1, 2, 2], sfac=25)
    cfv.text(f"{times[i]:.1f}", [3.5 * (i - 5) + 0.5, -2.5])

cfv.title("Snapshots (sec), magnification = 25")
ax = cfv.gca()
ax.set_axis_off()
cfv.show_and_wait()

# ----- End -------------------------------------------------------
