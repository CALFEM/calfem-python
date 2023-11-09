# example exd_beam2_m
# ----------------------------------------------------------------
# PURPOSE
#    Set up the fe-model and perform eigenvalue analysis
#    for a simple frame structure.
# ----------------------------------------------------------------

import numpy as np
import calfem.core as cfc
import calfem.vis_mpl as cfv

# ----- Generate the model ---------------------------------------
# ----- Material data --------------------------------------------

E = 3e10
Av = 0.1030e-2
Ah = 0.0764e-2
rho = 2500
Iv = 0.0171e-4
Ih = 0.00801e-4
ep1 = [E, Av, Iv, rho * Av]  # IPE100
ep2 = [E, Ah, Ih, rho * Ah]  # IPE80

# ----- Topology -------------------------------------------------

edof = np.array(
    [
        [1, 2, 3, 4, 5, 6],
        [4, 5, 6, 7, 8, 9],
        [7, 8, 9, 10, 11, 12],
        [10, 11, 12, 13, 14, 15],
    ]
)

# ----- List of coordinates --------------------------------------

coord = np.array(
    [
        [0.0, 0.0],
        [0.0, 1.5],
        [0.0, 3.0],
        [1.0, 3.0],
        [2.0, 3.0],
    ]
)

# ----- List of degrees of freedom -------------------------------

dof = np.array(
    [
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
        [10, 11, 12],
        [13, 14, 15],
    ]
)

# ----- Generate element matrices, assemble in global matrices ---

K = np.zeros([15, 15])
M = np.zeros([15, 15])

ex, ey = cfc.coordxtr(edof, coord, dof)
ep = np.array([ep1, ep1, ep2, ep2])

for elx, ely, eltopo, elprop in zip(ex, ey, edof, ep):
    Ke, Me = cfc.beam2de(elx, ely, elprop)
    cfc.assem(eltopo, K, Ke)
    cfc.assem(eltopo, M, Me)

# ----- Eigenvalue analysis --------------------------------------

b = np.array([1, 2, 3, 14])
La, Egv = cfc.eigen(K, M, b)
freq = np.sqrt(La) / (2 * np.pi)

if __name__ == "__main__":
    # ----- Draw a plot of the element mesh --------------------------

    cfv.figure(1, fig_size=(5.5, 4.5))
    cfv.eldraw2(ex, ey, [1, 2, 1])
    cfv.title("2-D Frame Structure")

    # ----- Plot one eigenmode ---------------------------------------

    cfv.figure(2, fig_size=(5.5, 4.5))
    cfv.eldraw2(ex, ey, [2, 3, 1])
    Edb = cfc.extract_ed(edof, Egv[:, 0])
    cfv.eldisp2(ex, ey, Edb, [1, 2, 2])
    cfv.title("The first eigenmode")
    cfv.text(f"{freq[0]:.2f}", [0.5, 1.75])
    ax = cfv.gca()
    ax.grid()

    # ----- Plot eight eigenmodes ------------------------------------

    cfv.figure(3, fig_size=(7, 5))
    for i in range(4):
        Edb = cfc.extract_ed(edof, Egv[:, i])
        ext = ex + i * 3
        cfv.eldraw2(ext, ey, [2, 3, 1])
        cfv.eldisp2(ext, ey, Edb, [1, 2, 2], sfac=0.5)
        cfv.text(f"{freq[i]:.2f}", [3 * i + 0.5, 1.5])
    eyt = ey - 4
    for i in range(4, 8):
        Edb = cfc.extract_ed(edof, Egv[:, i])
        ext = ex + (i - 4) * 3
        cfv.eldraw2(ext, eyt, [2, 3, 1])
        cfv.eldisp2(ext, eyt, Edb, [1, 2, 2], sfac=0.5)
        cfv.text(f"{freq[i]:.2f}", [3 * (i - 4) + 0.5, -2.5])
    cfv.title("The first eight eigenmodes [Hz]")
    ax = cfv.gca()
    ax.set_axis_off()
    cfv.show_and_wait()

# ----- End -------------------------------------------------------
