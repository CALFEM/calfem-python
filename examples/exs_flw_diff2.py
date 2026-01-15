# -*- coding: utf-8 -*-

# example exs8
# ----------------------------------------------------------------
# PURPOSE
#    Analysis of two dimensional diffusion
# ----------------------------------------------------------------

# REFERENCES
#     Karl-Gunnar Olsson 1995-10-08
#     Ola Dahlblom 2004-09-14
# ----------------------------------------------------------------

import numpy as np
import calfem.vis_mpl as cfv
import calfem.core as cfc
import calfem.utils as cfu

# ----- System matrices -----

coord = np.array([
    [0, 0],
    [0.025, 0],
    [0.05, 0],
    [0, 0.025],
    [0.025, 0.025],
    [0.05, 0.025],
    [0, 0.05],
    [0.025, 0.05],
    [0.05, 0.05],
    [0, 0.075],
    [0.025, 0.075],
    [0.05, 0.075],
    [0, 0.1],
    [0.025, 0.1],
    [0.05, 0.1],
])

dofs = np.arange(1, coord.shape[0] + 1).reshape(coord.shape[0], 1)

# ----- Element properties, topology and coordinates -----

ep = [1]

D = np.array([
    [1, 0], 
    [0, 1]
])

edof = np.array([
    [1, 2, 5, 4],
    [2, 3, 6, 5],
    [4, 5, 8, 7],
    [5, 6, 9, 8],
    [7, 8, 11, 10],
    [8, 9, 12, 11],
    [10, 11, 14, 13],
    [11, 12, 15, 14],
])

ex, ey = cfc.coordxtr(edof, coord, dofs)

K = np.zeros((15, 15))  # Global conductivity matrix
f = np.zeros((15, 1))   # Global source vector (no sources in this problem)

# ----- Generate FE-mesh -----

# clf; eldraw2(Ex,Ey,[1 3 0],Edof(:,1));
# disp('PRESS ENTER TO CONTINUE'); pause; clf;

# ----- Create and assemble element matrices -----

for elx, ely, etopo in zip(ex, ey, edof):
    Ke = cfc.flw2qe(elx, ely, ep, D)
    K = cfc.assem(etopo, K, Ke)

# for i in range(8):
#     Ke = cfc.flw2qe(ex[i], ey[i], ep, D)
#     K = cfc.assem(edof[i], K, Ke)

# ----- Solve equation system -----

bc_dofs = np.array([1, 2, 3, 4, 7, 10, 13, 14, 15])
bc_vals = np.array([0, 0, 0, 0, 0, 0, 0.5e-3, 1e-3, 1e-3])

a, r = cfc.solveq(K, f, bc_dofs, bc_vals)

cfu.disp_h2("Concentration at nodes [kg/m^3]:")
cfu.disp_array(a, ["a"])

cfu.disp_h2("Boundary fluxes at nodes [kg/m^2/s)]:")
cfu.disp_array(r, ["r"])


# ----- Compute element flux vector -----

ed = cfc.extract_eldisp(edof, a)

es = np.zeros((8, 2))
el_idx = 0

for elx, ely, eld in zip(ex, ey, ed):
    es_el, t = cfc.flw2qs(elx, ely, ep, D, eld)
    es[el_idx] = es_el

    el_idx += 1

# for i in range(8):
#     es[i], et = cfc.flw2qs(ex[i], ey[i], ep, D, ed[i])

cfu.disp_h2(f"Element flux vectors [kg/m^2/s]:")
cfu.disp_array(es, headers=["qx", "qy"])

cfu.disp_h2("Concentration field [*10^-3 kg/m^3]:")
cfu.disp("Pure water boundaries (DOFs 1-4,7,10): 0.000")
cfu.disp(f"Internal concentrations:")

a_dofs = [5, 6, 8, 9, 11, 12]
a_internal = a[a_dofs - np.ones(len(a_dofs), dtype=int)]
a_table = np.hstack((np.array(a_dofs, dtype=int).reshape(-1, 1), a_internal.reshape(-1, 1)))

cfu.disp_array(a_table, headers=["DOF", "Concentration"])

cfu.disp(f"Solution boundaries (DOFs 14,15): 1.000")

# ----- Draw flux vectors and contourlines -----

cfv.figure()
cfv.eldraw2(ex, ey, [1, 2, 1], range(1, ex.shape[0] + 1))
cfv.eliso2_mpl(ex, ey, ed)
cfv.show_and_wait()
