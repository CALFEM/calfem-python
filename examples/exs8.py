# -*- coding: utf-8 -*-

# example exs8
#----------------------------------------------------------------
# PURPOSE 
#    Analysis of two dimensional diffusion
#----------------------------------------------------------------

# REFERENCES
#     Karl-Gunnar Olsson 1995-10-08
#     Ola Dahlblom 2004-09-14
#----------------------------------------------------------------

import numpy as np
import calfem.vis_mpl as cfv
import calfem.core as cfc

# ----- System matrices -----

K = np.zeros((15,15))
f = np.zeros((15,1))
Coord = np.array([
    [0,     0    ],[0.025, 0    ],
    [0.05,  0    ],[0,     0.025],
    [0.025, 0.025],[0.05,  0.025],
    [0,     0.05 ],[0.025, 0.05 ],
    [0.05,  0.05 ],[0,     0.075],
    [0.025, 0.075],[0.05,  0.075],
    [0,     0.1  ],[0.025, 0.1  ],
    [0.05,  0.1  ]
])

Dof = np.array([
    [1 ],[2 ],[3 ],
    [4 ],[5 ],[6 ],
    [7 ],[8 ],[9 ],
    [10],[11],[12],
    [13],[14],[15]
])

# ----- Element properties, topology and coordinates -----

ep = np.array([1])
D = np.array([
    [1, 0],
    [0, 1]
])
Edof = np.array([
    [ 1, 2, 5, 4],
    [ 2, 3, 6, 5],
    [ 4, 5, 8, 7],
    [ 5, 6, 9, 8],
    [ 7, 8,11,10],
    [ 8, 9,12,11],
    [10,11,14,13],
    [11,12,15,14],
])
Ex,Ey = cfc.coordxtr(Edof,Coord,Dof)

# ----- Generate FE-mesh -----

#clf; eldraw2(Ex,Ey,[1 3 0],Edof(:,1));
#disp('PRESS ENTER TO CONTINUE'); pause; clf;

# ----- Create and assemble element matrices -----

for i in range(8):
    Ke = cfc.flw2qe(Ex[i],Ey[i],ep,D)
    K = cfc.assem(Edof[i],K,Ke)

# ----- Solve equation system -----

bcPrescr = np.array([1,2,3,4,7,10,13,14,15])
bcVal = np.array([0,0,0,0,0,0,0.5e-3,1e-3,1e-3])
a,r = cfc.solveq(K,f,bcPrescr,bcVal)

# ----- Compute element flux vector -----

Ed = cfc.extractEldisp(Edof,a)
Es = np.zeros((8,2))
for i in range(8):
    Es[i],Et = cfc.flw2qs(Ex[i],Ey[i],ep,D,Ed[i])

# ----- Draw flux vectors and contourlines -----

print(Ex)
print(Ey)
print(a)
print(Ed)

cfv.eldraw2(Ex, Ey, [1, 2, 1], range(1,Ex.shape[0]+1))
cfv.eliso2_mpl(Ex,Ey,Ed);
cfv.showAndWaitMpl()

#cfv.showAndWait()
#sfac=scalfact2(Ex,Ey,Es,0.5);
#eldraw2(Ex,Ey); 
#elflux2(Ex,Ey,Es,[1,4],sfac); 
#pltscalb2(sfac,[2e-2 0.06 0.01],4);
#disp('PRESS ENTER TO CONTINUE'); pause; clf;
#eldraw2(Ex,Ey,[1,3,0]); 
#eliso2(Ex,Ey,Ed,5,[1,4]);
#hold off; 
#echo off;

# ----------------- End --------------------------------

