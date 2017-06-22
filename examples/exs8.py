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

from numpy import *
from calfem.core import *

# ----- System matrices -----

K = zeros((15,15))
f = zeros((15,1))
Coord = array([
    [0,     0    ],[0.025, 0    ],
    [0.05,  0    ],[0,     0.025],
    [0.025, 0.025],[0.05,  0.025],
    [0,     0.05 ],[0.025, 0.05 ],
    [0.05,  0.05 ],[0,     0.075],
    [0.025, 0.075],[0.05,  0.075],
    [0,     0.1  ],[0.025, 0.1  ],
    [0.05,  0.1  ]
])

Dof = array([
    [1 ],[2 ],[3 ],
    [4 ],[5 ],[6 ],
    [7 ],[8 ],[9 ],
    [10],[11],[12],
    [13],[14],[15]
])

# ----- Element properties, topology and coordinates -----

ep = array([1])
D = array([
    [1, 0],
    [0, 1]
])
Edof = array([
    [ 1, 2, 5, 4],
    [ 2, 3, 6, 5],
    [ 4, 5, 8, 7],
    [ 5, 6, 9, 8],
    [ 7, 8,11,10],
    [ 8, 9,12,11],
    [10,11,14,13],
    [11,12,15,14],
])
Ex,Ey = coordxtr(Edof,Coord,Dof)

# ----- Generate FE-mesh -----

#clf; eldraw2(Ex,Ey,[1 3 0],Edof(:,1));
#disp('PRESS ENTER TO CONTINUE'); pause; clf;

# ----- Create and assemble element matrices -----

for i in range(8):
    Ke = flw2qe(Ex[i],Ey[i],ep,D)
    K = assem(Edof[i],K,Ke)

# ----- Solve equation system -----

bcPrescr = array([1,2,3,4,7,10,13,14,15])
bcVal = array([0,0,0,0,0,0,0.5e-3,1e-3,1e-3])
a,r = solveq(K,f,bcPrescr,bcVal)

# ----- Compute element flux vector -----

Ed = extractEldisp(Edof,a)
Es = zeros((8,2))
for i in range(8):
    Es[i],Et = flw2qs(Ex[i],Ey[i],ep,D,Ed[i])

# ----- Draw flux vectors and contourlines -----

#sfac=scalfact2(Ex,Ey,Es,0.5);
#eldraw2(Ex,Ey,[1,3,0]); 
#elflux2(Ex,Ey,Es,[1,4],sfac); 
#pltscalb2(sfac,[2e-2 0.06 0.01],4);
#disp('PRESS ENTER TO CONTINUE'); pause; clf;
#eldraw2(Ex,Ey,[1,3,0]); 
#eliso2(Ex,Ey,Ed,5,[1,4]);
#hold off; 
#echo off;

# ----------------- End --------------------------------

