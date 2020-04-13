# example exs6 
#----------------------------------------------------------------
# PURPOSE 
#    Analysis of a plane frame.
#----------------------------------------------------------------

# REFERENCES
#     G"oran Sandberg 94-03-08 
#     Karl-Gunnar Olsson 95-09-28
#     Anders Olsson 99-03-01
#     Ola Dahlblom 2004-09-14
#----------------------------------------------------------------

import numpy as np
import calfem.core as cfc
import calfem.utils as cfu
import calfem.vis_mpl as cfv

# ----- Topology -------------------------------------------------

Edof = np.array([
    [4,  5,  6, 1,  2,  3],
    [7,  8,  9, 10, 11, 12],
    [4,  5,  6,  7,  8,  9]      
])

# ----- Stiffness matrix K and load vector f ---------------------

K = np.matrix(np.zeros((12,12)))
f = np.matrix(np.zeros((12,1)))
f[3] = 2e+3

# ----- Element stiffness and element load matrices  -------------

E = 200e9
A1 = 2e-3
A2 = 6e-3
I1 = 1.6e-5
I2 = 5.4e-5

ep1 = np.array([E, A1, I1])
ep3 = np.array([E, A2, I2])
ex1 = np.array([0, 0])
ex2 = np.array([6, 6])
ex3 = np.array([0, 6])
ey1 = np.array([4, 0])
ey2 = np.array([4, 0])
ey3 = np.array([4, 4])
eq1 = np.array([0, 0])
eq2 = np.array([0, 0])
eq3 = np.array([0, -10e+3])

Ke1 = cfc.beam2e(ex1, ey1, ep1)
Ke2 = cfc.beam2e(ex2, ey2, ep1)
Ke3, fe3 = cfc.beam2e(ex3, ey3, ep3, eq3)

# ----- Assemble Ke into K ---------------------------------------

cfc.assem(Edof[0,:], K, Ke1);
cfc.assem(Edof[1,:], K, Ke2);
cfc.assem(Edof[2,:], K, Ke3, f, fe3);

# ----- Solve the system of equations and compute reactions ------

bc = np.array([1,2,3,10,11])
a, r = cfc.solveq(K,f,bc)

print("a = ")
print(a)
print("r = ")
print(r)

# ----- Section forces -------------------------------------------

Ed = cfc.extractEldisp(Edof,a);

es1, ed1, ec1 = cfc.beam2s(ex1, ey1, ep1, Ed[0,:], eq1, nep=21)
es2, ed2, ec2 = cfc.beam2s(ex2, ey2, ep1, Ed[1,:], eq2, nep=21)
es3, ed3, ec3 = cfc.beam2s(ex3, ey3, ep3, Ed[2,:], eq3, nep=21)

print("es1 = ")
print(es1)
print("es2 = ")
print(es2)
print("es3 = ")
print(es3)

# ----- Draw deformed frame ---------------------------------------

ex = np.array([
    ex1, ex2, ex3
])
print(ex)

ey = np.array([
    ey1, ey2, ey3
])
print(ey)

cfv.eldraw2(ex, ey)
#cfv.eldisp2(ex, ey, Ed)
cfv.showAndWait()

#figure(1)
#plotpar=[2 1 0];
#eldraw2(ex1,ey1,plotpar);
#eldraw2(ex2,ey2,plotpar);
#eldraw2(ex3,ey3,plotpar);
#sfac=scalfact2(ex3,ey3,Ed(3,:),0.1);
#plotpar=[1 2 1];
#eldisp2(ex1,ey1,Ed(1,:),plotpar,sfac);
#eldisp2(ex2,ey2,Ed(2,:),plotpar,sfac);
#eldisp2(ex3,ey3,Ed(3,:),plotpar,sfac);
#axis([-1.5 7.5 -0.5 5.5]); 
#pltscalb2(sfac,[1e-2 0.5 0]);
#axis([-1.5 7.5 -0.5 5.5]);
#title('displacements')
 
#----- Draw normal force diagram --------------------------------
 
#figure(2)
#plotpar=[2 1];
#sfac=scalfact2(ex1,ey1,es1(:,1),0.2);
#eldia2(ex1,ey1,es1(:,1),plotpar,sfac);
#eldia2(ex2,ey2,es2(:,1),plotpar,sfac);
#eldia2(ex3,ey3,es3(:,1),plotpar,sfac);
#axis([-1.5 7.5 -0.5 5.5]);
#pltscalb2(sfac,[3e4 1.5 0]);
#title('normal force')

#----- Draw shear force diagram ---------------------------------
 
#figure(3)
#plotpar=[2 1];
#sfac=scalfact2(ex3,ey3,es3(:,2),0.2);
#eldia2(ex1,ey1,es1(:,2),plotpar,sfac);
#eldia2(ex2,ey2,es2(:,2),plotpar,sfac);
#eldia2(ex3,ey3,es3(:,2),plotpar,sfac);
#axis([-1.5 7.5 -0.5 5.5]);
#pltscalb2(sfac,[3e4 0.5 0]);
#title('shear force') 

#----- Draw moment diagram --------------------------------------
 
#figure(4)
#plotpar=[2 1];
#sfac=scalfact2(ex3,ey3,es3(:,3),0.2);
#eldia2(ex1,ey1,es1(:,3),plotpar,sfac);
#eldia2(ex2,ey2,es2(:,3),plotpar,sfac);
#eldia2(ex3,ey3,es3(:,3),plotpar,sfac);
#axis([-1.5 7.5 -0.5 5.5]);
#pltscalb2(sfac,[3e4 0.5 0]);
#title('moment') 

