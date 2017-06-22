# example exs3 
# ----------------------------------------------------------------
# PURPOSE 
#    Analysis of a plane truss.
# ----------------------------------------------------------------

# REFERENCES
#     Ola Dahlblom 2004-09-07
#     Jonas Lindemann 2009-01-25
# ----------------------------------------------------------------

import numpy as np
import calfem.core as cfc

# ----- Topology matrix Edof -------------------------------------

Edof = np.array([
    [1,2,5,6],
    [5,6,7,8],
    [3,4,5,6]
])

# ----- Stiffness matrix K and load vector f ---------------------

K = np.matrix(np.zeros((8,8)))
f = np.matrix(np.zeros((8,1)))

# ----- Element properties ---------------------------------------
 
E = 2.0e11
A1 = 6.0e-4
A2 = 3.0e-4
A3 = 10.0e-4
ep1 = [E,A1]
ep2 = [E,A2]
ep3 = [E,A3]
 
#----- Element coordinates --------------------------------------

ex1 = np.array([0., 1.6])
ex2 = np.array([1.6, 1.6])
ex3 = np.array([0., 1.6])

ey1 = np.array([0., 0.])
ey2 = np.array([0., 1.2])
ey3 = np.array([1.2, 0.])
 
#----- Element stiffness matrices  ------------------------------

Ke1 = cfc.bar2e(ex1,ey1,ep1)	 
Ke2 = cfc.bar2e(ex2,ey2,ep2)
Ke3 = cfc.bar2e(ex3,ey3,ep3)	
 
#----- Assemble Ke into K ---------------------------------------

cfc.assem(Edof[0,:],K,Ke1)
cfc.assem(Edof[1,:],K,Ke2)
cfc.assem(Edof[2,:],K,Ke3)

print("Stiffness matrix K:")
print(K)
 
#----- Solve the system of equations ----------------------------

bc = np.array([1,2,3,4,7,8])
f[5] = -80e3
a, r = cfc.solveq(K,f,bc)

print("Displacements a:")
print(a)

print("Reaction forces r:")
print(r)

#----- Element forces -------------------------------------------

ed1 = cfc.extractEldisp(Edof[0,:],a);
N1 = cfc.bar2s(ex1,ey1,ep1,ed1)
ed2 = cfc.extractEldisp(Edof[1,:],a);
N2 = cfc.bar2s(ex2,ey2,ep2,ed2)
ed3 = cfc.extractEldisp(Edof[2,:],a);	
N3 = cfc.bar2s(ex3,ey3,ep3,ed3)

print("N1 = "+str(N1))
print("N2 = "+str(N2))
print("N3 = "+str(N3))
 
