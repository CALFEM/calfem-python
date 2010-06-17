# example exs1 
# ----------------------------------------------------------------
# PURPOSE 
#     Linear elastic spring analysis. Introduction to the basic 
#     steps in the finite element method.
# ----------------------------------------------------------------

# REFERENCES
#     P-E Austrell 1994-03-08 
#     K-G Olsson 1995-09-28
#     O Dahlblom 2004-09-06
#     J Lindemann 2009-01-25
# ----------------------------------------------------------------

from numpy import *
from pycalfem import *

# ----- Topology matrix Edof -------------------------------------

Edof = array([
    [1,2],
    [2,3],
    [2,3]
])

# ----- Stiffness matrix K and load vector f ---------------------

K=matrix(zeros((3,3)))
f=matrix(zeros((3,1)))

# ----- Element stiffness matrices  ------------------------------

k=1500.
ep1=k
ep2=2.*k
Ke1=spring1e(ep1)
Ke2=spring1e(ep2)

# ----- Assemble Ke into K ---------------------------------------

assem(Edof[0,:], K, Ke2)
assem(Edof[1,:], K, Ke1)
assem(Edof[2,:], K, Ke2)

print("Stiffness matrix K:")
print(K)

# ----- Solve the system of equations ----------------------------

bc = array([1,3])
f[1]=100

a, r = solveq(K, f, bc)

print("Displacements a:")
print(a)

print("Reaction forces Q:")
print(r)

# ----- Element forces -------------------------------------------

ed1=extract(Edof[0,:],a)
ed2=extract(Edof[1,:],a)
ed3=extract(Edof[2,:],a)

es1=spring1s(ep2,ed1)
es2=spring1s(ep1,ed2)
es3=spring1s(ep2,ed3)

print("N1 = "+str(es1))
print("N2 = "+str(es2))
print("N3 = "+str(es3))


