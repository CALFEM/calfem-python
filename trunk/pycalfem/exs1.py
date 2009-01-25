from numpy import *
from pycalfem import *

# ----- Topology matrix Edof -------------------------------------

# Edof=[1 1 2;
#       2 2 3;
#       3 2 3];
 
Edof = mat([
    [1,2],
    [2,3],
    [2,3]
])

print Edof

# ----- Stiffness matrix K and load vector f ---------------------

K=mat(zeros((3,3)))
f=mat(zeros((3,1)))
f[1]=100

print K
print f

# ----- Element stiffness matrices  ------------------------------

k=1500.
ep1=k
ep2=2.*k
Ke1=spring1e(ep1)
Ke2=spring1e(ep2)

assem(Edof[0,:], K, Ke2)
assem(Edof[1,:], K, Ke1)
assem(Edof[2,:], K, Ke2)

print K

bc = array([1,3])
a,Q = solveq(K, f, bc)

print a
print Q

ed1=extract(Edof[0,:],a)
ed2=extract(Edof[1,:],a)
ed3=extract(Edof[2,:],a)

ed = extract(Edof,a)

print ed

print ed1
print ed2
print ed3


