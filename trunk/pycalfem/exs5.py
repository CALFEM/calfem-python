# example exs5
# ----------------------------------------------------------------
# PURPOSE 
#    Analysis of a simply supported beam.
# ----------------------------------------------------------------

# REFERENCES
#     G"oran Sandberg 94-03-08 
#     Karl-Gunnar Olsson 95-09-28
#     Ola Dahlblom 2004-09-21
# ----------------------------------------------------------------

from numpy import *
from pycalfem import *

# ----- Topology -------------------------------------------------

Edof = array([
    [1, 2, 3, 4, 5, 6],
    [4, 5, 6, 7, 8, 9],
    [7, 8, 9, 10, 11, 12]
])

# ----- Stiffness matrix K and load vector f ---------------------

K=mat(zeros((12,12)))
f=mat(zeros((12,1)))
f[4]=-10000.

# ----- Element stiffness matrices  ------------------------------

E=2.1e11
A=45.3e-4
I=2510e-8
ep=array([E,A,I])
ex=array([0.,3.])
ey=array([0.,0.])

Ke=beam2e(ex,ey,ep)

print(Ke)

# ----- Assemble Ke into K ---------------------------------------

K=assem(Edof,K,Ke);

# ----- Solve the system of equations and compute support forces -

bc=array([1,2,11])
(a,r)=solveq(K,f,bc);

# ----- Section forces -------------------------------------------

Ed=extract(Edof,a);

es1,ed1,ec1=beam2s(ex,ey,ep,Ed[0,:],np=10)
es2,ed2,ec2=beam2s(ex,ey,ep,Ed[1,:],np=10)
es3,ed3,ec3=beam2s(ex,ey,ep,Ed[2,:],np=10)

# ----- Results --------------------------------------------------

print("a=")
print(a)
print("r=")
print(r)
print("es1=")
print(es1)
print("es2=")
print(es2)
print("es3=")
print(es3)

print("ed1=")
print(ed1)
print("ed2=")
print(ed2)
print("ed3=")
print(ed3)