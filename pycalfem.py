from numpy import *

def spring1e(ep):
    """Ke=spring1e(ep)
    -------------------------------------------------------------
     PURPOSE
      Compute element stiffness matrix for spring element.
    
     INPUT:  ep = [k];       spring stiffness or analog quantity
    
     OUTPUT: Ke :            stiffness matrix, dim(Ke)= 2 x 2
    -------------------------------------------------------------
    
     LAST MODIFIED: P-E Austrell 1994-11-02
     Copyright (c)  Division of Structural Mechanics and
                    Department of Solid Mechanics.
                    Lund Institute of Technology
    -------------------------------------------------------------"""
    k = ep
    Ke = mat([[k,-k],[-k,k]],'d')
    return Ke


def assem(edof,K,Ke,f=None,fe=None):
    edofArray = asarray(edof)
    for row in edofArray:
        idx = row-1
        K[ix_(idx,idx)] = K[ix_(idx,idx)] + Ke
        if (f!=None) and (fe!=None):
            f[ix_(idx)] = f[ix_(idx)] + fe
            
def solveq(K,f,bcPrescr,bcVal=None):
    nDofs = K.shape[0]
    bc = ones(nDofs, 'bool')    
    bcDofs = arange(nDofs)
    
    bc[ix_(bcPrescr-1)] = False
    bcDofs = bcDofs[bc]
    
    fsys = f[bcDofs]
    asys = linalg.solve(K[ix_((bcDofs),(bcDofs))], fsys);
    a = zeros([nDofs,1])
    return a

def hooke(ptype,E,v):

    # D=hooke(ptype,E,v)
    #-------------------------------------------------------------
    #  PURPOSE
    #   Calculate the material matrix for a linear
    #   elastic and isotropic material.
    #
    # INPUT:  ptype=1:  plane stress
    #               2:  plane strain
    #               3:  axisymmetry
    #               4:  three dimensional
    #
    #          E : Young's modulus
    #          v : Poissons const.
    #
    # OUTPUT: D : material matrix
    #-------------------------------------------------------------
    
    # LAST MODIFIED: M Ristinmaa 1995-10-25
    # Copyright (c)  Division of Structural Mechanics and
    #                Department of Solid Mechanics.
    #                Lund Institute of Technology
    #-------------------------------------------------------------
    
    if ptype == 1:
        D = E*matrix(
            [[1, v, 0],
             [v, 1, 0],
             [0, 0, (1-v)/2]]
            )/(1-v**2);
    elif ptype == 2:
        D = E/(1+v)*matrix(
            [[1-v, v, v, 0],
             [v, 1-v, v, 0],
             [v, v, 1-v, 0],
             [0, 0, 0, (1-2*v)/2]]
            )/(1-2*v)
    elif ptype == 3:
        D = E/(1+v)*matrix(
            [[1-v, v, v, 0],
             [v, 1-v, v, 0],
             [v, v, 1-v, 0],
             [0, 0, 0, (1-2*v)/2]]
            )/(1-2*v)
    elif ptype == 4:
        D = E*matrix(
            [[1-v, v, v, 0, 0, 0],
             [v, 1-v, v, 0, 0, 0],
             [v, v, 1-v, 0, 0, 0],
             [0, 0, 0, (1-2*v)/2, 0, 0],
             [0, 0, 0, 0, (1-2*v)/2, 0],
             [0, 0, 0, 0, 0, (1-2*v)/2]]
            )/(1+v)/(1-2*v)
    else:
        print "ptype not supported."
        
    return D

def plante(ex,ey,ep,D,eq=None):
    """
     Ke=plante(ex,ey,ep,D)
     [Ke,fe]=plante(ex,ey,ep,D,eq)
    -------------------------------------------------------------
     PURPOSE
      Calculate the stiffness matrix for a triangular plane stress
      or plane strain element.
    
     INPUT:  ex = [x1 x2 x3]         element coordinates
             ey = [y1 y2 y3]
     
             ep = [ptype t ]         ptype: analysis type
                                     t: thickness
     
             D                       constitutive matrix
    
             eq = [bx;               bx: body force x-dir
                   by]               by: body force y-dir
    
     OUTPUT: Ke : element stiffness matrix (6 x 6)
             fe : equivalent nodal forces (6 x 1)
    -------------------------------------------------------------
    
     LAST MODIFIED: M Ristinmaa 1995-10-25
     Copyright (c)  Division of Structural Mechanics and
                    Department of Solid Mechanics.
                    Lund Institute of Technology
    -------------------------------------------------------------
    """

    ptype = ep[0];
    t = ep[1];
    
    bx = 0.0
    by = 0.0
    
    if eq != None:
        bx = eq[0]
        by = eq[1]
        
    C = matrix(
        [[1, ex[0], ey[0], 0, 0, 0 ], 
         [0, 0, 0, 1, ex[0], ey[0] ],
         [1, ex[1], ey[1], 0, 0, 0 ],
         [0, 0, 0, 1, ex[1], ey[1] ],
         [1, ex[2], ey[2], 0, 0, 0 ],
         [0, 0, 0, 1, ex[2], ey[2] ]]
        )
    
    #print "C =", C
    
    A = 0.5*det(matrix([[1, ex[0], ey[0]],[1, ex[1], ey[1]],[1, ex[2], ey[2]]]))
    
    #print "A = ", A

    # --------- plane stress --------------------------------------
    
    if ptype == 1:
        B = matrix([[0,1,0,0,0,0],[0,0,0,0,0,1],[0,0,1,0,1,0]])*inv(C)
        
        #print "B = ", B
        
        colD = D.shape[1]
        
        if colD>3:
            Cm = inv(D)
            Dm = inv(Cm(ix_((0,1,3),(0,1,3))))
        else:
            Dm = D
            
        Ke = B.T*Dm*B*A*t
        fe = A/3*matrix([[bx,by,bx,by,bx,by]]).T*t
        
    return Ke, fe
       
#%--------- plane strain --------------------------------------       
#elseif ptype==2
#       B=[0 1 0 0 0 0
#          0 0 0 0 0 1
#          0 0 1 0 1 0]*inv(C);
#
#       colD=size(D,2);
#       if colD>3
#         Dm=D([1 2 4],[1 2 4]);
#       else
#         Dm=D;
#       end
#
#       Ke=B'*Dm*B*A*t;
#       fe=A/3*[bx by bx by bx by]'*t;
#       
#else
#   error('Error ! Check first argument, ptype=1 or 2 allowed')
#   return
#end
#%--------------------------end--------------------------------
#
#