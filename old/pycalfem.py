# -*- coding: iso-8859-15 -*-

from numpy import *

def spring1e(ep):
    """
    Compute element stiffness matrix for spring element.
    
    Parameters:
    
        ep = k          spring stiffness or analog quantity
        
    Returns:
    
        Ke              stiffness matrix, dim(Ke)= 2 x 2
        
    """
    k = ep
    return mat([[k,-k],[-k,k]],'d')

def spring1s(ep,ed):
    """
    Compute element force in spring element (spring1e).
    
    Parameters:
    
        ep = k          spring stiffness or analog quantity
        ed = [u1 u2]    element displacements
                        u1, u2: nodal displacements
                        
    Returns:
    
        es              element force [N]
    """
    k = ep
    return k*(ed[1]-ed[0]);   

def bar1e(ep):
    """
    Compute element stiffness matrix for spring element.
    
    Parameters:
    
        ep = k          spring stiffness or analog quantity
        
    Returns:
    
        Ke              stiffness matrix, dim(Ke)= 2 x 2
    
    """
    k = ep
    return mat([[k,-k],[-k,k]],'d')

def bar1s(ep,ed):
    """
    Compute element force in spring element (spring1e).
    
    Parameters:
    
        ep = k          spring stiffness or analog quantity
        ed = [u1 u2]    element displacements
                        u1, u2: nodal displacements
                        
    Returns:
    
        es              element force [N]
    
    """
    k = ep
    return k*(ed[1]-ed[0]);   

def bar2e(ex,ey,ep):
    """
    Compute the element stiffness matrix for two dimensional bar element.
    
    Parameters:
    
        ex = [x1 x2]
        ey = [y1 y2]    element node coordinates
    
        ep = [E A]      E: Young's modulus
                        A: Cross section area
                        
    Returns:
    
        Ke              stiffness matrix, dim(Ke)= 4 x 4
    
    """
    E=ep[0]
    A=ep[1]
    
    b = mat([[ex[1]-ex[0]],[ey[1]-ey[0]]])
    L = asscalar(sqrt(b.T*b))
    
    Kle = mat([[1.,-1.],[-1.,1.]])*E*A/L
    
    n = asarray(b.T/L).reshape(2,)
    
    G = mat([
        [n[0],n[1],0.,0.],
        [0.,0.,n[0],n[1]]
    ])
    
    return G.T*Kle*G

def bar2g(ex,ey,ep,N):
    """
    Compute element stiffness matrix for two dimensional geometric
    nonlinear bar element.
    
    Parameters:
    
        ex = [x1 x2]
        ey = [y1 y2]        element node coordinates
        
        ep = [E A]          E: Young's modulus
                            A: Cross sections area
                            
        N                   normal force
        
    Returns:
    
        Ke                  stiffness matrix, dim(Ke) = 4 x 4
        
    """
    E = ep[0]
    A = ep[1]
    
    b = mat([
        [ex[1]-ex[0]],
        [ey[1]-ey[0]]
    ])
    L = asscalar(sqrt(b.T*b))
    
    n = asarray(b.T/L).reshape(2,)

    G = mat([
        [ n[0], n[1], 0.,   0.  ],
        [-n[1], n[0], 0.,   0.  ],
        [ 0.,   0.,   n[0], n[1]],
        [ 0.,   0.,  -n[1], n[0]]
    ])
    
    Kle = E*A/L*mat([
        [ 1, 0,-1, 0],
        [ 0, 0, 0, 0],
        [-1, 0, 1, 0],
        [ 0, 0, 0, 0]
    ])+N/L*mat([
        [ 0, 0, 0, 0],
        [ 0, 1, 0,-1],
        [ 0, 0, 0, 0],
        [ 0,-1, 0, 1]
    ])

    return G.T*Kle*G

def bar2s(ex,ey,ep,ed):
    """
    Compute normal force in two dimensional bar element.
    
    Parameters:
    
        ex = [x1 x2]
        ey = [y1 y2]        element coordinates
    
        ep = [E A]          E : Young's modulus
                            A : Cross section area
    
        ed = [u1 u2 u3 u4]  element displacement vector
        
    Returns:
    
        es                  element force [N]
    
    """
    E=ep[0]
    A=ep[1]
    
    b = mat([[ex[1]-ex[0]],[ey[1]-ey[0]]])
    L = asscalar(sqrt(b.T*b))
    
    Kle = mat([[1.,-1.],[-1.,1.]])*E*A/L
    
    n = asarray(b.T/L).reshape(2,) 
    
    G = mat([
        [n[0],n[1],0.,0.],
        [0.,0.,n[0],n[1]]
    ])
    
    u=asmatrix(ed).T
    N=E*A/L*mat([[-1.,1.]])*G*u
    return asscalar(N)
    
def bar3e(ex,ey,ez,ep):
    """
    Compute element stiffness matrix for three dimensional bar element.
    
    Parameters:
    
        ex = [x1 x2]
        ey = [y1 y2]        element node coordinates
        ez = [z1 z2]
        
        ep = [E A]          E: Young's modulus
                            A: Cross sections area
        
    Returns:
    
        Ke                  stiffness matrix, dim(Ke) = 6 x 6
        
    """
    E = ep[0]
    A = ep[1]
    
    b = mat([
        [ex[1]-ex[0]],
        [ey[1]-ey[0]],
        [ez[1]-ez[0]]
    ])
    L = asscalar(sqrt(b.T*b))
    
    n = asarray(b.T/L).reshape(3)

    G = mat([
        [ n[0], n[1], n[2], 0.,   0.,   0.  ],
        [ 0.,   0.,   0.,   n[0], n[1], n[2]]
    ])
    
    Kle = E*A/L*mat([
        [ 1,-1],
        [-1, 1]
    ])

    return G.T*Kle*G

def bar3s(ex,ey,ez,ep,ed):
    """
    Compute normal force in three dimensional bar element.
    
    Parameters:
    
        ex = [x1 x2]
        ey = [y1 y2]        element node coordinates
        ez = [z1 z2]
        
        ep = [E A]          E: Young's modulus
                            A: Cross sections area
        
        ed = [u1 ... u6]    element displacements
        
    Returns:
    
        es = [N]            normal force

    """
    E = ep[0]
    A = ep[1]
    
    b = mat([
        [ex[1]-ex[0]],
        [ey[1]-ey[0]],
        [ez[1]-ez[0]]
    ])
    L = asscalar(sqrt(b.T*b))
    
    n = asarray(b.T/L).reshape(3)

    G = mat([
        [ n[0], n[1], n[2], 0.  , 0.  , 0.  ],
        [ 0.  , 0.  , 0.  , n[0], n[1], n[2]]
    ])
    
    Kle = E*A/L*mat([
        [ 1,-1],
        [-1, 1]
    ])

    u = asmatrix(ed).T
    N = E*A/L*mat([[-1.,1.]])*G*u

    return asscalar(N)

def beam2e(ex,ey,ep,eq=None):
    """
    Compute the stiffness matrix for a two dimensional beam element.
    
    Parameters:
     
        ex = [x1 x2]
        ey = [y1 y2]        element node coordinates
    
        ep = [E A I]        element properties
                            E: Young's modulus
                            A: Cross section area
                            I: Moment of inertia
    
        eq = [qx qy]        distributed loads, local directions
        
    Returns:
     
        Ke                  element stiffness matrix (6 x 6)
    
        fe                  element load vector (6 x 1)
    
    """

    b=mat([[ex[1]-ex[0]],[ey[1]-ey[0]]])
    L = asscalar(sqrt(b.T*b))
    n = asarray(b.T/L).reshape(2,) 
    
    E=ep[0]
    A=ep[1]
    I=ep[2]
    
    qx=0.
    qy=0.
    if eq!=None:
        qx=eq[0]
        qy=eq[1]
        
    Kle = mat([
        [E*A/L,      0.,          0.,    -E*A/L,    0.,        0.      ],
        [  0.,    12*E*I/L**3., 6*E*I/L**2.,    0., -12*E*I/L**3., 6*E*I/L**2. ],
        [  0.,    6*E*I/L**2.,  4*E*I/L,      0., -6*E*I/L**2.,  2*E*I/L   ],
        [-E*A/L,     0.,          0.,     E*A/L,    0.,        0.      ],
        [  0.,   -12*E*I/L**3.,-6*E*I/L**2.,    0.,  12*E*I/L**3.,-6*E*I/L**2. ],
        [  0.,    6*E*I/L**2.,  2*E*I/L,      0.,  -6*E*I/L**2., 4*E*I/L   ]
    ])
     
    fle=L*mat([qx/2, qy/2, qy*L/12, qx/2, qy/2, -qy*L/12]).T
     
    G=mat([
        [ n[0], n[1],  0.,    0.,    0.,   0.],
        [-n[1], n[0],  0.,    0.,    0.,   0.],
        [0.,    0.,    1.,    0.,    0.,   0.],
        [0.,    0.,    0.,   n[0],  n[1],  0.],
        [0.,    0.,    0.,  -n[1],  n[0],  0.],
        [0.,    0.,    0.,    0.,    0.,   1.]
    ])
    
    Ke=G.T*Kle*G
    fe=G.T*fle
    
    if eq==None:
        return Ke
    else:
        return Ke,fe
    
def beam2s(ex,ey,ep,ed,eq=None,np=None):
    """
    Compute section forces in two dimensional beam element (beam2e).
    
    Parameters:
 
        ex = [x1 x2]
        ey = [y1 y2]        element node coordinates

        ep = [E A I]        element properties,
                            E:  Young's modulus
                            A:  cross section area
                            I:  moment of inertia

        ed = [u1 ... u6]    element displacements

        eq = [qx qy]        distributed loads, local directions 

        n                   number of evaluation points ( default=2 )
        
    Returns:
          
        es = [ N1 V1 M1     section forces, local directions, in 
               N2 V2 M2     n points along the beam, dim(es)= n x 3
               .........]  
           
        edi = [ u1 v1       element displacements, local directions,
                u2 v2       in n points along the beam, dim(es)= n x 2
                .......]    

            eci = [ x1      local x-coordinates of the evaluation 
                    x2      points, (x1=0 and xn=L)
                    ...]
    
    """
    EA=ep[0]*ep[1]
    EI=ep[0]*ep[2]
    b=mat([
        [ex[1]-ex[0]],
        [ey[1]-ey[0]]
    ])
    
    L = asscalar(sqrt(b.T*b))
    n = asarray(b.T/L).reshape(2,)
    
    qx=0.
    qy=0.
    
    if eq!=None:
        qx=eq[0]
        qy=eq[1] 
      
    ne=2
    
    if np!=None:
        ne = np
        
    C=mat([
        [0.,   0.,   0.,    1.,   0.,   0.],
        [0.,   0.,   0.,    0.,   0.,   1.],
        [0.,   0.,   0.,    0.,   1.,   0.],
        [L,   0.,   0.,    1.,   0.,   0.],
        [0.,   L**3, L**2,   0.,   L,    1.],
        [0., 3*L**2, 2*L,   0.,   1.,   0.]
    ])
   
    G=mat([
        [ n[0], n[1],  0.,    0.,    0.,   0.],
        [-n[1], n[0],  0.,    0.,    0.,   0.],
        [0.,    0.,    1.,    0.,    0.,   0.],
        [0.,    0.,    0.,   n[0],  n[1],  0.],
        [0.,    0.,    0.,  -n[1],  n[0],  0.],
        [0.,    0.,    0.,    0.,    0.,   1.]
    ])
    
    M=ravel(C.I*(G*asmatrix(ed).T-matrix([0., 0., 0., -qx*L**2/(2*EA), qy*L**4/(24*EI), qy*L**3/(6*EI)]).T))
    A=matrix([M[0],M[3]]).T
    B=matrix([M[1],M[2],M[4],M[5]]).T
    
    x=asmatrix(arange(0.,L+L/(ne-1),L/(ne-1))).T
    zero=asmatrix(zeros([len(x)])).T
    one=asmatrix(ones([len(x)])).T
    
    u=concatenate((x,one),1)*A-power(x,2)*qx/(2*EA)
    du=concatenate((one,zero),1)*A-x*qx/EA
    v=concatenate((power(x,3),power(x,2),x,one),1)*B+power(x,4)*qy/(24*EI)
    d2v=concatenate((6*x,2*one,zero,zero),1)*B+power(x,2)*qy/(2*EI)
    d3v=concatenate((6*one,zero,zero,zero),1)*B+x*qy/EI
    
    N=EA*du
    M=EI*d2v
    V=-EI*d3v
    edi=concatenate((u,v),1)
    eci=x
    es=concatenate((N,V,M),1)
    
    return (es,edi,eci)

def beam2t(ex,ey,ep,eq=None):
    """
    Compute the stiffness matrix for a two dimensional elastic
    Timoshenko beam element.
    
    Parameters:
     
        ex = [x1 x2]
        ey = [y1 y2]        element node coordinates
    
        ep = [E G A I ks]   element properties
                              E: Young's modulus
                              G: Shear modulus
                              A: Cross section area
                              I: Moment of inertia
                             ks: Shear correction factor
    
        eq = [qx qy]        distributed loads, local directions
        
    Returns:
     
        Ke                  element stiffness matrix (6 x 6)
    
        fe                  element load vector (6 x 1)
    
    """

    b = mat([[ex[1]-ex[0]],[ey[1]-ey[0]]])
    L = asscalar(sqrt(b.T*b))
    n = asarray(b.T/L).reshape(2)
    
    E = ep[0]
    Gm = ep[1]
    A = ep[2]
    I = ep[3]
    ks = ep[4]
        
    qx = 0.
    qy = 0.
    if eq != None:
        qx = eq[0]
        qy = eq[1]
    
    m = (12/L**2)*(E*I/(Gm*A*ks))
    
    Kle = E/(1+m)*mat([
        [A*(1+m)/L,      0.,         0.,        -A*(1+m)/L,     0.,          0.      ],
        [0.,         12*I/L**3., 6*I/L**2.,         0.,    -12*I/L**3., 6*I/L**2.    ],
        [0.,         6*I/L**2.,  4*I*(1+m/4.)/L,    0.,    -6*I/L**2.,  2*I*(1-m/2)/L],
        [-A*(1+m)/L,     0.,         0.,         A*(1+m)/L,     0.,          0.      ],
        [0.,        -12*I/L**3.,-6*I/L**2.,         0.,     12*I/L**3.,-6*I/L**2.    ],
        [0.,         6*I/L**2.,  2*I*(1-m/2)/L,     0.,    -6*I/L**2.,  4*I*(1+m/4)/L]
    ])

    fle = L*mat([qx/2, qy/2, qy*L/12, qx/2, qy/2, -qy*L/12]).T
    
    G = mat([
        [ n[0], n[1],  0.,   0.,   0.,   0.],
        [-n[1], n[0],  0.,   0.,   0.,   0.],
        [  0.,   0.,   1.,   0.,   0.,   0.],
        [  0.,   0.,   0.,  n[0], n[1],  0.],
        [  0.,   0.,   0., -n[1], n[0],  0.],
        [  0.,   0.,   0.,   0.,   0.,   1.]
    ])
    
    Ke = G.T*Kle*G
    fe = G.T*fle
    
    if eq == None:
        return Ke
    else:
        return Ke,fe

def beam2ts(ex,ey,ep,ed,eq=None,np=None):
    """
    Compute section forces in two dimensional beam element (beam2e).
    
    Parameters:
 
        ex = [x1, x2]
        ey = [y1, y2]       element node coordinates

        ep = [E,G,A,I,ks]   element properties,
                              E:  Young's modulus
                              G:  shear modulus
                              A:  cross section area
                              I:  moment of inertia

        ed = [u1, ... ,u6]  element displacements

        eq = [qx, qy]       distributed loads, local directions 

        n                   number of evaluation points ( default=2 )
        
    Returns:
          
        es = [[N1,V1,M1],   section forces, local directions, in 
              [N2,V2,M2],   n points along the beam, dim(es)= n x 3
              ..........]  
    
        edi = [[u1,v1,teta1],   element displacements, local directions,
               [u2,v2,teta2],   and rotation of cross section at
               .............]   in n points along the beam, dim(es)= n x 2
    
    (Note! Rotation of the cross section is not equal to dv/dx for Timoshenko beam element)
    
        eci = [[x1],    local x-coordinates of the evaluation 
               [x2],    points, (x1=0 and xn=L)
               ....]
    
    """
    EA = ep[0]*ep[2]
    EI = ep[0]*ep[3]
    GAK = ep[1]*ep[2]*ep[4]
    alfa = EI/GAK
    
    b = mat([
        [ex[1]-ex[0]],
        [ey[1]-ey[0]]
    ])
    L = asscalar(sqrt(b.T*b))
    n = asarray(b.T/L).reshape(2)
    
    qx = 0.
    qy = 0.
    if eq != None:
        qx = eq[0]
        qy = eq[1] 
      
    ne = 2
    
    if np != None:
        ne = np
        
    C = mat([
        [ 0., 0.,              0.,   1., 0., 0.],
        [ 0., 0.,              0.,   0., 0., 1.],
        [ 0., 6*alfa,          0.,   0., 1., 0.],
        [ L,  0.,              0.,   1., 0., 0.],
        [ 0., L**3,            L**2, 0., L,  1.],
        [ 0., 3*(L**2+2*alfa), 2*L,  0., 1., 0.]
    ])
   
    G = mat([
        [ n[0], n[1], 0., 0.,   0.,   0.],
        [-n[1], n[0], 0., 0.,   0.,   0.],
        [ 0.,   0.,   1., 0.,   0.,   0.],
        [ 0.,   0.,   0., n[0], n[1], 0.],
        [ 0.,   0.,   0.,-n[1], n[0], 0.],
        [ 0.,   0.,   0., 0.,   0.,   1.]
    ])
    
    M = ravel(C.I*(G*asmatrix(ed).T-mat([0., 0., 0., -qx*L**2/(2*EA), qy*L**4/(24*EI)-qy*L**2/(2*GAK), qy*L**3/(6*EI)]).T))
    C2 = mat([M[0], M[3]]).T
    C4 = mat([M[1], M[2], M[4], M[5]]).T
    
    x = asmatrix(arange(0., L+L/(ne-1), L/(ne-1))).T
    zero = asmatrix(zeros([len(x)])).T
    one = asmatrix(ones([len(x)])).T
    
    u = concatenate((x,one),1)*C2-qx/(2*EA)*power(x,2)
    du = concatenate((one,zero),1)*C2-qx*x/EA
    
    v = concatenate((power(x,3),power(x,2),x,one),1)*C4+qy/(24*EI)*power(x,4)-qy/(2*GAK)*power(x,2)
    dv = concatenate((3*power(x,2),2*x,one,zero),1)*C4+qy*power(x,3)/(6*EI)-qy*x/GAK
    
    teta = concatenate((3*(power(x,2)+2*alfa*one),2*x,one,zero),1)*C4+qy*power(x,3)/(6*EI)
    dteta = concatenate((6*x,2*one,zero,zero),1)*C4+qy*power(x,2)/(2*EI)
    
    N = EA*du
    M = EI*dteta
    V = GAK*(dv-teta)
    
    es = concatenate((N,V,M),1)
    edi = concatenate((u,v,teta),1)
    eci = x

    if np != None:
        return es,edi,eci
    else:
        return es

def beam2w(ex,ey,ep,eq=None):
    """
    Compute the stiffness matrix for a two dimensional beam element
    on elastic foundation.
    
    Parameters:
 
        ex = [x1, x2]
        ey = [y1, y2]       element node coordinates

        ep = [E,A,I,ka,kt]  element properties,
                              E:  Young's modulus
                              A:  cross section area
                              I:  moment of inertia
                             ka:  axial foundation stiffness
                             kt:  transversal foundation stiffness

        eq = [qx, qy]       distributed loads, local directions

    Returns:

        Ke                  beam stiffness matrix (6 x 6)
        
        fe                  element load vector (6 x 1)
    """
    b = mat([[ex[1]-ex[0]],[ey[1]-ey[0]]])
    L = asscalar(sqrt(b.T*b))
    n = asarray(b/L).reshape(2)
    
    E,A,I,ka,kt = ep
    
    qx = 0
    qy = 0
    if eq != None:
        qx,qy = eq
    
    K1 = mat([
        [ E*A/L,  0,           0,         -E*A/L, 0,           0         ],
        [ 0,      12*E*I/L**3, 6*E*I/L**2, 0,    -12*E*I/L**3, 6*E*I/L**2],
        [ 0,      6*E*I/L**2,  4*E*I/L,    0,    -6*E*I/L**2,  2*E*I/L   ],
        [-E*A/L,  0,           0,          E*A/L, 0,           0         ],
        [ 0,     -12*E*I/L**3,-6*E*I/L**2, 0,     12*E*I/L**3,-6*E*I/L**2],
        [ 0,      6*E*I/L**2,  2*E*I/L,    0,    -6*E*I/L**2,  4*E*I/L   ]
        ])
    
    K2 = L/420*mat([
        [ 140*ka, 0,       0,         70*ka,  0,       0        ],
        [ 0,      156*kt,  22*kt*L,   0,      54*kt,  -13*kt*L  ],
        [ 0,      22*kt*L, 4*kt*L**2, 0,      13*kt*L,-3*kt*L**2],
        [ 70*ka,  0,       0,         140*ka, 0,       0        ],
        [ 0,      54*kt,   13*kt*L,   0,      156*kt, -22*kt*L  ],
        [ 0,     -13*kt*L,-3*kt*L**2, 0,     -22*kt*L, 4*kt*L**2]
    ])
    
    Kle = K1+K2
    fle = L*mat([qx/2, qy/2, qy*L/12, qx/2, qy/2, -qy*L/12]).T

    G = mat([
        [ n[0], n[1], 0, 0,    0,    0],
        [-n[1], n[0], 0, 0,    0,    0],
        [ 0,    0,    1, 0,    0,    0],
        [ 0,    0,    0, n[0], n[1], 0],
        [ 0,    0,    0,-n[1], n[0], 0],
        [ 0,    0,    0, 0,    0,    1]
    ])
    
    Ke = G.T*Kle*G
    fe = G.T*fle
    
    if eq != None:
        return Ke,fe
    else:
        return Ke
    
def beam2ws(ex,ey,ep,ed,eq=None):
    """
    Compute section forces in a two dimensional beam element
    on elastic foundation.
    
    Parameters:
 
        ex = [x1, x2]
        ey = [y1, y2]           element node coordinates

        ep = [E,A,I,ka,kt]      element properties,
                                  E:  Young's modulus
                                  A:  cross section area
                                  I:  moment of inertia
                                 ka:  axial foundation stiffness
                                 kt:  transversal foundation stiffness

        ed = [u1, ... ,u6]      element displacement vector

        eq = [qx, qy]           distributed loads, local directions

    Returns:

        es = [[N1, V1, M1],
              [N2, V2, M2]]     element forces, local direction
    """
    if asmatrix(ed).shape[0] > 1:
        print("Only one row is allowed in the ed matrix !!!")
        return

    b = mat([
        [ex[1]-ex[0]],
        [ey[1]-ey[0]]
    ])
    L = asscalar(sqrt(b.T*b))
    n = asarray(b/L).reshape(2,)
    
    E,A,I,ka,kt = ep
    
    qx = 0
    qy = 0
    if eq != None:
        qx,qy = eq
    
    K1 = mat([
        [ E*A/L, 0,           0,         -E*A/L, 0,           0         ],
        [ 0,     12*E*I/L**3, 6*E*I/L**2, 0,    -12*E*I/L**3, 6*E*I/L**2],
        [ 0,     6*E*I/L**2,  4*E*I/L,    0,    -6*E*I/L**2,  2*E*I/L   ],
        [-E*A/L, 0,           0,          E*A/L, 0,           0         ],
        [ 0,    -12*E*I/L**3,-6*E*I/L**2, 0,     12*E*I/L**3,-6*E*I/L**2],
        [ 0,     6*E*I/L**2,  2*E*I/L,    0,    -6*E*I/L**2,  4*E*I/L   ]
        ])
    
    K2 = L/420*mat([
        [ 140*ka, 0,       0,         70*ka,  0,       0        ],
        [ 0,      156*kt,  22*kt*L,   0,      54*kt,  -13*kt*L  ],
        [ 0,      22*kt*L, 4*kt*L**2, 0,      13*kt*L,-3*kt*L**2],
        [ 70*ka,  0,       0,         140*ka, 0,       0        ],
        [ 0,      54*kt,   13*kt*L,   0,      156*kt, -22*kt*L  ],
        [ 0,     -13*kt*L,-3*kt*L**2, 0,     -22*kt*L, 4*kt*L**2]
    ])
    
    Kle = K1+K2
    fle = L*mat([qx/2, qy/2, qy*L/12, qx/2, qy/2, -qy*L/12]).T
    
    G = mat([
        [ n[0], n[1], 0, 0,    0,    0],
        [-n[1], n[0], 0, 0,    0,    0],
        [ 0,    0,    1, 0,    0,    0],
        [ 0,    0,    0, n[0], n[1], 0],
        [ 0,    0,    0,-n[1], n[0], 0],
        [ 0,    0,    0, 0,    0,    1]
    ])

    P = Kle*G*asmatrix(ed).T-fle

    es = mat([
        [-P[0,0],-P[1,0],-P[2,0]],
        [ P[3,0], P[4,0], P[5,0]]
    ])
    
    return es

def beam2g(ex,ey,ep,N,eq=None):
    """
    Compute the element stiffness matrix for a two dimensional
    beam element with respect to geometric nonlinearity.
    
    Parameters:
 
        ex = [x1, x2]
        ey = [y1, y2]           element node coordinates

        ep = [E,A,I]            element properties;
                                  E:  Young's modulus
                                  A:  cross section area
                                  I:  moment of inertia

        N                       axial force in the beam

        eq                      distributed transverse load

    Returns:

        Ke                      element stiffness matrix (6 x 6)
        
        fe                      element load vector (6 x 1)
    """
    if eq != None:
        if size(eq) > 1:
            print("eq should be a scalar !!!")
            return
        else:
            q = eq[0]
    else:
        q = 0
    
    b = mat([
        [ex[1]-ex[0]],
        [ey[1]-ey[0]]
    ])
    L = asscalar(sqrt(b.T*b))
    n = asarray(b/L).reshape(2,)
    
    E,A,I = ep
    
    rho = -N*L**2/(pi**2*E*I)
    
    kL = pi*sqrt(abs(rho))+finfo(float).eps

    if rho > 0:
        f1 = (kL/2)/tan(kL/2)
        f2 = (1/12.)*kL**2/(1-f1)
        f3 = f1/4+3*f2/4
        f4 = -f1/2+3*f2/2
        f5 = f1*f2
        h = 6*(2/kL**2-(1+cos(kL))/(kL*sin(kL)))
    elif rho < 0:
        f1 = (kL/2)/tanh(kL/2)
        f2 = -(1/12.)*kL**2/(1-f1)
        f3 = f1/4+3*f2/4
        f4 = -f1/2+3*f2/2
        f5 = f1*f2
        h = -6*(2/kL**2-(1+cosh(kL))/(kL*sinh(kL)))
    else:
        f1 = f2 = f3 = f4 = f5 = h = 1

    Kle = mat([
        [ E*A/L, 0.,              0.,            -E*A/L, 0.,              0.            ],
        [ 0.,    12*E*I*f5/L**3., 6*E*I*f2/L**2., 0.,   -12*E*I*f5/L**3., 6*E*I*f2/L**2.],
        [ 0.,    6*E*I*f2/L**2.,  4*E*I*f3/L,     0.,   -6*E*I*f2/L**2.,  2*E*I*f4/L    ],
        [-E*A/L, 0.,              0.,             E*A/L, 0.,              0.            ],
        [ 0.,   -12*E*I*f5/L**3.,-6*E*I*f2/L**2., 0.,    12*E*I*f5/L**3.,-6*E*I*f2/L**2.],
        [ 0.,    6*E*I*f2/L**2.,  2*E*I*f4/L,     0.,   -6*E*I*f2/L**2.,  4*E*I*f3/L    ]
    ])

    fle = q*L*mat([0.,1/2.,L*h/12,0.,1/2.,-L*h/12]).T
    
    G = mat([
        [ n[0], n[1], 0, 0,    0,    0],
        [-n[1], n[0], 0, 0,    0,    0],
        [ 0,    0,    1, 0,    0,    0],
        [ 0,    0,    0, n[0], n[1], 0],
        [ 0,    0,    0,-n[1], n[0], 0],
        [ 0,    0,    0, 0,    0,    1]
    ])

    Ke = G.T*Kle*G
    fe = G.T*fle
    
    if eq != None:
        return Ke,fe
    else:
        return Ke
    
def beam2gs(ex,ey,ep,ed,N,eq=None):
    """
    Calculate section forces in a two dimensional nonlinear
    beam element.

    Parameters:
 
        ex = [x1, x2]
        ey = [y1, y2]           element node coordinates

        ep = [E,A,I]            element properties;
                                  E:  Young's modulus
                                  A:  cross section area
                                  I:  moment of inertia

        ed = [u1, ... ,u6]      element displacement vector

        N                       axial force

        eq = [qy]               distributed transverse load

    Returns:

        es = [[N1,V1,M1],       element forces, local directions
              [N2,V2,M2]]
    """
    if eq != None:
        eq = eq[0]
    else:
        eq = 0
    
    b = mat([
        [ex[1]-ex[0]],
        [ey[1]-ey[0]]
    ])
    L = asscalar(sqrt(b.T*b))
    n = asarray(b/L).reshape(2,)
    
    E,A,I = ep
    
    rho = -N*L**2/(pi**2*E*I)
    
    eps = 2.2204e-16
    kL = pi*sqrt(abs(rho))+eps

    if rho > 0:
        f1 = (kL/2)/tan(kL/2)
        f2 = (1/12.)*kL**2/(1-f1)
        f3 = f1/4+3*f2/4
        f4 = -f1/2+3*f2/2
        f5 = f1*f2
        h = 6*(2/kL**2-(1+cos(kL))/(kL*sin(kL)))
    elif rho < 0:
        f1 = (kL/2)/tanh(kL/2)
        f2 = -(1/12.)*kL**2/(1-f1)
        f3 = f1/4+3*f2/4
        f4 = -f1/2+3*f2/2
        f5 = f1*f2
        h = -6*(2/kL**2-(1+cosh(kL))/(kL*sinh(kL)))
    else:
        f1 = f2 = f3 = f4 = f5 = h = 1
    
    Kle = mat([
        [ E*A/L, 0,              0,            -E*A/L, 0,              0            ],
        [ 0,     12*E*I*f5/L**3, 6*E*I*f2/L**2, 0,    -12*E*I*f5/L**3, 6*E*I*f2/L**2],
        [ 0,     6*E*I*f2/L**2,  4*E*I*f3/L,    0,    -6*E*I*f2/L**2,  2*E*I*f4/L   ],
        [-E*A/L, 0,              0,             E*A/L, 0,              0            ],
        [ 0,    -12*E*I*f5/L**3,-6*E*I*f2/L**2, 0,     12*E*I*f5/L**3,-6*E*I*f2/L**2],
        [ 0,     6*E*I*f2/L**2,  2*E*I*f4/L,    0,    -6*E*I*f2/L**2,  4*E*I*f3/L   ]
    ])

    fle = eq*L*mat([0,1/2.,L*h/12,0,1/2.,-L*h/12]).T
    
    G = mat([
        [ n[0], n[1], 0, 0,    0,    0],
        [-n[1], n[0], 0, 0,    0,    0],
        [ 0,    0,    1, 0,    0,    0],
        [ 0,    0,    0, n[0], n[1], 0],
        [ 0,    0,    0,-n[1], n[0], 0],
        [ 0,    0,    0, 0,    0,    1]
    ])
    
    u = asmatrix(ed).T
    P = Kle*G*u-fle

    es = mat([
        [-P[0,0],-P[1,0],-P[2,0]],
        [ P[3,0], P[4,0], P[5,0]]
    ])
    
    return es

def beam2d(ex,ey,ep):
    """
    Calculate the stiffness matrix Ke, the mass matrix Me
    and the damping matrix Ce for a 2D elastic Bernoulli
    beam element.

    Parameters:
 
        ex = [x1, x2]
        ey = [y1, y2]           element node coordinates

        ep = [E,A,I,m,(a,b)]    element properties;
                                  E:  Young's modulus
                                  A:  cross section area
                                  I:  moment of inertia
                                  m:  mass per unit length
                                a,b:  damping coefficients,
                                      Ce=aMe+bKe

    Returns:

        Ke                      element stiffness matrix (6 x 6)
        Me                      element mass martix
        Ce                      element damping matrix, optional
    """
    b = mat([
        [ex[1]-ex[0]],
        [ey[1]-ey[0]]
    ])
    L = asscalar(sqrt(b.T*b))
    n = asarray(b/L).reshape(2,)
    
    a = 0
    b = 0
    if size(ep) == 4:
        E,A,I,m = ep
    elif size(ep) == 6:
        E,A,I,m,a,b = ep
    
    Kle = mat([
        [ E*A/L, 0,           0,         -E*A/L, 0,           0         ],
        [ 0,     12*E*I/L**3, 6*E*I/L**2, 0,    -12*E*I/L**3, 6*E*I/L**2],
        [ 0,     6*E*I/L**2,  4*E*I/L,    0,    -6*E*I/L**2,  2*E*I/L   ],
        [-E*A/L, 0,           0,          E*A/L, 0,           0         ],
        [ 0,    -12*E*I/L**3,-6*E*I/L**2, 0,     12*E*I/L**3,-6*E*I/L**2],
        [ 0,     6*E*I/L**2,  2*E*I/L,    0,    -6*E*I/L**2,  4*E*I/L   ]
    ])

    Mle = m*L/420*mat([
        [ 140, 0,    0,      70,  0,    0     ],
        [ 0,   156,  22*L,   0,   54,  -13*L  ],
        [ 0,   22*L, 4*L**2, 0,   13*L,-3*L**2],
        [ 70,  0,    0,      140, 0,    0     ],
        [ 0,   54,   13*L,   0,   156, -22*L  ],
        [ 0,  -13*L,-3*L**2, 0,  -22*L, 4*L**2]
    ])
    
    Cle = a*Mle+b*Kle

    G = mat([
        [ n[0], n[1], 0, 0,    0,    0],
        [-n[1], n[0], 0, 0,    0,    0],
        [ 0,    0,    1, 0,    0,    0],
        [ 0,    0,    0, n[0], n[1], 0],
        [ 0,    0,    0,-n[1], n[0], 0],
        [ 0,    0,    0, 0,    0,    1]
    ])
    
    Ke = G.T*Kle*G
    Me = G.T*Mle*G
    Ce = G.T*Cle*G
    
    if size(ep) == 4:
        return Ke,Me
    elif size(ep) == 6:
        return Ke,Me,Ce

def beam3e(ex,ey,ez,eo,ep,eq=None):
    """
    Calculate the stiffness matrix for a 3D elastic Bernoulli
    beam element.
    
    Parameters:
     
        ex = [x1 x2]
        ey = [y1 y2]
        ez = [z1 z2]            element node coordinates
        
        eo = [xz yz zz]         orientation of local z axis
        
        ep = [E G A Iy Iz Kv]   element properties
                                  E: Young's modulus
                                  G: Shear modulus
                                  A: Cross section area
                                 Iy: Moment of inertia, local y-axis
                                 Iz: Moment of inertia, local z-axis
                                 Kv: Saint-Venant's torsion constant
    
        eq = [qx qy qz qw]      distributed loads

    Returns:

        Ke                      beam stiffness matrix (12 x 12)

        fe                      equivalent nodal forces (12 x 1)

    """
    b = mat([
        [ex[1]-ex[0]],
        [ey[1]-ey[0]],
        [ez[1]-ez[0]]
        ])
    L = asscalar(sqrt(b.T*b))
    n1 = asarray(b.T/L).reshape(3,)
    
    eo = asmatrix(eo)
    lc = asscalar(sqrt(eo*eo.T))
    n3 = asarray(eo/lc).reshape(3,)
    
    E,Gs,A,Iy,Iz,Kv = ep

    qx = 0.
    qy = 0.
    qz = 0.
    qw = 0.
    if eq != None:
        qx,qy,qz,qw = eq

    a = E*A/L
    b = 12*E*Iz/L**3
    c = 6*E*Iz/L**2
    d = 12*E*Iy/L**3
    e = 6*E*Iy/L**2
    f = Gs*Kv/L
    g = 2*E*Iy/L
    h = 2*E*Iz/L

    Kle = mat([
        [ a, 0, 0, 0, 0,   0,  -a, 0, 0, 0, 0,   0  ],
        [ 0, b, 0, 0, 0,   c,   0,-b, 0, 0, 0,   c  ],
        [ 0, 0, d, 0,-e,   0,   0, 0,-d, 0,-e,   0  ],
        [ 0, 0, 0, f, 0,   0,   0, 0, 0,-f, 0,   0  ],
        [ 0, 0,-e, 0, 2*g, 0,   0, 0, e, 0, g,   0  ],
        [ 0, c, 0, 0, 0,   2*h, 0,-c, 0, 0, 0,   h  ],
        [-a, 0, 0, 0, 0,   0,   a, 0, 0, 0, 0,   0  ],
        [ 0,-b, 0, 0, 0,  -c,   0, b, 0, 0, 0,  -c  ],
        [ 0, 0,-d, 0, e,   0,   0, 0, d, 0, e,   0  ],
        [ 0, 0, 0,-f, 0,   0,   0, 0, 0, f, 0,   0  ],
        [ 0, 0,-e, 0, g,   0,   0, 0, e, 0, 2*g, 0  ],
        [ 0, c, 0, 0, 0,   h,   0,-c, 0, 0, 0,   2*h]
    ])

    fle = L/2*mat([qx, qy, qz, qw, -qz*L/6, qy*L/6, qx, qy, qz, qw, qz*L/6, -qy*L/6]).T

    n2 = array([0.,0.,0.])
    n2[0] = n3[1]*n1[2]-n3[2]*n1[1]
    n2[1] = -n1[2]*n3[0]+n1[0]*n3[2]
    n2[2] = n3[0]*n1[1]-n1[0]*n3[1]

    #An = append([n1,n2],[n3],0)

    G = mat([
        [ n1[0], n1[1], n1[2], 0,     0,     0,     0,     0,     0,     0,     0,     0    ],
        [ n2[0], n2[1], n2[2], 0,     0,     0,     0,     0,     0,     0,     0,     0    ],
        [ n3[0], n3[1], n3[2], 0,     0,     0,     0,     0,     0,     0,     0,     0    ],
        [ 0,     0,     0,     n1[0], n1[1], n1[2], 0,     0,     0,     0,     0,     0    ],
        [ 0,     0,     0,     n2[0], n2[1], n2[2], 0,     0,     0,     0,     0,     0    ],
        [ 0,     0,     0,     n3[0], n3[1], n3[2], 0,     0,     0,     0,     0,     0    ],
        [ 0,     0,     0,     0,     0,     0,     n1[0], n1[1], n1[2], 0,     0,     0    ],
        [ 0,     0,     0,     0,     0,     0,     n2[0], n2[1], n2[2], 0,     0,     0    ],
        [ 0,     0,     0,     0,     0,     0,     n3[0], n3[1], n3[2], 0,     0,     0    ],
        [ 0,     0,     0,     0,     0,     0,     0,     0,     0,     n1[0], n1[1], n1[2]],
        [ 0,     0,     0,     0,     0,     0,     0,     0,     0,     n2[0], n2[1], n2[2]],
        [ 0,     0,     0,     0,     0,     0,     0,     0,     0,     n3[0], n3[1], n3[2]]
    ])
    
    Ke = G.T*Kle*G
    fe = G.T*fle
    
    if eq == None:
        return Ke
    else:
        return Ke,fe

def beam3s(ex,ey,ez,eo,ep,ed,eq=None,n=None):
    """
    Calculate the variation of the section forces and displacements
    along a three-dimensional beam element.
    
    Parameters:
     
        ex = [x1 x2]                element node coordinates
        ey = [y1 y2]
        ez = [z1 z2]

        eo = [xz yz zz]             orientation of local z axis
        
        ep = [E G A Iy Iz Kv]       element properties
                                      E: Young's modulus
                                      G: Shear modulus
                                      A: Cross section area
                                     Iy: Moment of inertia, local y-axis
                                     Iz: Moment of inertia, local z-axis
                                     Kv: Saint-Venant's torsion constant

        ed                          the element displacement vector from the
                                    global coordinate system
    
        eq = [qx qy qz qw]          the disibuted axial, transversal and
                                    torsional loads

        n                           the number of point in which displacements
                                    and section forces are to be computed

    Returns:

        es = [[N1,Vy1,Vz1,T1,My1,Mz1],  section forces in n points along
              [N2,Vy2,Vz2,T2,My2,Mz2],  the local x-axis
              [..,...,...,..,...,...],
              [Nn,Vyn,Vzn,Tn,Myn,Mzn]]

        edi = [[u1,v1,w1,fi1],          displacements in n points along
               [u2,v2,w2,fi2],          the local x-axis
               [..,..,..,...],
               [un,vn,wn,fin]]

        eci = [[x1],                    local x-coordinates of the evaluation
               [x2],                    points
               [..],
               [xn]]

    """
    b = mat([
        [ex[1]-ex[0]],
        [ey[1]-ey[0]],
        [ez[1]-ez[0]]
    ])
    L = asscalar(sqrt(b.T*b))
    n1 = asarray(b.T/L).reshape(3,)
    
    eo = asmatrix(eo)
    lc = asscalar(sqrt(eo*eo.T))
    n3 = asarray(eo/lc).reshape(3,)

    EA = ep[0]*ep[2]
    EIy = ep[0]*ep[3]
    EIz = ep[0]*ep[4]
    GKv = ep[1]*ep[5]

    qx = 0.
    qy = 0.
    qz = 0.
    qw = 0.
    if eq != None:
        qx,qy,qz,qw = eq

    ne = 2
    if n != None:
        ne = n
        
    n2 = array([0.,0.,0.])
    n2[0] = n3[1]*n1[2]-n3[2]*n1[1]
    n2[1] = -n1[2]*n3[0]+n1[0]*n3[2]
    n2[2] = n3[0]*n1[1]-n1[0]*n3[1]

    G = mat([
        [ n1[0], n1[1], n1[2], 0,     0,     0,     0,     0,     0,     0,     0,     0    ],
        [ n2[0], n2[1], n2[2], 0,     0,     0,     0,     0,     0,     0,     0,     0    ],
        [ n3[0], n3[1], n3[2], 0,     0,     0,     0,     0,     0,     0,     0,     0    ],
        [ 0,     0,     0,     n1[0], n1[1], n1[2], 0,     0,     0,     0,     0,     0    ],
        [ 0,     0,     0,     n2[0], n2[1], n2[2], 0,     0,     0,     0,     0,     0    ],
        [ 0,     0,     0,     n3[0], n3[1], n3[2], 0,     0,     0,     0,     0,     0    ],
        [ 0,     0,     0,     0,     0,     0,     n1[0], n1[1], n1[2], 0,     0,     0    ],
        [ 0,     0,     0,     0,     0,     0,     n2[0], n2[1], n2[2], 0,     0,     0    ],
        [ 0,     0,     0,     0,     0,     0,     n3[0], n3[1], n3[2], 0,     0,     0    ],
        [ 0,     0,     0,     0,     0,     0,     0,     0,     0,     n1[0], n1[1], n1[2]],
        [ 0,     0,     0,     0,     0,     0,     0,     0,     0,     n2[0], n2[1], n2[2]],
        [ 0,     0,     0,     0,     0,     0,     0,     0,     0,     n3[0], n3[1], n3[2]]
    ])

    u = G*asmatrix(ed).T-array([    # u is the local element displacement
        [ 0               ],        # vector minus the particular solution
        [ 0               ],        # to the beam's diff.eq:s
        [ 0               ],
        [ 0               ],
        [ 0               ],
        [ 0               ],
        [-qx*L**2/(2*EA)  ],
        [ qy*L**4/(24*EIz)],
        [ qz*L**4/(24*EIy)],
        [-qw*L**2/(2*GKv) ],
        [-qz*L**3/(6*EIy) ],
        [ qy*L**3/(6*EIz) ]
    ])

    C = mat([
        [ 0, 1, 0,      0,    0, 0, 0,      0,    0, 0, 0, 0],
        [ 0, 0, 0,      0,    0, 1, 0,      0,    0, 0, 0, 0],
        [ 0, 0, 0,      0,    0, 0, 0,      0,    0, 1, 0, 0],
        [ 0, 0, 0,      0,    0, 0, 0,      0,    0, 0, 0, 1],
        [ 0, 0, 0,      0,    0, 0, 0,      0,   -1, 0, 0, 0],
        [ 0, 0, 0,      0,    1, 0, 0,      0,    0, 0, 0, 0],
        [ L, 1, 0,      0,    0, 0, 0,      0,    0, 0, 0, 0,],
        [ 0, 0, L**3,   L**2, L, 1, 0,      0,    0, 0, 0, 0],
        [ 0, 0, 0,      0,    0, 0, L**3,   L**2, L, 1, 0, 0],
        [ 0, 0, 0,      0,    0, 0, 0,      0,    0, 0, L, 1],
        [ 0, 0, 0,      0,    0, 0,-3*L**2,-2*L, -1, 0, 0, 0],
        [ 0, 0, 3*L**2, 2*L,  1, 0, 0,      0,    0, 0, 0, 0],
    ])

    m = linalg.inv(C)*u
    eci = zeros((ne,1))
    es = zeros((ne,6))
    edi = zeros((ne,4))
    for i in arange(ne):
        x = i*L/(ne-1)
        eci[i,0] = x
        es[i,:] = (mat([
            [ EA, 0, 0,       0,     0, 0, 0,       0,     0, 0, 0,   0],
            [ 0,  0,-6*EIz,   0,     0, 0, 0,       0,     0, 0, 0,   0],
            [ 0,  0, 0,       0,     0, 0,-6*EIy,   0,     0, 0, 0,   0],
            [ 0,  0, 0,       0,     0, 0, 0,       0,     0, 0, GKv, 0],
            [ 0,  0, 0,       0,     0, 0,-6*EIy*x,-2*EIy, 0, 0, 0,   0],
            [ 0,  0, 6*EIz*x, 2*EIz, 0, 0, 0,       0,     0, 0, 0,   0]
        ])*m+array([-qx*x,-qy*x,-qz*x,-qw*x,-qz*x**2/2,qy*x**2/2]).reshape(6,1)).T

        edi[i,:] = (mat([
            [ x, 1, 0,    0,    0, 0, 0,    0,    0, 0, 0, 0],
            [ 0, 0, x**3, x**2, x, 1, 0,    0,    0, 0, 0, 0],
            [ 0, 0, 0,    0,    0, 0, x**3, x**2, x, 1, 0, 0],
            [ 0, 0, 0,    0,    0, 0, 0,    0,    0, 0, x, 1]
        ])*m+array([-qx*x**2/(2*EA),qy*x**4/(24*EIz),qz*x**4/(24*EIy),-qw*x**2/(2*GKv)]).reshape(4,1)).T
    
    if n == None:
        return es
    else:
        return es,edi,eci
    
def flw2te(ex,ey,ep,D,eq=None):
    """
    Compute element stiffness (conductivity) matrix for a triangular field element.
    
    Parameters:
    
        ex = [x1 x2 x3]
        ey = [y1 y2 y3]     element coordinates
    
        ep = [t]            element thickness    

        D = [kxx kxy;
             kyx kyy]       constitutive matrix
    
             eq             heat supply per unit volume
             
    Returns:
    
        Ke                  element 'stiffness' matrix (3 x 3)

        fe                  element load vector (3 x 1)
    
    """
    t=ep[0];
    if eq==None:
        eq=0.
    
    exm = asmatrix(ex)
    eym = asmatrix(ey)
    C=asmatrix(hstack([ones((3,1)),exm.T,eym.T]))
    B=matrix([
        [0.,1.,0.],
        [0.,0.,1.]
    ])*C.I
    A=0.5*linalg.det(C)
  
    Ke=B.T*D*B*t*A
    fe=matrix([[1.,1.,1.]]).T*eq*A*t/3
       
    if eq==0.:
        return Ke
    else:
        return Ke, fe
    
def flw2ts(ex,ey,D,ed):
    """
    Compute flows or corresponding quantities in the triangular field element.
    
    Parameters:
    
        ex = [x1 x2 x3]
        ey = [y1 y2 y3]         element coordinates
                                 
             D = [kxx kxy
                  kyx kyy]      constitutive matrix
    
             ed =[u1 u2 u3]     u1,u2,u3: nodal values
                  .. .. ..;
                  
    Returns:
    
        es=[ qx qy ] 
             ... ..]                element flows
    
        et=[ gx gy ]
             ... ..]                element gradients
    
    """

    if len(ex.shape)>1:
        qs = zeros([ex.shape[0],2])
        qt = zeros([ex.shape[0],2])
        row = 0
        for exr, eyr, edr in zip(ex, ey, ed):
            exm = asmatrix(exr)
            eym = asmatrix(eyr)
            edm = asmatrix(edr)
            C=asmatrix(hstack([ones((3,1)),exm.T,eym.T]))
            B=matrix([
                [0.,1.,0.],
                [0.,0.,1.]
            ])*C.I

            qs[row,:]=(-D*B*edm.T).T
            qt[row,:]=(B*edm.T).T
            row += 1

        return qs, qt
    else:
        exm = asmatrix(ex)
        eym = asmatrix(ey)
        edm = asmatrix(ed)
        C=asmatrix(hstack([ones((3,1)),exm.T,eym.T]))
        B=matrix([
            [0.,1.,0.],
            [0.,0.,1.]
        ])*C.I

        qs=-D*B*edm.T
        qt=B*edm.T
    
        return qs.T, qt.T

def flw2qe(ex,ey,ep,D,eq=None):
    """
    Compute element stiffness (conductivity) matrix for a triangular field element.
    
    Parameters:
    
        ex = [x1, x2, x3, x4]
        ey = [y1, y2, y3, y4]   element coordinates
    
        ep = [t]                element thickness    

        D = [[kxx, kxy],
             [kyx, kyy]]        constitutive matrix
    
        eq                      heat supply per unit volume
             
    Returns:
    
        Ke                      element 'stiffness' matrix (4 x 4)

        fe                      element load vector (4 x 1)
    
    """
    xc = sum(ex)/4.
    yc = sum(ey)/4.

    K = zeros((5,5))
    f = zeros((5,1))
    
    if eq == None:
        k1 = flw2te([ex[0],ex[1],xc],[ey[0],ey[1],yc],ep,D)
        K = assem(array([1,2,5]),K,k1)
        k1 = flw2te([ex[1],ex[2],xc],[ey[1],ey[2],yc],ep,D)
        K = assem(array([2,3,5]),K,k1)
        k1 = flw2te([ex[2],ex[3],xc],[ey[2],ey[3],yc],ep,D)
        K = assem(array([3,4,5]),K,k1)
        k1 = flw2te([ex[3],ex[0],xc],[ey[3],ey[0],yc],ep,D)
        K = assem(array([4,1,5]),K,k1)
    else:
        k1,f1 = flw2te([ex[0],ex[1],xc],[ey[0],ey[1],yc],ep,D,eq)    
        K,f = assem(array([1,2,5]),K,k1,f,f1)
        k1,f1 = flw2te([ex[1],ex[2],xc],[ey[1],ey[2],yc],ep,D,eq)
        K,f = assem(array([2,3,5]),K,k1,f,f1)
        k1,f1 = flw2te([ex[2],ex[3],xc],[ey[2],ey[3],yc],ep,D,eq)
        K,f = assem(array([3,4,5]),K,k1,f,f1)
        k1,f1 = flw2te([ex[3],ex[0],xc],[ey[3],ey[0],yc],ep,D,eq)
        K,f = assem(array([4,1,5]),K,k1,f,f1)
    Ke1,fe1 = statcon(K,f,array([5]));

    Ke = Ke1
    fe = fe1
    
    if eq == None:
        return Ke
    else:
        return Ke,fe

def flw2qs(ex,ey,ep,D,ed,eq=None):
    """
    Compute flows or corresponding quantities in the
    quadrilateral field element.
    
    Parameters:
    
        ex = [x1, x2, x3, x4]
        ey = [y1, y2, y3, y4]      element coordinates
    
        ep = [t]                   element thickness    

        D = [[kxx, kxy],
             [kyx, kyy]]           constitutive matrix

        ed = [[u1, u2, u3, u4],
              [.., .., .., ..]]    u1,u2,u3,u4: nodal values
    
        eq                         heat supply per unit volume
             
    Returns:
    
        es = [[qx, qy],
              [.., ..]]            element flows

        et = [[gx, gy],
              [.., ..]]            element gradients
    
    """
    K = zeros((5,5))
    f = zeros((5,1))
    
    xm = sum(ex)/4
    ym = sum(ey)/4
    
    if eq == None:
        q = 0
    else:
        q = eq
    
    En = array([
        [1,2,5],
        [2,3,5],
        [3,4,5],
        [4,1,5]
    ])
    ex1 = array([ex[0],ex[1],xm])
    ey1 = array([ey[0],ey[1],ym])
    ex2 = array([ex[1],ex[2],xm])
    ey2 = array([ey[1],ey[2],ym])
    ex3 = array([ex[2],ex[3],xm])
    ey3 = array([ey[2],ey[3],ym])
    ex4 = array([ex[3],ex[0],xm])
    ey4 = array([ey[3],ey[0],ym])
    
    if eq == None:
        k1 = flw2te(ex1,ey1,ep,D)
        K = assem(En[0],K,k1)
        k1 = flw2te(ex2,ey2,ep,D)
        K = assem(En[1],K,k1)
        k1 = flw2te(ex3,ey3,ep,D)
        K = assem(En[2],K,k1)
        k1 = flw2te(ex4,ey4,ep,D)
        K = assem(En[3],K,k1)
    else:
        k1,f1 = flw2te(ex1,ey1,ep,D,q)
        K,f = assem(En[0],K,k1,f,f1)
        k1,f1 = flw2te(ex2,ey2,ep,D,q)
        K,f = assem(En[1],K,k1,f,f1)
        k1,f1 = flw2te(ex3,ey3,ep,D,q)
        K,f = assem(En[2],K,k1,f,f1)
        k1,f1 = flw2te(ex4,ey4,ep,D,q)
        K,f = assem(En[3],K,k1,f,f1)
    
    if rank(ed)==1:
        ed = array([ed])
    ni,nj = shape(ed)

    a = zeros((5,ni))
    for i in range(ni):
        a[ix_(range(5),[i])],r = asarray(solveq(K,f,arange(1,5),ed[i]))

    s1,t1 = flw2ts(ex1,ey1,D,a[ix_(En[0,:]-1,arange(ni))].T)
    s2,t2 = flw2ts(ex2,ey2,D,a[ix_(En[1,:]-1,arange(ni))].T)
    s3,t3 = flw2ts(ex3,ey3,D,a[ix_(En[2,:]-1,arange(ni))].T)
    s4,t4 = flw2ts(ex4,ey4,D,a[ix_(En[3,:]-1,arange(ni))].T)
    
    es = (s1+s2+s3+s4)/4.
    et = (t1+t2+t3+t4)/4.
    
    return es,et

def flw2i4e(ex,ey,ep,D,eq=None):
    """
    Compute element stiffness (conductivity)
    matrix for 4 node isoparametric field element

    Parameters:
        
        ex = [x1 x2 x3 x4]  element coordinates
        ey = [y1 y2 y3 y4]

        ep = [t ir]         thickness and integration rule

        D  = [[kxx kxy],
              [kyx kyy]]    constitutive matrix

        eq                  heat supply per unit volume

    Returns:
        Ke                  element 'stiffness' matrix (4 x 4)
        fe                  element load vector (4 x 1)

    """
    t = ep[0]
    ir = ep[1]
    ngp = ir*ir
    
    if eq == None:
        q = 0
    else:
        q = eq

    if ir == 1:
        g1 = 0.0
        w1 = 2.0
        gp = mat([g1,g1])
        w = mat([w1,w1])
    elif ir == 2:
        g1 = 0.577350269189626
        w1 = 1
        gp = mat([
            [-g1,-g1],
            [ g1,-g1],
            [-g1, g1],
            [ g1, g1]
        ])
        w = mat([
            [ w1, w1],
            [ w1, w1],
            [ w1, w1],
            [ w1, w1]
        ])
    elif ir == 3:
        g1 = 0.774596669241483
        g2 = 0.
        w1 = 0.555555555555555
        w2 = 0.888888888888888
        gp = mat([
            [-g1,-g1],
            [-g2,-g1],
            [ g1,-g1],
            [-g1, g2],
            [ g2, g2],
            [ g1, g2],
            [-g1, g1],
            [ g2, g1],
            [ g1, g1]
        ])
        w = mat([
            [ w1, w1],
            [ w2, w1],
            [ w1, w1],
            [ w1, w2],
            [ w2, w2],
            [ w1, w2],
            [ w1, w1],
            [ w2, w1],
            [ w1, w1]
        ])
    else:
        print("Used number of integration points not implemented")
    wp = multiply(w[:,0],w[:,1])
    
    xsi = gp[:,0]
    eta = gp[:,1]
    r2 = ngp*2

    N = multiply((1-xsi),(1-eta))/4.
    N = append(N,multiply((1+xsi),(1-eta))/4.,axis=1)
    N = append(N,multiply((1+xsi),(1+eta))/4.,axis=1)
    N = append(N,multiply((1-xsi),(1+eta))/4.,axis=1)
    
    dNr = mat(zeros((r2,4)))
    dNr[0:r2:2,0] = -(1-eta)/4.
    dNr[0:r2:2,1] = (1-eta)/4.
    dNr[0:r2:2,2] = (1+eta)/4.
    dNr[0:r2:2,3] = -(1+eta)/4.
    dNr[1:r2+1:2,0] = -(1-xsi)/4.
    dNr[1:r2+1:2,1] = -(1+xsi)/4.
    dNr[1:r2+1:2,2] = (1+xsi)/4.
    dNr[1:r2+1:2,3] = (1-xsi)/4.

    Ke1 = mat(zeros((4,4)))
    fe1 = mat(zeros((4,1)))
    JT = dNr*mat([ex,ey]).T

    for i in range(ngp):
        indx = array([2*(i+1)-1,2*(i+1)])
        detJ = linalg.det(JT[indx-1,:])
        if detJ < 10*finfo(float).eps:
            print("Jacobideterminanten lika med noll!")
        JTinv = linalg.inv(JT[indx-1,:])
        B = JTinv*dNr[indx-1,:]
        Ke1 = Ke1+B.T*D*B*detJ*asscalar(wp[i])
        fe1 = fe1+N[i,:].T*detJ*wp[i]

    if eq == None:
        return Ke1*t
    else:
        return Ke1*t,fe1*t*eq

def flw2i4s(ex,ey,ep,D,ed):
    """
    Compute flows or corresponding quantities in the
    4 node isoparametric element.
    
    Parameters:
        
        ex = [x1 x2 x3 x4]         element coordinates
        ey = [y1 y2 y3 y4]

        ep = [t ir]                thickness and integration rule

        D  = [[kxx kxy],
              [kyx kyy]]           constitutive matrix

        ed = [u1, u2, u3, u4]      u1,u2,u3,u4: nodal values

    Returns:
        es = [[qx, qy],
              [.., ..]]             element flows

        et = [[qx, qy],
              [... ..]]             element gradients

        eci=[[ix1, iy1],            Gauss point location vector
             [...  ...],            nint: number of integration points
             [ix(nint), iy(nint)]

    """
    t = ep[0]
    ir = ep[1]
    ngp = ir*ir

    if ir == 1:
        g1 = 0.0
        w1 = 2.0
        gp = mat([g1,g1])
        w = mat([w1,w1])
    elif ir == 2:
        g1 = 0.577350269189626
        w1 = 1
        gp = mat([
            [-g1,-g1],
            [ g1,-g1],
            [-g1, g1],
            [ g1, g1]
        ])
        w = mat([
            [ w1, w1],
            [ w1, w1],
            [ w1, w1],
            [ w1, w1]
        ])
    elif ir == 3:
        g1 = 0.774596669241483
        g2 = 0.
        w1 = 0.555555555555555
        w2 = 0.888888888888888
        gp = mat([
            [-g1,-g1],
            [-g2,-g1],
            [ g1,-g1],
            [-g1, g2],
            [ g2, g2],
            [ g1, g2],
            [-g1, g1],
            [ g2, g1],
            [ g1, g1]
        ])
        w = mat([
            [ w1, w1],
            [ w2, w1],
            [ w1, w1],
            [ w1, w2],
            [ w2, w2],
            [ w1, w2],
            [ w1, w1],
            [ w2, w1],
            [ w1, w1]
        ])
    else:
        print("Used number of integration points not implemented")
    wp = multiply(w[:,0],w[:,1])

    xsi = gp[:,0]
    eta = gp[:,1]
    r2 = ngp*2

    N = multiply((1-xsi),(1-eta))/4.
    N = append(N,multiply((1+xsi),(1-eta))/4.,axis=1)
    N = append(N,multiply((1+xsi),(1+eta))/4.,axis=1)
    N = append(N,multiply((1-xsi),(1+eta))/4.,axis=1)
    
    dNr = mat(zeros((r2,4)))
    dNr[0:r2:2,0] = -(1-eta)/4.
    dNr[0:r2:2,1] = (1-eta)/4.
    dNr[0:r2:2,2] = (1+eta)/4.
    dNr[0:r2:2,3] = -(1+eta)/4.
    dNr[1:r2+1:2,0] = -(1-xsi)/4.
    dNr[1:r2+1:2,1] = -(1+xsi)/4.
    dNr[1:r2+1:2,2] = (1+xsi)/4.
    dNr[1:r2+1:2,3] = (1-xsi)/4.

    eci = N*mat([ex,ey]).T
    if rank(ed) == 1:
        ed = array([ed])
    red,ced = shape(ed)
    JT = dNr*mat([ex,ey]).T
    
    es = mat(zeros((ngp*red,2)))
    et = mat(zeros((ngp*red,2)))
    for i in range(ngp):
        indx = array([2*(i+1)-1,2*(i+1)])
        detJ = linalg.det(JT[indx-1,:])
        if detJ < 10*finfo(float).eps:
            print("Jacobideterminanten lika med noll!")
        JTinv = linalg.inv(JT[indx-1,:])
        B = JTinv*dNr[indx-1,:]
        p1 = -D*B*ed.T
        p2 = B*ed.T
        es[i:ngp*red:ngp,:] = p1.T
        et[i:ngp*red:ngp,:] = p2.T
    
    return es,et,eci

def flw2i8e(ex,ey,ep,D,eq=None):
    """
    Compute element stiffness (conductivity)
    matrix for 8 node isoparametric field element.
    
    Parameters:
    
        ex = [x1, ..., x8]      element coordinates
        ey = [y1, ..., y8]
        
        ep = [t, ir]            thickness and integration rule

        D = [[kxx, kxy],
             [kyx, kyy]]        constitutive matrix

        eq                      heat supply per unit volume

    Returns:
    
        Ke                      element 'stiffness' matrix (8 x 8)
        fe                      element load vector (8 x 1)

    """
    t = ep[0]
    ir = ep[1]
    ngp = ir*ir

    if eq == None:
        q = 0
    else:
        q = eq

    if ir == 1:
        g1 = 0.0
        w1 = 2.0
        gp = mat([g1,g1])
        w = mat([w1,w1])
    elif ir == 2:
        g1 = 0.577350269189626
        w1 = 1
        gp = mat([
            [-g1,-g1],
            [ g1,-g1],
            [-g1, g1],
            [ g1, g1]
        ])
        w = mat([
            [ w1, w1],
            [ w1, w1],
            [ w1, w1],
            [ w1, w1]
        ])
    elif ir == 3:
        g1 = 0.774596669241483
        g2 = 0.
        w1 = 0.555555555555555
        w2 = 0.888888888888888
        gp = mat([
            [-g1,-g1],
            [-g2,-g1],
            [ g1,-g1],
            [-g1, g2],
            [ g2, g2],
            [ g1, g2],
            [-g1, g1],
            [ g2, g1],
            [ g1, g1]
        ])
        w = mat([
            [ w1, w1],
            [ w2, w1],
            [ w1, w1],
            [ w1, w2],
            [ w2, w2],
            [ w1, w2],
            [ w1, w1],
            [ w2, w1],
            [ w1, w1]
        ])
    else:
        print("Used number of integration points not implemented")
    wp = multiply(w[:,0],w[:,1])

    xsi = gp[:,0]
    eta = gp[:,1]
    r2 = ngp*2

    N = multiply(multiply(-(1-xsi),(1-eta)),(1+xsi+eta))/4.
    N = append(N,multiply(multiply(-(1+xsi),(1-eta)),(1-xsi+eta))/4.,axis=1)
    N = append(N,multiply(multiply(-(1+xsi),(1+eta)),(1-xsi-eta))/4.,axis=1)
    N = append(N,multiply(multiply(-(1-xsi),(1+eta)),(1+xsi-eta))/4.,axis=1)
    N = append(N,multiply((1-multiply(xsi,xsi)),(1-eta))/2.,axis=1)
    N = append(N,multiply((1+xsi),(1-multiply(eta,eta)))/2.,axis=1)
    N = append(N,multiply((1-multiply(xsi,xsi)),(1+eta))/2.,axis=1)
    N = append(N,multiply((1-xsi),(1-multiply(eta,eta)))/2.,axis=1)

    dNr = mat(zeros((r2,8)))
    dNr[0:r2:2,0] = -(-multiply((1-eta),(1+xsi+eta))+multiply((1-xsi),(1-eta)))/4.
    dNr[0:r2:2,1] = -(multiply((1-eta),(1-xsi+eta))-multiply((1+xsi),(1-eta)))/4.
    dNr[0:r2:2,2] = -(multiply((1+eta),(1-xsi-eta))-multiply((1+xsi),(1+eta)))/4.
    dNr[0:r2:2,3] = -(-multiply((1+eta),(1+xsi-eta))+multiply((1-xsi),(1+eta)))/4.
    dNr[0:r2:2,4] = -multiply(xsi,(1-eta))
    dNr[0:r2:2,5] = (1-multiply(eta,eta))/2.
    dNr[0:r2:2,6] = -multiply(xsi,(1+eta))
    dNr[0:r2:2,7] = -(1-multiply(eta,eta))/2.
    dNr[1:r2+1:2,0] = -(-multiply((1-xsi),(1+xsi+eta))+multiply((1-xsi),(1-eta)))/4.
    dNr[1:r2+1:2,1] = -(-multiply((1+xsi),(1-xsi+eta))+multiply((1+xsi),(1-eta)))/4.
    dNr[1:r2+1:2,2] = -(multiply((1+xsi),(1-xsi-eta))-multiply((1+xsi),(1+eta)))/4.
    dNr[1:r2+1:2,3] = -(multiply((1-xsi),(1+xsi-eta))-multiply((1-xsi),(1+eta)))/4.
    dNr[1:r2+1:2,4] = -(1-multiply(xsi,xsi))/2.
    dNr[1:r2+1:2,5] = -multiply(eta,(1+xsi))
    dNr[1:r2+1:2,6] = (1-multiply(xsi,xsi))/2.
    dNr[1:r2+1:2,7] = -multiply(eta,(1-xsi))

    Ke1 = mat(zeros((8,8)))
    fe1 = mat(zeros((8,1)))
    JT = dNr*mat([ex,ey]).T

    for i in range(ngp):
        indx = array([2*(i+1)-1,2*(i+1)])
        detJ = linalg.det(JT[indx-1,:])
        if detJ < 10*finfo(float).eps:
            print("Jacobideterminanten lika med noll!")
        JTinv = linalg.inv(JT[indx-1,:])
        B = JTinv*dNr[indx-1,:]
        Ke1 = Ke1+B.T*D*B*detJ*asscalar(wp[i])
        fe1 = fe1+N[i,:].T*detJ*wp[i]
    
    if eq != None:
        return Ke1*t,fe1*t*q
    else:
        return Ke1*t

def flw2i8s(ex,ey,ep,D,ed):
    """
    Compute flows or corresponding quantities in the
    8 node isoparametric element.
    
    Parameters:
        
        ex = [x1,x2,x3....,x8]     element coordinates
        ey = [y1,y2,y3....,y8]

        ep = [t,ir]                thickness and integration rule

        D  = [[kxx,kxy],
              [kyx,kyy]]           constitutive matrix

        ed = [u1,....,u8]          u1,....,u8: nodal values

    Returns:
        es = [[qx,qy],
              [..,..]]             element flows

        et = [[qx,qy],
              [..,..]]             element gradients

        eci=[[ix1,iy1],            Gauss point location vector
             [...,...],            nint: number of integration points
             [ix(nint),iy(nint)]]

    """
    t = ep[0]
    ir = ep[1]
    ngp = ir*ir

    if ir == 1:
        g1 = 0.0
        w1 = 2.0
        gp = mat([g1,g1])
        w = mat([w1,w1])
    elif ir == 2:
        g1 = 0.577350269189626
        w1 = 1
        gp = mat([
            [-g1,-g1],
            [ g1,-g1],
            [-g1, g1],
            [ g1, g1]
        ])
        w = mat([
            [ w1, w1],
            [ w1, w1],
            [ w1, w1],
            [ w1, w1]
        ])
    elif ir == 3:
        g1 = 0.774596669241483
        g2 = 0.
        w1 = 0.555555555555555
        w2 = 0.888888888888888
        gp = mat([
            [-g1,-g1],
            [-g2,-g1],
            [ g1,-g1],
            [-g1, g2],
            [ g2, g2],
            [ g1, g2],
            [-g1, g1],
            [ g2, g1],
            [ g1, g1]
        ])
        w = mat([
            [ w1, w1],
            [ w2, w1],
            [ w1, w1],
            [ w1, w2],
            [ w2, w2],
            [ w1, w2],
            [ w1, w1],
            [ w2, w1],
            [ w1, w1]
        ])
    else:
        print("Used number of integration points not implemented")
    wp = multiply(w[:,0],w[:,1])

    xsi = gp[:,0]
    eta = gp[:,1]
    r2 = ngp*2

    N = multiply(multiply(-(1-xsi),(1-eta)),(1+xsi+eta))/4.
    N = append(N,multiply(multiply(-(1+xsi),(1-eta)),(1-xsi+eta))/4.,axis=1)
    N = append(N,multiply(multiply(-(1+xsi),(1+eta)),(1-xsi-eta))/4.,axis=1)
    N = append(N,multiply(multiply(-(1-xsi),(1+eta)),(1+xsi-eta))/4.,axis=1)
    N = append(N,multiply((1-multiply(xsi,xsi)),(1-eta))/2.,axis=1)
    N = append(N,multiply((1+xsi),(1-multiply(eta,eta)))/2.,axis=1)
    N = append(N,multiply((1-multiply(xsi,xsi)),(1+eta))/2.,axis=1)
    N = append(N,multiply((1-xsi),(1-multiply(eta,eta)))/2.,axis=1)

    dNr = mat(zeros((r2,8)))
    dNr[0:r2:2,0] = -(-multiply((1-eta),(1+xsi+eta))+multiply((1-xsi),(1-eta)))/4.
    dNr[0:r2:2,1] = -(multiply((1-eta),(1-xsi+eta))-multiply((1+xsi),(1-eta)))/4.
    dNr[0:r2:2,2] = -(multiply((1+eta),(1-xsi-eta))-multiply((1+xsi),(1+eta)))/4.
    dNr[0:r2:2,3] = -(-multiply((1+eta),(1+xsi-eta))+multiply((1-xsi),(1+eta)))/4.
    dNr[0:r2:2,4] = -multiply(xsi,(1-eta))
    dNr[0:r2:2,5] = (1-multiply(eta,eta))/2.
    dNr[0:r2:2,6] = -multiply(xsi,(1+eta))
    dNr[0:r2:2,7] = -(1-multiply(eta,eta))/2.
    dNr[1:r2+1:2,0] = -(-multiply((1-xsi),(1+xsi+eta))+multiply((1-xsi),(1-eta)))/4.
    dNr[1:r2+1:2,1] = -(-multiply((1+xsi),(1-xsi+eta))+multiply((1+xsi),(1-eta)))/4.
    dNr[1:r2+1:2,2] = -(multiply((1+xsi),(1-xsi-eta))-multiply((1+xsi),(1+eta)))/4.
    dNr[1:r2+1:2,3] = -(multiply((1-xsi),(1+xsi-eta))-multiply((1-xsi),(1+eta)))/4.
    dNr[1:r2+1:2,4] = -(1-multiply(xsi,xsi))/2.
    dNr[1:r2+1:2,5] = -multiply(eta,(1+xsi))
    dNr[1:r2+1:2,6] = (1-multiply(xsi,xsi))/2.
    dNr[1:r2+1:2,7] = -multiply(eta,(1-xsi))

    eci = N*mat([ex,ey]).T
    if rank(ed) == 1:
        ed = array([ed])
    red,ced = shape(ed)
    JT = dNr*mat([ex,ey]).T
    
    es = mat(zeros((ngp*red,2)))
    et = mat(zeros((ngp*red,2)))
    for i in range(ngp):
        indx = array([2*(i+1)-1,2*(i+1)])
        detJ = linalg.det(JT[indx-1,:])
        if detJ < 10*finfo(float).eps:
            print("Jacobideterminanten lika med noll!")
        JTinv = linalg.inv(JT[indx-1,:])
        B = JTinv*dNr[indx-1,:]
        p1 = -D*B*ed.T
        p2 = B*ed.T
        es[i:ngp*red:ngp,:] = p1.T
        et[i:ngp*red:ngp,:] = p2.T

    return es,et,eci

def flw3i8e(ex,ey,ez,ep,D,eq=None):
    """
    Compute element stiffness (conductivity)
    matrix for 8 node isoparametric field element.
    
    Parameters:
    
        ex = [x1,x2,x3,...,x8]
        ey = [y1,y2,y3,...,y8]      element coordinates
        ez = [z1,z2,z3,...,z8]

        ep = [ir]                   Ir: Integration rule

        D = [[kxx,kxy,kxz],
             [kyx,kyy,kyz],
             [kzx,kzy,kzz]]         constitutive matrix

        eq                          heat supply per unit volume

    Output:

        Ke                          element 'stiffness' matrix (8 x 8)
        fe                          element load vector (8 x 1)

    """
    ir = ep[0]
    ngp = ir*ir*ir
    
    if eq == None:
        q = 0
    else:
        q = eq
    
    if ir == 2:
        g1 = 0.577350269189626
        w1 = 1
        gp = mat([
            [-1,-1,-1],
            [ 1,-1,-1],
            [ 1, 1,-1],
            [-1, 1,-1],
            [-1,-1, 1],
            [ 1,-1, 1],
            [ 1, 1, 1],
            [-1, 1, 1]
        ])*g1
        w = mat(ones((8,3)))*w1
    elif ir == 3:
        g1 = 0.774596669241483
        g2 = 0.
        w1 = 0.555555555555555
        w2 = 0.888888888888888
        gp = mat(zeros((27,3)))
        w = mat(zeros((27,3)))
        I1 = array([-1,0,1,-1,0,1,-1,0,1])
        I2 = array([0,-1,0,0,1,0,0,1,0])
        gp[:,0] = mat([I1,I1,I1]).reshape(27,1)*g1
        gp[:,0] = mat([I2,I2,I2]).reshape(27,1)*g2+gp[:,0]
        I1 = abs(I1)
        I2 = abs(I2)
        w[:,0] = mat([I1,I1,I1]).reshape(27,1)*w1
        w[:,0] = mat([I2,I2,I2]).reshape(27,1)*w2+w[:,0]
        I1 = array([-1,-1,-1,0,0,0,1,1,1])
        I2 = array([0,0,0,1,1,1,0,0,0])
        gp[:,1] = mat([I1,I1,I1]).reshape(27,1)*g1
        gp[:,1] = mat([I2,I2,I2]).reshape(27,1)*g2+gp[:,1]
        I1 = abs(I1)
        I2 = abs(I2)
        w[:,1] = mat([I1,I1,I1]).reshape(27,1)*w1
        w[:,1] = mat([I2,I2,I2]).reshape(27,1)*w2+w[:,1]
        I1 = array([-1,-1,-1,-1,-1,-1,-1,-1,-1])
        I2 = array([0,0,0,0,0,0,0,0,0])
        I3 = abs(I1)
        gp[:,2] = mat([I1,I2,I3]).reshape(27,1)*g1
        gp[:,2] = mat([I2,I3,I2]).reshape(27,1)*g2+gp[:,2]
        w[:,2] = mat([I3,I2,I3]).reshape(27,1)*w1
        w[:,2] = mat([I2,I3,I2]).reshape(27,1)*w2+w[:,2]
    else:
        print("Used number of integration points not implemented")
        return

    wp = multiply(multiply(w[:,0],w[:,1]),w[:,2])

    xsi = gp[:,0]
    eta = gp[:,1]
    zet = gp[:,2]
    r2 = ngp*3
    
    N = multiply(multiply((1-xsi),(1-eta)),(1-zet))/8.
    N = append(N,multiply(multiply((1+xsi),(1-eta)),(1-zet))/8.,axis=1)
    N = append(N,multiply(multiply((1+xsi),(1+eta)),(1-zet))/8.,axis=1)
    N = append(N,multiply(multiply((1-xsi),(1+eta)),(1-zet))/8.,axis=1)
    N = append(N,multiply(multiply((1-xsi),(1-eta)),(1+zet))/8.,axis=1)
    N = append(N,multiply(multiply((1+xsi),(1-eta)),(1+zet))/8.,axis=1)
    N = append(N,multiply(multiply((1+xsi),(1+eta)),(1+zet))/8.,axis=1)
    N = append(N,multiply(multiply((1-xsi),(1+eta)),(1+zet))/8.,axis=1)
    
    dNr = mat(zeros((r2,8)))
    dNr[0:r2:3,0]= multiply(-(1-eta),(1-zet))
    dNr[0:r2:3,1]= multiply((1-eta),(1-zet))
    dNr[0:r2:3,2]= multiply((1+eta),(1-zet))
    dNr[0:r2:3,3]= multiply(-(1+eta),(1-zet))
    dNr[0:r2:3,4]= multiply(-(1-eta),(1+zet))
    dNr[0:r2:3,5]= multiply((1-eta),(1+zet))
    dNr[0:r2:3,6]= multiply((1+eta),(1+zet))
    dNr[0:r2:3,7]= multiply(-(1+eta),(1+zet))
    dNr[1:r2+1:3,0] = multiply(-(1-xsi),(1-zet))
    dNr[1:r2+1:3,1] = multiply(-(1+xsi),(1-zet))
    dNr[1:r2+1:3,2] = multiply((1+xsi),(1-zet))
    dNr[1:r2+1:3,3] = multiply((1-xsi),(1-zet))
    dNr[1:r2+1:3,4] = multiply(-(1-xsi),(1+zet))
    dNr[1:r2+1:3,5] = multiply(-(1+xsi),(1+zet))
    dNr[1:r2+1:3,6] = multiply((1+xsi),(1+zet))
    dNr[1:r2+1:3,7] = multiply((1-xsi),(1+zet))
    dNr[2:r2+2:3,0] = multiply(-(1-xsi),(1-eta))
    dNr[2:r2+2:3,1] = multiply(-(1+xsi),(1-eta))
    dNr[2:r2+2:3,2] = multiply(-(1+xsi),(1+eta))
    dNr[2:r2+2:3,3] = multiply(-(1-xsi),(1+eta))
    dNr[2:r2+2:3,4] = multiply((1-xsi),(1-eta))
    dNr[2:r2+2:3,5] = multiply((1+xsi),(1-eta))
    dNr[2:r2+2:3,6] = multiply((1+xsi),(1+eta))
    dNr[2:r2+2:3,7] = multiply((1-xsi),(1+eta))
    dNr = dNr/8.

    Ke1 = mat(zeros((8,8)))
    fe1 = mat(zeros((8,1)))
    JT = dNr*mat([ex,ey,ez]).T
    
    for i in range(ngp):
        indx = array([3*(i+1)-2,3*(i+1)-1,3*(i+1)])
        detJ = linalg.det(JT[indx-1,:])
        if detJ < 10*finfo(float).eps:
            print("Jacobideterminanten lika med noll!")
        JTinv = linalg.inv(JT[indx-1,:])
        B = JTinv*dNr[indx-1,:]
        Ke1 = Ke1+B.T*D*B*detJ*asscalar(wp[i])
        fe1 = fe1+N[i,:].T*detJ*wp[i]

    if eq != None:
        return Ke1,fe1*q
    else:
        return Ke1

def flw3i8s(ex,ey,ez,ep,D,ed):
    """
    Compute flows or corresponding quantities in the
    8 node (3-dim) isoparametric field element.
    
    Parameters:
    
        ex = [x1,x2,x3,...,x8]
        ey = [y1,y2,y3,...,y8]              element coordinates
        ez = [z1,z2,z3,...,z8]

        ep = [ir]                           Ir: Integration rule

        D = [[kxx,kxy,kxz],
             [kyx,kyy,kyz],
             [kzx,kzy,kzz]]                 constitutive matrix

        ed = [[u1,....,u8],                 element nodal values
              [..,....,..]]

    Output:

        es = [[qx,qy,qz],
              [..,..,..]]                   element flows(s)

        et = [[qx,qy,qz],                   element gradients(s)
              [..,..,..]]

        eci = [[ix1,ix1,iz1],               location vector
               [...,...,...],               nint: number of integration points
               [ix(nint),iy(nint),iz(nint)]]

    """
    ir = ep[0]
    ngp = ir*ir*ir
    
    if ir == 2:
        g1 = 0.577350269189626
        w1 = 1
        gp = mat([
            [-1,-1,-1],
            [ 1,-1,-1],
            [ 1, 1,-1],
            [-1, 1,-1],
            [-1,-1, 1],
            [ 1,-1, 1],
            [ 1, 1, 1],
            [-1, 1, 1]
        ])*g1
        w = mat(ones((8,3)))*w1
    elif ir == 3:
        g1 = 0.774596669241483
        g2 = 0.
        w1 = 0.555555555555555
        w2 = 0.888888888888888
        gp = mat(zeros((27,3)))
        w = mat(zeros((27,3)))
        I1 = array([-1,0,1,-1,0,1,-1,0,1])
        I2 = array([0,-1,0,0,1,0,0,1,0])
        gp[:,0] = mat([I1,I1,I1]).reshape(27,1)*g1
        gp[:,0] = mat([I2,I2,I2]).reshape(27,1)*g2+gp[:,0]
        I1 = abs(I1)
        I2 = abs(I2)
        w[:,0] = mat([I1,I1,I1]).reshape(27,1)*w1
        w[:,0] = mat([I2,I2,I2]).reshape(27,1)*w2+w[:,0]
        I1 = array([-1,-1,-1,0,0,0,1,1,1])
        I2 = array([0,0,0,1,1,1,0,0,0])
        gp[:,1] = mat([I1,I1,I1]).reshape(27,1)*g1
        gp[:,1] = mat([I2,I2,I2]).reshape(27,1)*g2+gp[:,1]
        I1 = abs(I1)
        I2 = abs(I2)
        w[:,1] = mat([I1,I1,I1]).reshape(27,1)*w1
        w[:,1] = mat([I2,I2,I2]).reshape(27,1)*w2+w[:,1]
        I1 = array([-1,-1,-1,-1,-1,-1,-1,-1,-1])
        I2 = array([0,0,0,0,0,0,0,0,0])
        I3 = abs(I1)
        gp[:,2] = mat([I1,I2,I3]).reshape(27,1)*g1
        gp[:,2] = mat([I2,I3,I2]).reshape(27,1)*g2+gp[:,2]
        w[:,2] = mat([I3,I2,I3]).reshape(27,1)*w1
        w[:,2] = mat([I2,I3,I2]).reshape(27,1)*w2+w[:,2]
    else:
        print("Used number of integration points not implemented")
        return

    wp = multiply(multiply(w[:,0],w[:,1]),w[:,2])
    
    xsi = gp[:,0]
    eta = gp[:,1]
    zet = gp[:,2]
    r2 = ngp*3
    
    N = multiply(multiply((1-xsi),(1-eta)),(1-zet))/8.
    N = append(N,multiply(multiply((1+xsi),(1-eta)),(1-zet))/8.,axis=1)
    N = append(N,multiply(multiply((1+xsi),(1+eta)),(1-zet))/8.,axis=1)
    N = append(N,multiply(multiply((1-xsi),(1+eta)),(1-zet))/8.,axis=1)
    N = append(N,multiply(multiply((1-xsi),(1-eta)),(1+zet))/8.,axis=1)
    N = append(N,multiply(multiply((1+xsi),(1-eta)),(1+zet))/8.,axis=1)
    N = append(N,multiply(multiply((1+xsi),(1+eta)),(1+zet))/8.,axis=1)
    N = append(N,multiply(multiply((1-xsi),(1+eta)),(1+zet))/8.,axis=1)
    
    dNr = mat(zeros((r2,8)))
    dNr[0:r2:3,0]= multiply(-(1-eta),(1-zet))
    dNr[0:r2:3,1]= multiply((1-eta),(1-zet))
    dNr[0:r2:3,2]= multiply((1+eta),(1-zet))
    dNr[0:r2:3,3]= multiply(-(1+eta),(1-zet))
    dNr[0:r2:3,4]= multiply(-(1-eta),(1+zet))
    dNr[0:r2:3,5]= multiply((1-eta),(1+zet))
    dNr[0:r2:3,6]= multiply((1+eta),(1+zet))
    dNr[0:r2:3,7]= multiply(-(1+eta),(1+zet))
    dNr[1:r2+1:3,0] = multiply(-(1-xsi),(1-zet))
    dNr[1:r2+1:3,1] = multiply(-(1+xsi),(1-zet))
    dNr[1:r2+1:3,2] = multiply((1+xsi),(1-zet))
    dNr[1:r2+1:3,3] = multiply((1-xsi),(1-zet))
    dNr[1:r2+1:3,4] = multiply(-(1-xsi),(1+zet))
    dNr[1:r2+1:3,5] = multiply(-(1+xsi),(1+zet))
    dNr[1:r2+1:3,6] = multiply((1+xsi),(1+zet))
    dNr[1:r2+1:3,7] = multiply((1-xsi),(1+zet))
    dNr[2:r2+2:3,0] = multiply(-(1-xsi),(1-eta))
    dNr[2:r2+2:3,1] = multiply(-(1+xsi),(1-eta))
    dNr[2:r2+2:3,2] = multiply(-(1+xsi),(1+eta))
    dNr[2:r2+2:3,3] = multiply(-(1-xsi),(1+eta))
    dNr[2:r2+2:3,4] = multiply((1-xsi),(1-eta))
    dNr[2:r2+2:3,5] = multiply((1+xsi),(1-eta))
    dNr[2:r2+2:3,6] = multiply((1+xsi),(1+eta))
    dNr[2:r2+2:3,7] = multiply((1-xsi),(1+eta))
    dNr = dNr/8.

    eci = N*mat([ex,ey,ez]).T
    if rank(ed) == 1:
        ed = array([ed])
        red,ced = shape(ed)
    JT = dNr*mat([ex,ey,ez]).T

    es = mat(zeros((ngp*red,3)))
    et = mat(zeros((ngp*red,3)))
    for i in range(ngp):
        indx = array([3*(i+1)-2,3*(i+1)-1,3*(i+1)])
        detJ = linalg.det(JT[indx-1,:])
        if detJ < 10*finfo(float).eps:
            print("Jacobideterminanten lika med noll!")
        JTinv = linalg.inv(JT[indx-1,:])
        B = JTinv*dNr[indx-1,:]
        p1 = -D*B*ed.T
        p2 = B*ed.T
        es[i:ngp*red:ngp,:] = p1.T
        et[i:ngp*red:ngp,:] = p2.T

    return es,et,eci

def plante(ex,ey,ep,D,eq=None):
    """
    Calculate the stiffness matrix for a triangular plane stress or plane strain element.
    
    Parameters:
    
        ex = [x1,x2,x3]         element coordinates
        ey = [y1,y2,y3]
     
        ep = [ptype,t]          ptype: analysis type
                                t: thickness
     
        D                       constitutive matrix
    
        eq = [[bx],               bx: body force x-dir
              [by]]               by: body force y-dir
              
    Returns:
    
        Ke                      element stiffness matrix (6 x 6)
        fe                      equivalent nodal forces (6 x 1) (if eq is given)

    """

    ptype,t = ep
    
    bx = 0.0
    by = 0.0
    
    if not eq is None:
        bx = eq[0]
        by = eq[1]
        
    C = mat([
        [1, ex[0], ey[0], 0,     0,     0], 
        [0,     0,     0, 1, ex[0], ey[0]],
        [1, ex[1], ey[1], 0,     0,     0],
        [0,     0,     0, 1, ex[1], ey[1]],
        [1, ex[2], ey[2], 0,     0,     0],
        [0,     0,     0, 1, ex[2], ey[2]]
        ])
    
    A = 0.5*linalg.det(mat([
        [1, ex[0], ey[0]],
        [1, ex[1], ey[1]],
        [1, ex[2], ey[2]]
        ]))
    
    # --------- plane stress --------------------------------------
    
    if ptype == 1:
        B = mat([
            [0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 1],
            [0, 0, 1, 0, 1, 0]
            ])*linalg.inv(C)
        
        colD = D.shape[1]
        
        if colD > 3:
            Cm = linalg.inv(D)
            Dm = linalg.inv(Cm[ix_((0,1,3),(0,1,3))])
        else:
            Dm = D
            
        Ke = B.T*Dm*B*A*t
        fe = A/3*mat([bx,by,bx,by,bx,by]).T*t
        
        if eq is None:
            return Ke
        else:
            return Ke,fe.T
       
    #--------- plane strain --------------------------------------       
    
    elif ptype == 2:
        B = mat([
            [0, 1, 0, 0, 0, 0,],
            [0, 0, 0, 0, 0, 1,],
            [0, 0, 1, 0, 1, 0,]
            ])*linalg.inv(C)

        colD = D.shape[1]
        
        if colD > 3:
            Dm = D[ix_((0,1,3),(0,1,3))]
        else:
            Dm = D

        Ke = B.T*Dm*B*A*t
        fe = A/3*mat([bx,by,bx,by,bx,by]).T*t
        
        if eq == None:
            return Ke
        else:
            return Ke,fe.T

    else:
        print("Error ! Check first argument, ptype=1 or 2 allowed")
        if eq == None:
            return None
        else:
            return None,None

def plants(ex,ey,ep,D,ed):
    """
    Calculate element normal and shear stress for a
    triangular plane stress or plane strain element.
    
    INPUT:  ex = [x1 x2 x3]         element coordinates
           ey = [y1 y2 y3]
    
           ep = [ptype t ]         ptype: analysis type
                                   t: thickness
    
           D                       constitutive matrix
    
           ed =[u1 u2 ...u6        element displacement vector
                ......     ]       one row for each element
    
    OUTPUT: es = [ sigx sigy [sigz] tauxy   element stress matrix
                 ......                 ]  one row for each element
    
           et = [ epsx epsy [epsz] gamxy   element strain matrix
                 ......                 ]  one row for each element
    """

    ptype=ep[0]
    
    if ndim(ex) == 1:
        ex = array([ex])
    if ndim(ey) == 1:
        ey = array([ey])
    if ndim(ed) == 1:
        ed = array([ed])

    rowed=ed.shape[0]
    rowex=ex.shape[0]
    
    # --------- plane stress --------------------------------------
    
    if ptype==1:
        
        colD = D.shape[1]

        if colD>3:
            Cm = linalg.inv(D)
            Dm = linalg.inv(Cm[ix_((0,1,3),(0,1,3))])
        else:
            Dm = D
            
        incie=0

        if rowex==1:
            incie=0
        else:
            incie=1
      
        et=zeros([rowed,colD])
        es=zeros([rowed,colD])
        
        ie=0
        
        for i in range(rowed):
            C = matrix(
                [[1, ex[ie,0], ey[ie,0], 0, 0, 0 ], 
                 [0, 0, 0, 1, ex[ie,0], ey[ie,0] ],
                 [1, ex[ie,1], ey[ie,1], 0, 0, 0 ],
                 [0, 0, 0, 1, ex[ie,1], ey[ie,1] ],
                 [1, ex[ie,2], ey[ie,2], 0, 0, 0 ],
                 [0, 0, 0, 1, ex[ie,2], ey[ie,2] ]]
                )
            
            B = matrix([
                [0,1,0,0,0,0],
                [0,0,0,0,0,1],
                [0,0,1,0,1,0]])*linalg.inv(C)
            
            ee=B*asmatrix(ed[ie,:]).T

            if colD>3:
                ss=zeros([colD,1])
                ss[[0,1,3]]=Dm*ee
                ee=Cm*ss
            else:
                ss=Dm*ee
    
            et[ie,:] = ee.T
            es[ie,:] = ss.T
    
            ie = ie + incie
            
        return es, et
    
def plantf(ex,ey,ep,es):
    """
    Compute internal element force vector in a triangular element
    in plane stress or plane strain. 

    Parameters:

        ex = [x1,x2,x3]                 node coordinates
        ey = [y1,y2,y3]

        ep = [ptype,t]                  ptype: analysis type
                                        t: thickness

        es = [[sigx,sigy,[sigz],tauxy]  element stress matrix
              [  ......              ]] one row for each element

    OUTPUT:

        fe = [[f1],[f2],...,[f8]]       internal force vector

    """

    ptype,t = ep

    colD = es.shape[1]

    #--------- plane stress --------------------------------------

    if ptype == 1:

        C = mat([
            [ 1, ex[0], ey[0], 0, 0,     0    ],
            [ 0, 0,     0,     1, ex[0], ey[0]],
            [ 1, ex[1], ey[1], 0, 0,     0    ],
            [ 0, 0,     0,     1, ex[1], ey[1]],
            [ 1, ex[2], ey[2], 0, 0,     0    ],
            [ 0, 0,     0,     1, ex[2], ey[2]]
        ])

        A = 0.5*linalg.det(mat([
            [ 1, ex[0], ey[0]],
            [ 1, ex[1], ey[1]],
            [ 1, ex[2], ey[2]]
        ]))

        B = mat([
            [ 0, 1, 0, 0, 0, 0],
            [ 0, 0, 0, 0, 0, 1],
            [ 0, 0, 1, 0, 1, 0]
        ])*linalg.inv(C)
    
        if colD > 3:
            stress = asmatrix(es[ix_((0,1,3))])
        else:
            stress = asmatrix(es)

        ef = (A*t*B.T*stress.T).T
        
        return reshape(asarray(ef),6)

    #--------- plane strain --------------------------------------
    
    elif ptype == 2:

        C = mat([
            [ 1, ex[0], ey[0], 0, 0,     0    ],
            [ 0, 0,     0,     1, ex[0], ey[0]],
            [ 1, ex[1], ey[1], 0, 0,     0    ],
            [ 0, 0,     0,     1, ex[1], ey[1]],
            [ 1, ex[2], ey[2], 0, 0,     0    ],
            [ 0, 0,     0,     1, ex[2], ey[2]]
        ])

        A = 0.5*linalg.det(mat([
            [ 1, ex[0], ey[0]],
            [ 1, ex[1], ey[1]],
            [ 1, ex[2], ey[2]]
        ]))

        B = mat([
            [ 0, 1, 0, 0, 0, 0],
            [ 0, 0, 0, 0, 0, 1],
            [ 0, 0, 1, 0, 1, 0]
        ])*linalg.inv(C)

        if colD > 3:
            stress = asmatrix(es[ix_((1,2,4))])
        else:
            stress = asmatrix(es)

        ef = (A*t*B.T*stress.T).T
        
        return reshape(asarray(ef),6)
  
    else:
        print("Error ! Check first argument, ptype=1 or 2 allowed")
        return None

def platre(ex,ey,ep,D,eq=None):
    """
    Calculate the stiffness matrix for a rectangular plate element.
    NOTE! Element sides must be parallel to the coordinate axis.
    
    Parameters:

        ex = [x1,x2,x3,x4]          element coordinates
        ey = [y1,y2,y3,y4]

        ep = [t]                    thicknes

        D                           constitutive matrix for
                                    plane stress         

        eq = [qz]                   load/unit area
    Returns:

        Ke                          element stiffness matrix (12 x 12)
        fe                          equivalent nodal forces (12 x 1)

    """
    Lx = (ex[2]-ex[0]).astype(float)
    Ly = (ey[2]-ey[0]).astype(float)
    t = ep[0]

    D = t**3/12.*D

    A1 = Ly/(Lx**3)
    A2 = Lx/(Ly**3)
    A3 = 1/Lx/Ly
    A4 = Ly/(Lx**2)
    A5 = Lx/(Ly**2)
    A6 = 1/Lx
    A7 = 1/Ly
    A8 = Ly/Lx
    A9 = Lx/Ly

    C1 = 4*A1*D[0,0]+4*A2*D[1,1]+2*A3*D[0,1]+5.6*A3*D[2,2]
    C2 = -4*A1*D[0,0]+2*A2*D[1,1]-2*A3*D[0,1]-5.6*A3*D[2,2]
    C3 = 2*A1*D[0,0]-4*A2*D[1,1]-2*A3*D[0,1]-5.6*A3*D[2,2]
    C4 = -2*A1*D[0,0]-2*A2*D[1,1]+2*A3*D[0,1]+5.6*A3*D[2,2]
    C5 = 2*A5*D[1,1]+A6*D[0,1]+0.4*A6*D[2,2]
    C6 = 2*A4*D[0,0]+A7*D[0,1]+0.4*A7*D[2,2]
    
    C7 = 2*A5*D[1,1]+0.4*A6*D[2,2]
    C8 = 2*A4*D[0,0]+0.4*A7*D[2,2]
    C9 = A5*D[1,1]-A6*D[0,1]-0.4*A6*D[2,2]
    C10 = A4*D[0,0]-A7*D[0,1]-0.4*A7*D[2,2]
    C11 = A5*D[1,1]-0.4*A6*D[2,2]
    C12 = A4*D[0,0]-0.4*A7*D[2,2]
    
    C13 = 4/3.*A9*D[1,1]+8/15.*A8*D[2,2]
    C14 = 4/3.*A8*D[0,0]+8/15.*A9*D[2,2]
    C15 = 2/3.*A9*D[1,1]-8/15.*A8*D[2,2]
    C16 = 2/3.*A8*D[0,0]-8/15.*A9*D[2,2]
    C17 = 2/3.*A9*D[1,1]-2/15.*A8*D[2,2]
    C18 = 2/3.*A8*D[0,0]-2/15.*A9*D[2,2]
    C19 = 1/3.*A9*D[1,1]+2/15.*A8*D[2,2]
    C20 = 1/3.*A8*D[0,0]+2/15.*A9*D[2,2]
    C21 = D[0,1]

    Keq = mat(zeros((12,12)))
    Keq[0,0:13] = C1,C5,-C6,C2,C9,-C8,C4,C11,-C12,C3,C7,-C10
    Keq[1,1:13] = C13,-C21,C9,C15,0,-C11,C19,0,-C7,C17,0
    Keq[2,2:13] = C14,C8,0,C18,C12,0,C20,-C10,0,C16
    Keq[3,3:13] = C1,C5,C6,C3,C7,C10,C4,C11,C12
    Keq[4,4:13] = C13,C21,-C7,C17,0,-C11,C19,0
    Keq[5,5:13] = C14,C10,0,C16,-C12,0,C20
    Keq[6,6:13] = C1,-C5,C6,C2,-C9,C8
    Keq[7,7:13] = C13,-C21,-C9,C15,0
    Keq[8,8:13] = C14,-C8,0,C18
    Keq[9,9:13] = C1,-C5,-C6
    Keq[10,10:13] = C13,C21
    Keq[11,11] = C14
    Keq = Keq.T+Keq-diag(diag(Keq))

    if eq != None:
        q = eq
        R1 = q*Lx*Ly/4
        R2 = q*Lx*Ly**2/24
        R3 = q*Ly*Lx**2/24

        feq = mat([R1,R2,-R3,R1,R2,R3,R1,-R2,R3,R1,-R2,-R3])

    if eq != None:
        return Keq,feq
    else:
        return Keq

        
def planqe(ex,ey,ep,D,eq=None):
    """
    Calculate the stiffness matrix for a quadrilateral
    plane stress or plane strain element.

    Parameters:
        ex=[x1 x2 x3 x4]    element coordinates
        ey=[y1 y2 y3 y4]
                                
        ep = [ptype, t]     ptype: analysis type
                            t: element thickness 

        D                   constitutive matrix

        eq = [bx;           bx: body force in x direction
              by]           by: body force in y direction

    OUTPUT: Ke :  element stiffness matrix (8 x 8)
            fe : equivalent nodal forces (row array)
    """
    K=zeros((10,10))
    f=zeros((10,1))
    
    xm=sum(ex)/4.
    ym=sum(ey)/4.
    
    b1 = eq if eq is not None else array([[0],[0]])
    
    ke1, fe1 = plante(array([ex[0], ex[1], xm]), array([ey[0], ey[1], ym]), ep, D, b1)
    K, f = assem(array([1, 2, 3, 4, 9, 10]), K, ke1, f, fe1)
    ke1, fe1 = plante(array([ex[1], ex[2], xm]), array([ey[1], ey[2], ym]), ep, D, b1)
    K, f = assem(array([3, 4, 5, 6, 9, 10]), K, ke1, f, fe1)
    ke1, fe1 = plante(array([ex[2], ex[3], xm]), array([ey[2], ey[3], ym]), ep, D, b1)
    K, f = assem(array([5, 6, 7, 8, 9, 10]), K, ke1, f, fe1)
    ke1, fe1 = plante(array([ex[3], ex[0], xm]), array([ey[3], ey[0], ym]), ep, D, b1)
    K, f = assem(array([7, 8, 1, 2, 9, 10]), K, ke1, f, fe1)
    Ke, fe = statcon(K, f, array([[9],[10]]))
    
    if eq == None:
        return Ke
    else:
        return Ke,fe
        

def planqs(ex,ey,ep,D,ed,eq=None):
    """
    Calculate element normal and shear stress for a quadrilateral 
    plane stress or plane strain element.
    
    Parameters:
            ex = [x1 x2 x3 x4]      element coordinates
            ey = [y1 y2 y3 y4]
    
            ep = [ptype, t]         ptype: analysis type
                                    t:  thickness
                                   
            D                       constitutive matrix
    
            ed = [u1 u2 ..u8]       element displacement vector
    
            eq = [[bx]               bx: body force in x direction
                  [by]]              by: body force in y direction
    
    OUTPUT: es = [ sigx sigy (sigz) tauxy]    element stress array
            et = [ epsx epsy (epsz) gamxy]    element strain array
    """
    
    if ex.shape != (4,) or ey.shape != (4,) or ed.shape != (8,):
        raise ValueError('Error ! PLANQS: only one element at the time (ex, ey, ed must be a row arrays)')

    K = zeros((10,10))
    f = zeros((10,1))
    
    xm = sum(ex)/4.
    ym = sum(ey)/4.
        
    b1 = eq if eq is not None else array([[0],[0]])
    
    ex1 = array([ex[0], ex[1], xm]) 
    ey1 = array([ey[0], ey[1], ym])
    ex2 = array([ex[1], ex[2], xm]) 
    ey2 = array([ey[1], ey[2], ym])
    ex3 = array([ex[2], ex[3], xm]) 
    ey3 = array([ey[2], ey[3], ym])
    ex4 = array([ex[3], ex[0], xm]) 
    ey4 = array([ey[3], ey[0], ym])
    
    ke1, fe1 = plante(ex1, ey1, ep, D, b1)
    K, f = assem(array([1, 2, 3, 4, 9, 10]), K, ke1, f, fe1)
    ke1,fe1 = plante(ex2, ey2, ep, D, b1)
    K, f = assem(array([3, 4, 5, 6, 9, 10]), K, ke1, f, fe1)
    ke1, fe1 = plante(ex3, ey3, ep, D, b1)
    K, f = assem(array([5, 6, 7, 8, 9, 10]), K, ke1, f, fe1)
    ke1, fe1 = plante(ex4, ey4, ep, D, b1)
    K, f = assem(array([7, 8, 1, 2, 9, 10]), K, ke1, f, fe1)
    
    A1 = 0.5 * linalg.det( hstack([ones((3,1)), mat(ex1).T, mat(ey1).T]) )
    A2 = 0.5 * linalg.det( hstack([ones((3,1)), mat(ex2).T, mat(ey2).T]) )
    A3 = 0.5 * linalg.det( hstack([ones((3,1)), mat(ex3).T, mat(ey3).T]) )
    A4 = 0.5 * linalg.det( hstack([ones((3,1)), mat(ex4).T, mat(ey4).T]) )
    Atot = A1+A2+A3+A4;
    
    
    a, _ = solveq(K, f, array(range(1,9)), ed)
        
#    ni = ed.shape[0]
#    a = mat(empty((10,ni)))
#    for i in range(ni):
#        a[:,i] = solveq(K, f, array(range(1,9)), ed[i,:])[0]
#        #a = hstack([a, solveq(K, f, hstack([matrix(range(1,9)).T, ed[i,:].T]) ) ])
    
    s1, t1 = plants(ex1, ey1, ep, D, hstack([a[[0, 1, 2, 3, 8, 9], :].T]) );
    s2, t2 = plants(ex2, ey2, ep, D, hstack([a[[2, 3, 4, 5, 8, 9], :].T]) );
    s3, t3 = plants(ex3, ey3, ep, D, hstack([a[[4, 5, 6, 7, 8, 9], :].T]) );
    s4, t4 = plants(ex4, ey4, ep, D, hstack([a[[6, 7, 0, 1, 8, 9], :].T]) );
    
    es = (s1*A1+s2*A2+s3*A3+s4*A4)/Atot;
    et = (t1*A1+t2*A2+t3*A3+t4*A4)/Atot;
    
    return es[0], et[0] #[0] because these are 1-by-3 arrays and we want row arrays out.


def plani4e(ex,ey,ep,D,eq=None):
    """
    Calculate the stiffness matrix for a 4 node isoparametric
    element in plane strain or plane stress.
    
    Parameters:
        ex = [x1 ...   x4]  element coordinates. Row array
        ey = [y1 ...   y4]
                                
        ep =[ptype, t, ir]  ptype: analysis type
                            t : thickness
                            ir: integration rule
    
        D                   constitutive matrix
    
        eq = [bx; by]       bx: body force in x direction
                            by: body force in y direction
                                Any array with 2 elements acceptable
    
    Returns:
        Ke : element stiffness matrix (8 x 8)
        fe : equivalent nodal forces (8 x 1)
    """
    ptype=ep[0] 
    t=ep[1]  
    ir=ep[2]  
    ngp=ir*ir
    if eq == None:
        q = zeros((2,1))
    else:
        q = reshape(eq, (2,1))
#--------- gauss points --------------------------------------    
    if ir == 1:
        g1 = 0.0
        w1 = 2.0
        gp = mat([g1,g1])
        w = mat([w1,w1])
    elif ir == 2:
        g1 = 0.577350269189626
        w1 = 1
        gp = mat([
            [-g1,-g1],
            [ g1,-g1],
            [-g1, g1],
            [ g1, g1]])
        w = mat([
            [ w1, w1],
            [ w1, w1],
            [ w1, w1],
            [ w1, w1]])
    elif ir == 3:
        g1 = 0.774596669241483
        g2 = 0.
        w1 = 0.555555555555555
        w2 = 0.888888888888888
        gp = mat([
            [-g1,-g1],
            [-g2,-g1],
            [ g1,-g1],
            [-g1, g2],
            [ g2, g2],
            [ g1, g2],
            [-g1, g1],
            [ g2, g1],
            [ g1, g1]])
        w = mat([
            [ w1, w1],
            [ w2, w1],
            [ w1, w1],
            [ w1, w2],
            [ w2, w2],
            [ w1, w2],
            [ w1, w1],
            [ w2, w1],
            [ w1, w1]])
    else:
        print("Used number of integration points not implemented")
    wp = multiply(w[:,0],w[:,1])
    xsi = gp[:,0]
    eta = gp[:,1]
    r2 = ngp*2
    # Shape Functions
    N = multiply((1-xsi),(1-eta))/4.
    N = append(N,multiply((1+xsi),(1-eta))/4.,axis=1)
    N = append(N,multiply((1+xsi),(1+eta))/4.,axis=1)
    N = append(N,multiply((1-xsi),(1+eta))/4.,axis=1)
    
    dNr = mat(zeros((r2,4)))
    dNr[0:r2:2,0] = -(1-eta)/4.
    dNr[0:r2:2,1] = (1-eta)/4.
    dNr[0:r2:2,2] = (1+eta)/4.
    dNr[0:r2:2,3] = -(1+eta)/4.
    dNr[1:r2+1:2,0] = -(1-xsi)/4.
    dNr[1:r2+1:2,1] = -(1+xsi)/4.
    dNr[1:r2+1:2,2] = (1+xsi)/4.
    dNr[1:r2+1:2,3] = (1-xsi)/4.

#
    Ke1 = mat(zeros((8,8)))
    fe1 = mat(zeros((8,1)))
    JT = dNr*mat([ex,ey]).T 
    # --------- plane stress --------------------------------------
    if ptype==1:
        colD=shape(D)[0]
        if colD>3:
            Cm=linalg.inv(D)
            Dm=linalg.inv(Cm[ ix_([0,1,3],[0,1,3]) ])
        else:
            Dm=D          
#
        B=matrix(zeros((3,8)))
        N2=matrix(zeros((2,8)))
        for i in range(ngp):
            indx = array([2*(i+1)-1,2*(i+1)])
            detJ = linalg.det(JT[indx-1,:])
            if detJ < 10*finfo(float).eps:
                print("Jacobideterminant equal or less than zero!")
            JTinv = linalg.inv(JT[indx-1,:])  
            dNx=JTinv*dNr[indx-1,:]
#   
            index_array_even=array([0,2,4,6])
            index_array_odd=array([1,3,5,7])
#    
            counter=0    
            for index in index_array_even:      
                B[0,index] = dNx[0,counter]
                B[2,index] = dNx[1,counter]
                N2[0,index]=N[i,counter]
                counter=counter+1
#
            counter=0    
            for index in index_array_odd:
                B[1,index]   = dNx[1,counter]
                B[2,index]   = dNx[0,counter]
                N2[1,index]  =N[i,counter]
                counter=counter+1
#   
            Ke1 = Ke1+B.T*Dm*B*detJ*asscalar(wp[i])*t
            fe1 = fe1 + N2.T * q * detJ * asscalar(wp[i]) * t

        return Ke1,fe1
#--------- plane strain --------------------------------------
    elif ptype==2:
#      
        colD=shape(D)[0]
        if colD>3:
            Dm = D[ix_([0,1,3],[0,1,3])]
        else:
            Dm = D
#
        B=matrix(zeros((3,8)))
        N2=matrix(zeros((2,8)))
        for i in range(ngp):
            indx = array([2*(i+1)-1,2*(i+1)])
            detJ = linalg.det(JT[indx-1,:])
            if detJ < 10*finfo(float).eps:
                print("Jacobideterminant equal or less than zero!")
            JTinv = linalg.inv(JT[indx-1,:])  
            dNx=JTinv*dNr[indx-1,:]
#   
            index_array_even=array([0,2,4,6])
            index_array_odd=array([1,3,5,7])
#    
            counter=0    
            for index in index_array_even:
#
                B[0,index] = dNx[0,counter]
                B[2,index] = dNx[1,counter]
                N2[0,index]=N[i,counter]
#
                counter=counter+1
#
            counter=0    
            for index in index_array_odd:
                B[1,index]   = dNx[1,counter]
                B[2,index]   = dNx[0,counter]
                N2[1,index]  =N[i,counter]
                counter=counter+1
#   
            Ke1 = Ke1 + B.T * Dm * B * detJ * asscalar(wp[i]) * t
            fe1 = fe1+N2.T*q*detJ*asscalar(wp[i])*t
        return Ke1,fe1
    else:
        print("Error ! Check first argument, ptype=1 or 2 allowed")
        
        
def assem(edof,K,Ke,f=None,fe=None):
    """
    Assemble element matrices Ke ( and fe ) into the global
    stiffness matrix K ( and the global force vector f )
    according to the topology matrix edof.
    
    Parameters:
    
        edof        dof topology array
        K           the global stiffness matrix
        Ke          element stiffness matrix
        f           the global force vector
        fe          element force vector
        
    Output parameters:
    
        K           the new global stiffness matrix
        f           the new global force vector
        fe          element force vector
    
    """
    
    if rank(edof) == 1:
        idx = edof-1
        K[ix_(idx,idx)] = K[ix_(idx,idx)] + Ke
        if (not f is None) and (not fe is None):
            f[ix_(idx)] = f[ix_(idx)] + fe
    else:
        for row in edof:
            idx = row-1
            K[ix_(idx,idx)] = K[ix_(idx,idx)] + Ke
            if (not f is None) and (not fe is None):
                f[ix_(idx)] = f[ix_(idx)] + fe
            
    if f is None:
        return K
    else:
        return K,f
            
def solveq(K,f,bcPrescr,bcVal=None):
    """
    Solve static FE-equations considering boundary conditions.
    
    Parameters:
    
        K           global stiffness matrix, dim(K)= nd x nd
        f           global load vector, dim(f)= nd x 1
    
        bcPrescr    1-dim integer array containing prescribed dofs.
        bcVal       1-dim float array containing prescribed values.
                    If not given all prescribed dofs are assumed 0.
        
    Returns:
    
        a           solution including boundary values
        Q           reaction force vector
                    dim(a)=dim(Q)= nd x 1, nd : number of dof's
    
    """    
    
    nDofs = K.shape[0]
    nPdofs = bcPrescr.shape[0]
    
    if bcVal is None:
        bcVal = zeros([nPdofs],'d')
    
    bc = ones(nDofs, 'bool')    
    bcDofs = arange(nDofs)
    
    bc[ix_(bcPrescr-1)] = False
    bcDofs = bcDofs[bc]
    
    fsys = f[bcDofs]-K[ix_((bcDofs),(bcPrescr-1))]*asmatrix(bcVal).reshape(nPdofs,1)
    asys = linalg.solve(K[ix_((bcDofs),(bcDofs))], fsys);
    
    a = zeros([nDofs,1])
    a[ix_(bcPrescr-1)] = asmatrix(bcVal).reshape(nPdofs,1)
    a[ix_(bcDofs)] = asys
    
    Q=K*asmatrix(a)-f
    
    return (asmatrix(a),Q)
    
def extractEldisp(edof,a):
    """
    Extract element displacements from the global displacement
    vector according to the topology matrix edof.
    
    Parameters:
    
        a           the global displacement vector
        edof        dof topology array
    
    Returns:
    
        ed:     element displacement array
    
    """

    ed = None
    
    if rank(edof)==1:
        nDofs = len(edof)
        ed = zeros([nDofs])
        idx = edof-1
        ed[:] = a[ix_(idx)].T
    else:
        nElements = edof.shape[0]
        nDofs = edof.shape[1]
        ed = zeros([nElements,nDofs])
        i=0
        for row in edof:
            idx = row-1
            ed[i,:]=a[ix_(idx)].T
            i+=1
        
    return ed

def statcon(K,f,cd):
    """
    Condensation of static FE-equations according to the vector cd.

    Parameters:
    
        K                       global stiffness matrix, dim(K) = nd x nd
        f                       global load vector, dim(f)= nd x 1

        cd                      vector containing dof's to be eliminated
                                dim(cd)= nc x 1, nc: number of condensed dof's
    Returns:
    
        K1                      condensed stiffness matrix,
                                dim(K1)= (nd-nc) x (nd-nc)
        f1                      condensed load vector, dim(f1)= (nd-nc) x 1
    """
    nd,nd = shape(K)
    cd = (cd-1).flatten()
  
    aindx = arange(nd)
    aindx = delete(aindx,cd,0)
    bindx = cd

    Kaa = mat(K[ix_(aindx,aindx)])
    Kab = mat(K[ix_(aindx,bindx)])
    Kbb = mat(K[ix_(bindx,bindx)])

    fa = mat(f[aindx])
    fb = mat(f[bindx])
    
    K1 = Kaa-Kab*Kbb.I*Kab.T
    f1 = fa-Kab*Kbb.I*fb
    
    return K1,f1

def c_mul(a, b):
    return eval(hex((long(a) * b) & 0xFFFFFFFF)[:-1])

def dofHash(dof):
    if len(dof)==1:
        return dof[0]
    value = 0x345678
    for item in dof:
        value = c_mul(1000003, value) ^ hash(item)
    value = value ^ len(dof)
    if value == -1:
        value = -2
    return value

def createdofs(nCoords,nDof):
    """
    Create dof array [nCoords x nDof]
    """
    return arange(nCoords*nDof).reshape(nCoords,nDof)+1

def coordxtr(edof,coords,dofs):
    """
    Create element coordinate matrices ex, ey, ez from edof
    coord and dofs matrices.
    
    Parameters:
    
        edof            [nel x (nen * nnd)], nnd = number of node dofs
        coords          [ncoords x ndims],   ndims = node dimensions
        dofs            [ncoords x nnd]
        
    Returns:
    
        ex              if ndims = 1
        ex, ey          if ndims = 2
        ex, ey, ez      if ndims = 3
    """
    
    # Create dictionary with dof indices
    
    dofDict = {}
    nDofs = size(dofs,1)
    nElements = size(edof,0)
    nDimensions = size(coords,1)
    nElementDofs = size(edof,1)
    
    nElementNodes = int(nElementDofs/nDofs)
    
    idx = 0
    for dof in dofs:
        dofDict[dofHash(dof)] = idx
        idx += 1
              
    # Loop over edof and extract element coords
    
    ex = zeros((nElements,nElementNodes))
    ey = zeros((nElements,nElementNodes))
    ez = zeros((nElements,nElementNodes))
    
    elementIdx = 0
    for etopo in edof:
        for i in range(nElementNodes):
            i0 = i*nDofs
            i1 = i*nDofs+nDofs-1
            dof = []
            if i0==i1:
                dof = [etopo[i*nDofs]]
            else:
                dof = etopo[i*nDofs:(i*nDofs+nDofs)]
            
            nodeCoord = coords[dofDict[dofHash(dof)]]
            
            if nDimensions>=1:
                ex[elementIdx,i] = nodeCoord[0]
            if nDimensions>=2:
                ey[elementIdx,i] = nodeCoord[1]
            if nDimensions>=3:
                ez[elementIdx,i] = nodeCoord[2]
            
        elementIdx += 1
        
    if nDimensions==1:
        return ex
    
    if nDimensions==2:
        return ex, ey
    
    if nDimensions==3:
        return ex, ey, ez

def hooke(ptype,E,v):
    """
    Calculate the material matrix for a linear
    elastic and isotropic material.
    
    Parameters:
    
        ptype=  1:  plane stress
                2:  plane strain
                3:  axisymmetry
                4:  three dimensional
    
        E           Young's modulus
        v           Poissons const.
        
    Returns:
    
        D           material matrix
    
    """
   
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
        print("ptype not supported.")
        
    return D

def effmises(es,ptype):
    """
    Calculate effective von mises stresses.
    
    Parameters:
        
        es
    
        ptype=  1:  plane stress
                2:  plane strain
                3:  axisymmetry
                4:  three dimensional
    
       es = [[sigx,sigy,[sigz],tauxy]  element stress matrix
              [  ......              ]] one row for each element
              
    Returns:
    
        eseff  = [eseff_0 .. eseff_nel-1]
    
    """
    
    nel = size(es,0)
    escomps = size(es, 1)
    
    eseff = zeros([nel])
     
    if ptype == 1:
        sigxx = es[:,0]
        sigyy = es[:,1]
        sigxy = es[:,2]
        eseff = sqrt(sigxx*sigxx+sigyy*sigyy-sigxx*sigyy+3*sigxy*sigxy)
        return eseff
    
def stress2nodal(eseff, edof):
    """
    Convert element effective stresses to nodal effective
    stresses.
    
    Parameters:
        
        eseff  = [eseff_0 .. eseff_nel-1] 
        edof   = [dof topology array]
    
    Returns:
    
        ev:     element value array [[ev_0_0 ev_0_1 ev_0_nen-1 ]
                                      ..
                                      ev_nel-1_0 ev_nel-1_1 ev_nel-1_nen-1]
                      
    """
    
    values = zeros(edof.max())
    elnodes = size(edof,1) / 2
    
    for etopo, eleseff in zip(edof, eseff):
        values[etopo-1] = values[etopo-1] + eleseff / elnodes
        
    evtemp = extractEldisp(edof,values)
    ev = evtemp[:,range(0,elnodes*2,2)]
                     
    return ev