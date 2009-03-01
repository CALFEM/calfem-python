from numpy import *

haveMatplotLib = True

try:
    from matplotlib.pyplot import *
except:
    haveMatplotLib = False

def spring1e(ep):
    """-------------------------------------------------------------
     PURPOSE
      Compute element stiffness matrix for spring element.
    
     INPUT:  ep = k;       spring stiffness or analog quantity
    
     OUTPUT: Ke:           stiffness matrix, dim(Ke)= 2 x 2
    -------------------------------------------------------------"""
    k = ep
    return mat([[k,-k],[-k,k]],'d')

def spring1s(ep,ed):
    """-------------------------------------------------------------
     PURPOSE
      Compute element force in spring element (spring1e).
    
     INPUT:  ep = k        spring stiffness or analog quantity
             ed = [u1 u2]  element displacements
                           u1, u2: nodal displacements
                                
    
     OUTPUT: es  = [N]       element force
    -------------------------------------------------------------"""
    k = ep
    return k*(ed[1]-ed[0]);   

def bar1e(ep):
    """-------------------------------------------------------------
     PURPOSE
      Compute element stiffness matrix for spring element.
    
     INPUT:  ep = k;       spring stiffness or analog quantity
    
     OUTPUT: Ke:           stiffness matrix, dim(Ke)= 2 x 2
    -------------------------------------------------------------"""
    k = ep
    return mat([[k,-k],[-k,k]],'d')

def bar1s(ep,ed):
    """-------------------------------------------------------------
     PURPOSE
      Compute element force in spring element (spring1e).
    
     INPUT:  ep = k        spring stiffness or analog quantity
             ed = [u1 u2]  element displacements
                           u1, u2: nodal displacements
                                
    
     OUTPUT: es  = [N]       element force
    -------------------------------------------------------------"""
    k = ep
    return k*(ed[1]-ed[0]);   

def bar2e(ex,ey,ep):
    """----------------------------------------------------------------------
    PURPOSE
     Compute the element stiffness matrix for two dimensional bar element.
    
    INPUT:  ex = [x1 x2];
            ey = [y1 y2];      element node coordinates
    
            ep = [E A]         E: Young's modulus
                               A: Cross section area
    
    OUTPUT: Ke : stiffness matrix, dim(Ke)= 4 x 4
    ----------------------------------------------------------------------"""
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

def bar2s(ex,ey,ep,ed):
    """-------------------------------------------------------------
     PURPOSE
      Compute normal force in two dimensional bar element.
    
     INPUT:  ex = [x1 x2]
             ey = [y1 y2]        element coordinates
    
             ep = [E A]          E : Young's modulus
                                 A : Cross section area
    
             ed : [u1 u2 u3 u4]  element displacement vector
    
     OUTPUT: es = [N]            element force 
    -------------------------------------------------------------"""
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

def beam2e(ex,ey,ep,eq=None):
    """---------------------------------------------------------------------
        PURPOSE
         Compute the stiffness matrix for a two dimensional beam element. 
     
        INPUT:  ex = [x1 x2]
                ey = [y1 y2]       element node coordinates
    
                ep = [E A I]       element properties
                                      E: Young's modulus
                                      A: Cross section area
                                      I: Moment of inertia
    
                eq = [qx qy]       distributed loads, local directions
     
        OUTPUT: Ke : element stiffness matrix (6 x 6)
    
                fe : element load vector (6 x 1)
    --------------------------------------------------------------------"""

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
    """---------------------------------------------------------------------
    PURPOSE
      Compute section forces in two dimensional beam element (beam2e). 
 
    INPUT:  ex = [x1 x2]
            ey = [y1 y2]     element node coordinates

            ep = [E A I]     element properties,
                              E:  Young's modulus
                              A:  cross section area
                              I:  moment of inertia

            ed = [u1 ... u6] element displacements

            eq = [qx qy]     distributed loads, local directions 

            n : number of evaluation points ( default=2 )
          
    OUTPUT: es = [ N1 V1 M1 ;  section forces, local directions, in 
                   N2 V2 M2 ;  n points along the beam, dim(es)= n x 3
                   .........]  
           
            edi = [ u1 v1 ;    element displacements, local directions,
                    u2 v2 ;    in n points along the beam, dim(es)= n x 2
                   .......]    

            eci = [ x1  ;      local x-coordinates of the evaluation 
                    x2 ;       points, (x1=0 and xn=L)
                    ...]
    -------------------------------------------------------------------------"""
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
    
def flw2te(ex,ey,ep,D,eq=None):
    """Ke=flw2te(ex,ey,ep,D)
    [Ke,fe]=flw2te(ex,ey,ep,D,eq)
    -------------------------------------------------------------
     PURPOSE
      Compute element stiffness (conductivity) matrix for a 
      triangular field element.
    
     INPUT:  ex = [x1 x2 x3]
             ey = [y1 y2 y3]      element coordinates
    
             ep = [t]  	       element thickness  	 
                                 
             D = [kxx kxy;
                  kyx kyy]        constitutive matrix
    
             eq                   heat supply per unit volume
    
     OUTPUT: Ke :  element 'stiffness' matrix (3 x 3)
    
             fe :  element load vector (3 x 1)
    -------------------------------------------------------------"""
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
    """[es,et]=flw2ts(ex,ey,D,ed)
    -------------------------------------------------------------
     PURPOSE
      Compute flows or corresponding quantities in the
      triangular field element.
    
     INPUT:  ex = [x1 x2 x3]
             ey = [y1 y2 y3]         element coordinates
                                 
             D = [kxx kxy;
                  kyx kyy]           constitutive matrix
    
             ed =[u1 u2 u3]          u1,u2,u3: nodal values
                  .. .. ..;
    
    
     OUTPUT: es=[ qx qy ] 
                  ... ..]            element flows
    
             et=[ gx gy ]
                  ... ..]            element gradients
    -------------------------------------------------------------"""

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

def assem(edof,K,Ke,f=None,fe=None):
    """-------------------------------------------------------------
    PURPOSE
     Assemble element matrices Ke ( and fe ) into the global
     stiffness matrix K ( and the global force vector f )
     according to the topology matrix edof.
    
    INPUT: edof:       dof topology array
           K :         the global stiffness matrix
           Ke:         element stiffness matrix
           f :         the global force vector
           fe:         element force vector
    
    OUTPUT: K :        the new global stiffness matrix
            f :        the new global force vector
    -------------------------------------------------------------"""
    
    if rank(edof) == 1:
        idx = edof-1
        K[ix_(idx,idx)] = K[ix_(idx,idx)] + Ke
        if (f!=None) and (fe!=None):
            f[ix_(idx)] = f[ix_(idx)] + fe
    else:
        for row in edof:
            idx = row-1
            K[ix_(idx,idx)] = K[ix_(idx,idx)] + Ke
            if (f!=None) and (fe!=None):
                f[ix_(idx)] = f[ix_(idx)] + fe
            
    return K
            
def solveq(K,f,bcPrescr,bcVal=None):
    """-------------------------------------------------------------
    PURPOSE
     Solve static FE-equations considering boundary conditions.
    
    INPUT:  K : global stiffness matrix, dim(K)= nd x nd
            f : global load vector, dim(f)= nd x 1
    
            bc :    1-dim integer array containing prescribed dofs.
                    
            bcVal:  1-dim float array containing prescribed values.
    
    OUTPUT:  a : solution including boundary values
             Q : reaction force vector
                 dim(a)=dim(Q)= nd x 1, nd : number of dof's
    -------------------------------------------------------------"""    
    
    nDofs = K.shape[0]
    nPdofs = bcPrescr.shape[0]
    
    if bcVal==None:
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
    
def extract(edof,a):
    """-------------------------------------------------------------
    PURPOSE
     Extract element displacements from the global displacement
     vector according to the topology matrix edof.
    
    INPUT:  a:      the global displacement vector
            edof:   dof topology array
    
    OUTPUT: ed:     element displacement array
    -------------------------------------------------------------"""

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

def c_mul(a, b):
    return eval(hex((long(a) * b) & 0xFFFFFFFFL)[:-1])

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
    return arange(nCoords*nDof).reshape(nCoords,nDof)+1

def coordxtr(edof,coords,dofs):
    
    # Create dictionary with dof indices
    
    dofDict = {}
    nDofs = size(dofs,1)
    nElements = size(edof,0)
    nDimensions = size(coords,1)
    nElementDofs = size(edof,1)
    
    nElementNodes = nElementDofs/nDofs
    
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
    """D=hooke(ptype,E,v)
    -------------------------------------------------------------
      PURPOSE
       Calculate the material matrix for a linear
       elastic and isotropic material.
    
     INPUT:  ptype=1:  plane stress
                   2:  plane strain
                   3:  axisymmetry
                   4:  three dimensional
    
              E : Young's modulus
              v : Poissons const.
    
     OUTPUT: D : material matrix
    -------------------------------------------------------------"""
   
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

def eldraw2(ex,ey,plotpar=None,elnum=None):
    if rank(ex)==1:
        nen = ex.shape[0]
        nel = 1
    else:
        nen = ex.shape[1]
        nel = ex.shape[0]
        
    if nen > 2:
        ex = hstack([ex,ex[:,0].reshape(nel,1)])
        ey = hstack([ey,ey[:,0].reshape(nel,1)])
        
    plot(ex.transpose(),ey.transpose(), color="blue")
    gca().set_aspect("equal")
    
def elmargin(scale=0.2):
    a = gca()
    xlim = a.get_xlim()
    ylim = a.get_ylim()
    xs = xlim[1]-xlim[0]
    ys = ylim[1]-ylim[0]
    a.set_xlim([xlim[0]-xs*scale,xlim[1]+xs*scale])
    a.set_ylim([ylim[0]-ys*scale,ylim[1]+ys*scale])
    
def scalfact2(ex,ey,ed,rat=0.2):
    """[sfac]=scalfact2(ex,ey,ed,rat)
    -------------------------------------------------------------
     PURPOSE 
       Determine scale factor for drawing computational results, such as 
       displacements, section forces or flux.
    
     INPUT
        ex,ey:  element node coordinates
                       
        ed:     element displacement matrix or section force matrix
    
        rat: relation between illustrated quantity and element size. 
        If not specified, 0.2 is used.
        
    -------------------------------------------------------------"""

    nen = -1
    if ex.shape != ey.shape:
        print "ex and ey shapes do not match."
        return 1.0
    
    dlmax = 0.
    edmax = 1.
    
    print rank(ex)

    if rank(ex)==1:
        nen = ex.shape[0]
        nel = 1
        dxmax=ex.T.max()-ex.T.min()
        dymax=ey.T.max()-ey.T.min()
        dlmax=max(dxmax,dymax);
        edmax=abs(ed).max();
    else:
        nen = ex.shape[1]
        nel = ex.shape[0]
        dxmax=ex.T.max()-ex.T.min()
        dymax=ey.T.max()-ey.T.min()
        dlmax=max(dxmax,dymax);
        edmax=abs(ed).max();
        
    k = rat
    return k*dlmax/edmax

def eldisp2(ex,ey,ed,rat=0.2):
    """[sfac]=scalfact2(ex,ey,ed,rat)
    -------------------------------------------------------------
     PURPOSE 
       Determine scale factor for drawing computational results, such as 
       displacements, section forces or flux.
    
     INPUT
        ex,ey:  element node coordinates
                       
        ed:     element displacement matrix or section force matrix
    
        rat: relation between illustrated quantity and element size. 
        If not specified, 0.2 is used.
        
    -------------------------------------------------------------"""

    nen = -1
    if ex.shape != ey.shape:
        print "ex and ey shapes do not match."
        return 1.0
    
    dlmax = 0.
    edmax = 1.
    
    print rank(ex)

    if rank(ex)==1:
        nen = ex.shape[0]
        nel = 1
        dxmax=ex.T.max()-ex.T.min()
        dymax=ey.T.max()-ey.T.min()
        dlmax=max(dxmax,dymax);
        edmax=abs(ed).max();
    else:
        nen = ex.shape[1]
        nel = ex.shape[0]
        dxmax=ex.T.max()-ex.T.min()
        dymax=ey.T.max()-ey.T.min()
        dlmax=max(dxmax,dymax);
        edmax=abs(ed).max();
        
    k = rat
    return k*dlmax/edmax
    
    