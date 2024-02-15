import numpy as np
import calfem.core as cfc

def soli8e(ex, ey, ez, ep, D, eqp=None):
    """
    Ke=soli8e(ex,ey,ez,ep,D)
    [Ke,fe]=soli8e(ex,ey,ez,ep,D,eq)
    -------------------------------------------------------------
    PURPOSE
    Calculate the stiffness matrix for a 8 node (brick)
    isoparametric element.

    INPUT:   ex = [x1 x2 x3 ... x8]
            ey = [y1 y2 y3 ... y8]  element coordinates
            ez = [z1 z2 z3 ... z8]

            ep = [ir]               ir integration rule

            D                       constitutive matrix

            eq = [bx; by; bz]       bx: body force in x direction
                                    by: body force in y direction
                                    bz: body force in z direction

    OUTPUT: Ke : element stiffness matrix
            fe : equivalent nodal forces 
    -------------------------------------------------------------

    LAST MODIFIED: M Ristinmaa   1995-10-25
    Copyright (c)  Division of Structural Mechanics and
                Division of Solid Mechanics.
                Lund University
    -------------------------------------------------------------    
    """
    #ir=ep(1),  ngp=ir*ir*ir,
    ir = ep[0]
    ngp = ir*ir*ir

    if eqp==None:
        eq = np.zeros((3,1))
    else:
        eq = eqp
    
    if ir==1:
        g1 = 0.0
        w1 = 2.0
        gp = np.array([ g1, g1, g1 ]).reshape(1,3)
        w = np.array([ w1, w1, w1 ]).reshape(1,3)
    elif ir==2:
        g1=0.577350269189626
        w1=1
        gp = np.zeros((8,3))
        w = np.zeros((8,3))
        gp[:,0] = np.array([-1, 1, 1,-1,-1, 1, 1,-1])*g1 
        w[:,0]  = np.array([ 1, 1, 1, 1, 1, 1, 1, 1])*w1
        gp[:,1] = np.array([-1,-1, 1, 1,-1,-1, 1, 1])*g1
        w[:,1]  = np.array([ 1, 1, 1, 1, 1, 1, 1, 1])*w1
        gp[:,2] = np.array([-1,-1,-1,-1, 1, 1, 1, 1])*g1
        w[:,2]  = np.array([ 1, 1, 1, 1, 1, 1, 1, 1])*w1
    else:
        g1=0.774596669241483,
        g2=0.0
        w1=0.555555555555555
        w2=0.888888888888888

        gp = np.zeros((27,3))
        w = np.zeros((27,3))
        
        I1=np.array([-1, 0, 1,-1, 0, 1,-1, 0, 1]).reshape(1,9)
        I2=np.array([ 0,-1, 0, 0, 1, 0, 0, 1, 0]).reshape(1,9)
        #gp[:,0]=[I1 I1 I1]'*g1
        #gp[:,0]=[I2 I2 I2]'*g2+gp[:,0]
        gp[:,0]=np.concatenate((I1, I1, I1), axis=1)*g1
        gp[:,0]=np.concatenate((I2, I2, I2), axis=1)*g2 + gp[:,0]
        I1=np.abs(I1)
        I2=np.abs(I2)
        #w[:,0]=[I1 I1 I1]'*w1
        #w[:,0]=[I2 I2 I2]'*w2+w[:,0]
        w[:,0]=np.concatenate((I1, I1, I1), axis=1)*w1
        w[:,0]=np.concatenate((I2, I2, I2), axis=1)*w2 + w[:,0]

        #I1=[-1,-1,-1, 0, 0, 0, 1, 1, 1]'
        #I2=[ 0, 0, 0, 1, 1, 1, 0, 0, 0]'
        I1=np.array([-1,-1,-1, 0, 0, 0, 1, 1, 1]).reshape(1,9)
        I2=np.array([ 0, 0, 0, 1, 1, 1, 0, 0, 0]).reshape(1,9)
        #gp[:,1]=[I1 I1 I1]'*g1
        #gp[:,1]=[I2 I2 I2]'*g2+gp[:,1]
        gp[:,1]=np.concatenate((I1, I1, I1), axis=1)*g1
        gp[:,1]=np.concatenate((I2, I2, I2), axis=1)*g2 + gp[:,1]

        I1=np.abs(I1)
        I2=np.abs(I2)

        #w[:,1]=[I1 I1 I1]'*w1
        #w[:,1]=[I2 I2 I2]'*w2+w[:,1]
        w[:,1]=np.concatenate((I1, I1, I1), axis=1)*w1
        w[:,1]=np.concatenate((I2, I2, I2), axis=1)*w2 + w[:,1]
        
        #I1=[-1,-1,-1,-1,-1,-1,-1,-1,-1]'
        #I2=[ 0, 0, 0, 0, 0, 0, 0, 0, 0]'
        I1=np.array([-1,-1,-1,-1,-1,-1,-1,-1,-1]).reshape(1,9)
        I2=np.array([ 0, 0, 0, 0, 0, 0, 0, 0, 0]).reshape(1,9)
        I3=np.abs(I1)

        #gp[:,2]=[I1 I2 I3]'*g1
        #gp[:,2]=[I2 I3 I2]'*g2+gp[:,2)
        gp[:,2]=np.concatenate((I1, I2, I3), axis=1)*g1
        gp[:,2]=np.concatenate((I2, I3, I2), axis=1)*g2 + gp[:,2]

        #w[:,2]=[I3 I2 I3]'*w1
        #w[:,2]=[I2 I3 I2]'*w2+w[:,2)
        w[:,2]=np.concatenate((I3, I2, I3), axis=1)*w1
        w[:,2]=np.concatenate((I2, I3, I2), axis=1)*w2 + w[:,2]

    #wp=w(:,1).*w(:,2).*w(:,3);
    wp = w[:,0]*w[:,1]*w[:,2]

    #xsi=gp(:,1);  eta=gp(:,2); zet=gp(:,3);  r2=ngp*3;
    xsi = gp[:,0]
    eta = gp[:,1]
    zet = gp[:,2]
    r2 = ngp*3

    N = np.zeros((ngp,8))
    dNr = np.zeros((r2,8))

    N[:,0]=(1-xsi)*(1-eta)*(1-zet)/8
    N[:,1]=(1+xsi)*(1-eta)*(1-zet)/8
    N[:,2]=(1+xsi)*(1+eta)*(1-zet)/8
    N[:,3]=(1-xsi)*(1+eta)*(1-zet)/8
    N[:,4]=(1-xsi)*(1-eta)*(1+zet)/8
    N[:,5]=(1+xsi)*(1-eta)*(1+zet)/8
    N[:,6]=(1+xsi)*(1+eta)*(1+zet)/8
    N[:,7]=(1-xsi)*(1+eta)*(1+zet)/8

    dNr[0:r2+1:3,0]=-(1-eta)*(1-zet);  dNr[0:r2+1:3,1]= (1-eta)*(1-zet)
    dNr[0:r2+1:3,2]= (1+eta)*(1-zet);  dNr[0:r2+1:3,3]=-(1+eta)*(1-zet)
    dNr[0:r2+1:3,4]=-(1-eta)*(1+zet);  dNr[0:r2+1:3,5]= (1-eta)*(1+zet)
    dNr[0:r2+1:3,6]= (1+eta)*(1+zet);  dNr[0:r2+1:3,7]=-(1+eta)*(1+zet)
    dNr[1:r2+2:3,0]=-(1-xsi)*(1-zet);  dNr[1:r2+2:3,1]=-(1+xsi)*(1-zet)
    dNr[1:r2+2:3,2]= (1+xsi)*(1-zet);  dNr[1:r2+2:3,3]= (1-xsi)*(1-zet)
    dNr[1:r2+2:3,4]=-(1-xsi)*(1+zet);  dNr[1:r2+2:3,5]=-(1+xsi)*(1+zet)
    dNr[1:r2+2:3,6]= (1+xsi)*(1+zet);  dNr[1:r2+2:3,7]= (1-xsi)*(1+zet)
    dNr[2:r2+3:3,0]=-(1-xsi)*(1-eta);  dNr[2:r2+3:3,1]=-(1+xsi)*(1-eta)
    dNr[2:r2+3:3,2]=-(1+xsi)*(1+eta);  dNr[2:r2+3:3,3]=-(1-xsi)*(1+eta)
    dNr[2:r2+3:3,4]= (1-xsi)*(1-eta);  dNr[2:r2+3:3,5]= (1+xsi)*(1-eta)
    dNr[2:r2+3:3,6]= (1+xsi)*(1+eta);  dNr[2:r2+3:3,7]= (1-xsi)*(1+eta)

    dNr = dNr/8.0

    Ke = np.zeros((24,24))
    fe = np.zeros((24,1))

    # [81 x 8] x [8 x 3]

    ex = np.asarray(ex).reshape((8,1))
    ey = np.asarray(ey).reshape((8,1))
    ez = np.asarray(ez).reshape((8,1))

    JT = dNr@np.concatenate((ex, ey, ez), axis=1)

    eps = np.finfo(float).eps

    for i in range(ngp):
        indx = [ i*3, i*3+1, i*3+2]
        detJ = np.linalg.det(JT[indx,:])
        if detJ<10*eps:
            print('Jacobideterminant equal or less than zero!')
        JTinv=np.linalg.inv(JT[indx,:])
        dNx=JTinv@dNr[indx,:]

        B=np.zeros((6,24))     
        N2=np.zeros((3,24))

        B[0,0:24:3] = dNx[0,:]
        B[1,1:25:3] = dNx[1,:]
        B[2,2:26:3] = dNx[2,:]
        B[3,0:24:3] = dNx[1,:]
        B[3,1:25:3] = dNx[0,:]
        B[4,0:24:3] = dNx[2,:]
        B[4,2:26:3] = dNx[0,:]
        B[5,1:25:3] = dNx[2,:]
        B[5,2:26:3] = dNx[1,:]

        N2[0,0:24:3] = N[i,:]
        N2[1,1:25:3] = N[i,:]
        N2[2,2:26:3] = N[i,:]

        Ke = Ke + (np.transpose(B)@D@B)*detJ*wp[i]
        fe = fe + (np.transpose(N2)@eq)*detJ*wp[i]

    if eqp!=None:
        return Ke, fe
    else:
        return Ke

if __name__ == "__main__":

    np.set_printoptions(linewidth=np.inf, precision=4)

    D = cfc.hooke(4, 2.1e9, 0.35)

    ex = [0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0]
    ey = [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0]
    ez = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]

    Ke = soli8e(ex, ey, ez, [1], D)
    print(Ke)
