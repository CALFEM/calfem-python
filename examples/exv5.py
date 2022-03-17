# -*- coding: utf-8 -*-
"""
3D example using Vedo, plate elements

@author: Andreas Åmand
"""
"""
import os
import sys

os.system('clear')
sys.path.append("../")

from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import calfem.core as cfc
import numpy as np
#import vis_vedo as cfvv
import vis_vedo_no_qt as cfvv
from PyQt5 import Qt
"""
import os
import sys

os.system('clear')
sys.path.append("../")
sys.path.append("../../../calfem-python-develop/calfem")

from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
#import calfem.core as cfc
import core as cfc
import numpy as np
#import vis_vedo as cfvv
import vis_vedo_no_qt as cfvv
from PyQt5 import Qt

#import calfem.geometry as cfg
#import calfem.mesh as cfm
#import calfem.vis as cfv
import geometry as cfg
import mesh as cfm
import vis as cfv
import utils as cfu
#import calfem.utils as cfu









def platrs(ex,ey,ep,D,ed):
    # [es,et]=platrs(ex,ey,ep,D,ed)
    #-------------------------------------------------------------
    # PURPOSE
    #  Calculate element normal and shear force for a
    #  rectangular plate element.
    #
    # INPUT: ex = [x1 x2 x3 x4]     element coordinates
    #        ey = [y1 y2 y3 y4]
    #
    #        ep = [t]               thickness
    #
    #        D                      constitutive matrix for
    #                               plane stress
    #
    #        ed = [u1 u2......u12;  element displacement vector
    #             .............   ] one row for each element
    #
    # OUTPUT: es = [ Mxx Myy Mxy Vxz Vyz;   element force matrix
    #                    ......          ]  one row for each element
    #         et = [kxx,kyy,kxy]       curvature in global coordinates
    #-------------------------------------------------------------

    # LAST MODIFIED: Anders Olsson   1999-03-01
    # Copyright (c)  Division of Structural Mechanics and
    #                Department of Solid Mechanics.
    #                Lund Institute of Technology
    #-------------------------------------------------------------
    #
    Lx=ex[2]-ex[0]
    Ly=ey[2]-ey[0]
    t=ep[0]
    #
    D=((t**3)/12)*D
    #
    A1=D[1,1]/2/Ly
    A2=D[0,0]/2/Lx
    A3=D[0,1]/2/Ly
    A4=D[0,1]/2/Lx
    A5=D[2,2]/2/Ly
    A6=D[2,2]/2/Lx
    A7=4*D[2,2]/Lx/Ly

    B1=6*D[1,1]/Ly/Ly/Ly
    B2=6*D[0,0]/Lx/Lx/Lx
    B3=-3*D[1,1]/Ly/Ly
    B4=3*D[0,0]/Lx/Lx
    B5=D[0,1]/Lx/Ly
    #

    mx=A3*(-ed[1]-ed[4]+ed[7]+ed[10])+A2*(ed[2]-ed[5]-ed[8]+ed[11])
    my=A1*(-ed[1]-ed[4]+ed[7]+ed[10])+A4*(ed[2]-ed[5]-ed[8]+ed[11])
    mxy=A6*(ed[1]-ed[4]-ed[7]+ed[10])+A5*(-ed[2]-ed[5]+ed[8]+ed[11])+A7*(ed[0]-ed[3]+ed[6]-ed[9])

    #mx=A3*(-ed[:,1]-ed[:,4]+ed[:,7]+ed[:,10])+A2*(ed[:,2]-ed[:,5]-ed[:,8]+ed[:,11])
    #my=A1*(-ed[:,1]-ed[:,4]+ed[:,7]+ed[:,10])+A4*(ed[:,2]-ed[:,5]-ed[:,8]+ed[:,11])
    #mxy=A6*(ed[:,1]-ed[:,4]-ed[:,7]+ed[:,10])+A5*(-ed[:,2]-ed[:,5]+ed[:,8]+ed[:,11])+A7*(ed[:,0]-ed[:,3]+ed[:,6]-ed[:,9])

    m1=0.5*(mx+my)+np.sqrt(0.25*(mx-my)**2+mxy**2)
    m2=0.5*(mx+my)-np.sqrt(0.25*(mx-my)**2+mxy**2)
    alfa=0.5*180/np.pi*np.arctan2(mxy,(mx-my)/2)

    #m1=0.5*(mx+my)+np.sqrt(0.25*(mx-my)**2+mxy**2)
    #m2=0.5*(mx+my)-np.sqrt(0.25*(mx-my)**2+mxy**2)
    #alfa=0.5*180/np.pi*np.arctan2(mxy,(mx-my)/2)

    vx=B5*(-ed[1]+ed[4]-ed[7]+ed[10])+B4*(ed[2]+ed[5]+ed[8]+ed[11])+B2*(-ed[0]+ed[3]+ed[6]-ed[9])
    vy=B3*(ed[1]+ed[4]+ed[7]+ed[10])+B5*(ed[2]-ed[5]+ed[8]-ed[11])+B1*(-ed[0]-ed[3]+ed[6]+ed[9])
    
    #vx=B5*(-ed[:,1]+ed[:,4]-ed[:,7]+ed[:,10])+B4*(ed[:,2]+ed[:,5]+ed[:,8]+ed[:,11])+B2*(-ed[:,0]+ed[:,3]+ed[:,6]-ed[:,9])
    #vy=B3*(ed[:,1]+ed[:,4]+ed[:,7]+ed[:,10])+B5*(ed[:,2]-ed[:,5]+ed[:,8]-ed[:,11])+B1*(-ed[:,0]-ed[:,3]+ed[:,6]+ed[:,9])

             
    es=np.transpose(np.array([mx, my, mxy, vx, vy]))

    #et=-np.linalg.inv(D)*[mx,my,mxy]
    #et=np.transpose(et)

    #print('shape',es.shape[0],es.shape[1])

    return es
#--------------------------end--------------------------------








d=0.1;
t=0.05

ncoord_x = 30+1;
ncoord_y = 30+1;


ncoord_init = ncoord_x*ncoord_y;

coord = np.zeros([ncoord_init,2]);
row = 0;


for y in range(ncoord_y):
    for x in range(ncoord_x):
        coord[row,:] = [x*d,y*d];
        row = row+1;

ncoord = np.size(coord,0);



dof = np.zeros([ncoord,3]);

it = 1;
dofs = [0,1,2]
for row in range(ncoord):
    for col in dofs:
        dof[row,col] = it;
        it = it + 1;
    #dof[row] = row;


ndof = np.size(dof,0)*np.size(dof,1);



nel_x = (ncoord_x-1);
nel_y = (ncoord_y-1);

#edof = np.zeros([nel_x*nel_y*nel_z,8*3]);
edof = np.zeros((nel_x*nel_y,4*3));
bc = np.zeros((ncoord_y*2*3*2-12,1));

x_step = 1;
y_step = ncoord_x;

it = 0;
bc_it = 0;
node = 0;


for row in range(nel_y):
    for el in range(nel_x):

        #edof[it,0] = it;
        edof[it,0:3] = dof[node,:];
        #edof[it,0] = dof[node];
        
        if el == 0:
            bc[bc_it,0] = int(dof[node,0])
            bc_it = bc_it + 1;
            bc[bc_it,0] = int(dof[node,1])
            bc_it = bc_it + 1;
            bc[bc_it,0] = int(dof[node,2])
            bc_it = bc_it + 1;
        
        node = node+x_step;
        edof[it,3:6] = dof[node,:];
        #edof[it,1] = dof[node];
        print(bc_it)
        print(node)
        if el == nel_x-1:
            bc[bc_it,0] = int(dof[node,0])
            bc_it = bc_it + 1;
            bc[bc_it,0] = int(dof[node,1])
            bc_it = bc_it + 1;
            bc[bc_it,0] = int(dof[node,2])
            bc_it = bc_it + 1;
        
        node = node+y_step;
        edof[it,6:9] = dof[node,:];
        #edof[it,2] = dof[node];
        
        if el == nel_x-1:
            bc[bc_it,0] = int(dof[node,0])
            bc_it = bc_it + 1;
            bc[bc_it,0] = int(dof[node,1])
            bc_it = bc_it + 1;
            bc[bc_it,0] = int(dof[node,2])
            bc_it = bc_it + 1;
        
        node = node-x_step;
        edof[it,9:12] = dof[node,:];
        #edof[it,3] = dof[node];
        
        if el == 0:
            bc[bc_it,0] = int(dof[node,0])
            bc_it = bc_it + 1;
            bc[bc_it,0] = int(dof[node,1])
            bc_it = bc_it + 1;
            bc[bc_it,0] = int(dof[node,2])
            bc_it = bc_it + 1;
        

        if el == nel_x-1:
            #node = node-z_step-y_step-2*x_step+y_step;
            #node = node-z_step-2*x_step;
            node = node-y_step+2
        else:
            node = node+x_step-y_step;
        
        it = it+1;


#coord = np.delete(coord,84,0)
#dof = np.delete(dof,np.size(dof,0)-2,0)
#edof = np.delete(edof,65,0)
#edof = np.delete(edof,66,0)
#edof = np.delete(edof,77,0)
#edof = np.delete(edof,78,0)

#print(edof)

edof = np.int_(edof)


#print(ex)






nnode = np.size(coord, axis = 0)
ndof = np.size(dof, axis = 0)*np.size(dof, axis = 1)
nel = np.size(edof, axis = 0)

ep=[t]

E=25*1000000000
v=0.2
eq=-250*1000

D = cfc.hooke(1,E,v);

ex, ey = cfc.coordxtr(edof,coord,dof)

#print(ex[0])

#K = np.int_(np.zeros((ndof,ndof)))
K = np.zeros((ndof,ndof))
#f = np.int_(np.zeros((ndof,1)))
f = np.zeros((ndof,1))

#print(f)

for eltopo, elx, ely in zip(edof, ex, ey):
    Ke,fe = cfc.platre(elx, ely, ep, D, eq)
    # Transposing fe due to error in platre
    fe = np.transpose(fe)
    cfc.assem(eltopo, K, Ke, f, fe)
"""
for i in range(nel):
    print(ex[i])
    print(ey[i])
    print(ep)
    print(D)
    print(eq)
    Ke, fe = cfc.platre(ex[i], ey[i], ep, D, eq)
    print(fe)
    #print(edof[i,:])
    #print(K)
    K, f = cfc.assem(edof[i,:],K,Ke,f,fe)
"""
bc = bc[:,0]
#bc[bc_it] = 
#for i in range(np.size(bc)):
    #print(int(bc[i,0]))
    #bc = int(bc[i,0])
bc = bc.astype(np.int64)
#print(bc)
extra_bcs = np.array([15*3+1,945*3+1])
bc = np.concatenate((bc, extra_bcs))
#np.insert(bc, 14*3+1, 944*3+1)
#np.append(bc, 14*3+1)
#np.append(bc, 944*3+1)
#bc.append(14*3+1)
#bc.append(944*3+1)
#print(bc)
#bcVal = np.zeros((1,np.size(bc)))
#print(bcVal[0])

a,r = cfc.solveq(K, f, bc)
#print('a',a)

ed = cfc.extract_eldisp(edof,a)
#print(ed)

es = np.zeros((nel,5))

#es = platrs(ex,ey,ep,D,ed)

for i in range(nel):
    es[i,:] = platrs(ex[i],ey[i],ep,D,ed[i])

vM = np.zeros((nel,1))

for i in range(nel):
    sigma_xx = es[i,0]/t
    sigma_yy = es[i,1]/t
    sigma_xy = es[i,2]/t

    sigma_1 = (sigma_xx+sigma_yy)/2 + np.sqrt( ((sigma_xx+sigma_yy)**2)/4 + sigma_xy**2)
    sigma_2 = (sigma_xx+sigma_yy)/2 - np.sqrt( ((sigma_xx+sigma_yy)**2)/4 + sigma_xy**2)

    vM[i] = np.sqrt(sigma_1**2 -sigma_1*sigma_2 + sigma_2**2)

#print(es)

# Node 6 & 162 med 12 element per axel

# Node 15 & 945 med 12 element per axel


cfvv.draw_mesh(edof,coord,dof,6,scale=0.005)
#cfvv.show_and_wait()
cfvv.figure(2)
def_scale = 5
mesh = cfvv.draw_displaced_mesh(edof,coord,dof,6,a,scale=0.002,def_scale=def_scale,scalars=vM/1000000,scalar_title='von Mises [MPa]')
cfvv.add_rulers()
cfvv.add_scalar_bar('von Mises [MPa]')
cfvv.add_text(f'Defomation scale: {def_scale}')

#Start Calfem-vedo visualization
#cfvv.show_and_wait()


# Animation

cfvv.figure(3)

cfvv.animation(edof,coord,dof,6,a,vM/1000000,def_scale=def_scale,dt=250,steps=20)

# For exporting animation
#cfvv.animation(edof,coord,dof,6,a,vM/1000000,def_scale=def_scale,dt=250,steps=20,export=True,file='export/exv5/anim/exv5',scalar_title='von Mises [MPa]')





cfvv.add_scalar_bar('von Mises [MPa]')
cfvv.add_text(f'Defomation scale: {def_scale}')

#Start Calfem-vedo visualization
cfvv.show_and_wait()


#For exporting mesh
#cfvv.export_vtk('export/exv5/exv5', mesh)









