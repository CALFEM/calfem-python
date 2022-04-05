# -*- coding: utf-8 -*-
"""
3D example using Vedo, plate elements

@author: Andreas Åmand
"""

import numpy as np
import calfem.core as cfc
import calfem.vedo as cfv

# Adaptation of platre from version 3.4 for MATLAB
# This function is not in the Python version currently
# Only es is included as it's the only needed output
def platrs(ex,ey,ep,D,ed):
    Lx=ex[2]-ex[0]
    Ly=ey[2]-ey[0]
    t=ep[0]
    
    D=((t**3)/12)*D
    
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

    mx=A3*(-ed[1]-ed[4]+ed[7]+ed[10])+A2*(ed[2]-ed[5]-ed[8]+ed[11])
    my=A1*(-ed[1]-ed[4]+ed[7]+ed[10])+A4*(ed[2]-ed[5]-ed[8]+ed[11])
    mxy=A6*(ed[1]-ed[4]-ed[7]+ed[10])+A5*(-ed[2]-ed[5]+ed[8]+ed[11])+A7*(ed[0]-ed[3]+ed[6]-ed[9])

    m1=0.5*(mx+my)+np.sqrt(0.25*(mx-my)**2+mxy**2)
    m2=0.5*(mx+my)-np.sqrt(0.25*(mx-my)**2+mxy**2)
    alfa=0.5*180/np.pi*np.arctan2(mxy,(mx-my)/2)

    vx=B5*(-ed[1]+ed[4]-ed[7]+ed[10])+B4*(ed[2]+ed[5]+ed[8]+ed[11])+B2*(-ed[0]+ed[3]+ed[6]-ed[9])
    vy=B3*(ed[1]+ed[4]+ed[7]+ed[10])+B5*(ed[2]-ed[5]+ed[8]-ed[11])+B1*(-ed[0]-ed[3]+ed[6]+ed[9])
             
    es=np.transpose(np.array([mx, my, mxy, vx, vy]))

    return es

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

ndof = np.size(dof,0)*np.size(dof,1);

nel_x = (ncoord_x-1);
nel_y = (ncoord_y-1);

edof = np.zeros((nel_x*nel_y,4*3));
bc = np.zeros((ncoord_y*2*3*2-12,1));

x_step = 1;
y_step = ncoord_x;

it = 0;
bc_it = 0;
node = 0;

for row in range(nel_y):
    for el in range(nel_x):

        edof[it,0:3] = dof[node,:];
        
        if el == 0:
            bc[bc_it,0] = int(dof[node,0])
            bc_it = bc_it + 1;
            bc[bc_it,0] = int(dof[node,1])
            bc_it = bc_it + 1;
            bc[bc_it,0] = int(dof[node,2])
            bc_it = bc_it + 1;
        
        node = node+x_step;
        edof[it,3:6] = dof[node,:];
        if el == nel_x-1:
            bc[bc_it,0] = int(dof[node,0])
            bc_it = bc_it + 1;
            bc[bc_it,0] = int(dof[node,1])
            bc_it = bc_it + 1;
            bc[bc_it,0] = int(dof[node,2])
            bc_it = bc_it + 1;
        
        node = node+y_step;
        edof[it,6:9] = dof[node,:];
        
        if el == nel_x-1:
            bc[bc_it,0] = int(dof[node,0])
            bc_it = bc_it + 1;
            bc[bc_it,0] = int(dof[node,1])
            bc_it = bc_it + 1;
            bc[bc_it,0] = int(dof[node,2])
            bc_it = bc_it + 1;
        
        node = node-x_step;
        edof[it,9:12] = dof[node,:];
        
        if el == 0:
            bc[bc_it,0] = int(dof[node,0])
            bc_it = bc_it + 1;
            bc[bc_it,0] = int(dof[node,1])
            bc_it = bc_it + 1;
            bc[bc_it,0] = int(dof[node,2])
            bc_it = bc_it + 1;
        
        if el == nel_x-1:
            node = node-y_step+2
        else:
            node = node+x_step-y_step;
        
        it = it+1;

edof = np.int_(edof)

nnode = np.size(coord, axis = 0)
ndof = np.size(dof, axis = 0)*np.size(dof, axis = 1)
nel = np.size(edof, axis = 0)

ep=[t]

E=25*1000000000
v=0.2
eq=-250*1000

D = cfc.hooke(1,E,v);


ex, ey = cfc.coordxtr(edof,coord,dof)

K = np.zeros((ndof,ndof))
f = np.zeros((ndof,1))

for eltopo, elx, ely in zip(edof, ex, ey):
    Ke,fe = cfc.platre(elx, ely, ep, D, eq)
    # Transposing fe due to error in platre
    fe = np.transpose(fe)
    cfc.assem(eltopo, K, Ke, f, fe)

bc = bc[:,0]
bc = bc.astype(np.int64)
extra_bcs = np.array([15*3+1,945*3+1])
bc = np.concatenate((bc, extra_bcs))

a,r = cfc.solveq(K, f, bc)

ed = cfc.extract_eldisp(edof,a)
es = np.zeros((nel,5))

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
cfv.draw_mesh(edof,coord,dof,6,scale=0.005)

cfv.figure(2)
def_scale = 5
mesh = cfv.draw_displaced_mesh(edof,coord,dof,6,a,scale=0.002,def_scale=def_scale,scalars=vM/1000000,scalar_title='von Mises [MPa]')
cfv.add_rulers()
cfv.add_scalar_bar('von Mises [MPa]')
cfv.add_text(f'Defomation scale: {def_scale}')

#Start Calfem-vedo visualization
cfv.show_and_wait()

cfv.figure(3)

cfv.animation(edof,coord,dof,6,a,vM/1000000,def_scale=def_scale,dt=250,steps=20,export=True,file='export/exv5/anim/exv5',scalar_title='von Mises [MPa]')
cfv.add_scalar_bar('von Mises [MPa]')
cfv.add_text(f'Defomation scale: {def_scale}')

#Start Calfem-vedo visualization
cfv.show_and_wait()

# Export the mesh
cfv.export_vtk('export/exv5/exv5', mesh)

# For not exporting animation
#cfv.animation(edof,coord,dof,6,a,vM/1000000,def_scale=def_scale,dt=250,steps=20)

