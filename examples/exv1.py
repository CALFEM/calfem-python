# -*- coding: utf-8 -*-
"""
Example using Vedo, spring elements

@author: Andreas Åmand
"""

import os
import sys

os.system('clear')
sys.path.append("../")

import numpy as np # Ska vara line 8 -> 6 tillkommer
import calfem.core as cfc
import vis_vedo_no_qt as cfvv

coord = np.array([
    [0],
    [0.5],
    [1],
    [1.5]
])

dof = np.array([
    [1],
    [2],
    [3],
    [4]
])

edof = np.array([
    [1, 2],
    [2, 3],
    [3, 4]
])

k = 1000
ep = [3*k, k, 8*k]

ndof = dof.shape[0]*dof.shape[1]
#ndof = np.size(dof, axis = 0)*np.size(dof, axis = 1)
nel = edof.shape[0]
#nel = np.size(edof, axis = 0)

K = np.zeros([ndof,ndof])

for i in range(nel):
    Ke = cfc.spring1e(ep[i])
    K = cfc.assem(edof[i],K,Ke)

f = np.zeros([ndof,1])
f[3,0] = 500 #Newton

bcPrescr = np.array([1])
a,r = cfc.solveq(K, f, bcPrescr)

#ed = cfc.extractEldisp(edof,a)

#es = np.zeros((3,1))
#for i in range(nel):
#    es[i] = cfc.spring1s(ep[i],ed[i])

#print(es)

print('coord',coord)

cfvv.figure(1,flat=True)
cfvv.draw_mesh(edof,coord,dof,1)
mesh = cfvv.draw_displaced_mesh(edof,coord,dof,1,a,offset=[0,0.2,0],render_nodes=True)
#cfvv.add_scalar_bar('Normal force [kN]')
cfvv.add_text_3D('k=3 kN/m',[0.15,-0.1,0],size=0.03)
cfvv.add_text_3D('k=1 kN/m',[0.65,-0.1,0],size=0.03)
cfvv.add_text_3D('k=8 kN/m',[1.15,-0.1,0],size=0.03)
cfvv.add_text_3D('F_x =500 N',[1.55,-0.02,0],size=0.03)
#cfvv.add_rulers()
#cfvv.draw_displaced_geometry(edof,coord,dof,5,a,scale=0.5,alpha=0.5)

#Start Calfem-vedo visualization
#cfvv.show_and_wait()


# Animation

cfvv.figure(2)
steps = 20
cfvv.add_text(f'Looping bewteen undef. & def. state w/ {steps} steps',pos='top-middle')
cfvv.animation(edof,coord,dof,1,a,loop=True,steps=20,dt=0)

# For exporting animation
#cfvv.animation(edof,coord,dof,1,a,loop=True,steps=20,dt=0,export=True,file='export/exv1/anim/exv1')

#Start Calfem-vedo visualization
cfvv.show_and_wait()



# For exporting mesh
#cfvv.export_vtk('export/exv1/exv1', mesh)






