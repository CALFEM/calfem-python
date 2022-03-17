# -*- coding: utf-8 -*-
"""
3D example using Vedo, flow elements

@author: Andreas Åmand
"""

import os
import sys

os.system('clear')
sys.path.append("../")
sys.path.append("../../../calfem-python-develop/calfem")

#from PyQt5.QtWidgets import *
#from PyQt5.QtCore import *
#import calfem.core as cfc

import numpy as np
#import vis_vedo as cfvv
import vis_vedo_no_qt as cfvv
#from PyQt5 import Qt

#import calfem.geometry as cfg
#import calfem.mesh as cfm
#import calfem.vis as cfv
import core as cfc
import geometry as cfg
import mesh as cfm
import vis as cfv
import utils as cfu
#import calfem.utils as cfu



# x=0.4, y=1, z=0.4

g = cfg.Geometry()
"""
g.point([0.0, 0.0, 0.0]) # point 0
g.point([0.0, 0.0, 0.4]) # point 1
g.point([0.0, 1.0, 0.0]) # point 2
g.point([0.0, 1.0, 0.4]) # point 3
g.point([0.4, 0.0, 0.0]) # point 4
g.point([0.4, 0.0, 0.4]) # point 5
g.point([0.4, 1.0, 0.0]) # point 6
g.point([0.4, 1.0, 0.4]) # point 7


g.spline([0, 1]) # line 0
g.spline([1, 5]) # line 1
g.spline([5, 4]) # line 2
g.spline([4, 0]) # line 2

g.spline([2, 3]) # line 0
g.spline([3, 7]) # line 1
g.spline([7, 6]) # line 2
g.spline([6, 2]) # line 2

g.spline([0, 2]) # line 0
g.spline([1, 3]) # line 1
g.spline([5, 7]) # line 2
g.spline([4, 6]) # line 2

g.surface([0, 1, 5, 4])
g.surface([2, 3, 7, 6])

g.surface([0, 1, 3, 2])
g.surface([1, 5, 7, 3])
g.surface([5, 4, 6, 7])
g.surface([4, 6, 2, 0])

g.volume([0,1,2,3,4,5])
"""

l = 1
h = 0.4
w = 0.4

n_el_x = 4
n_el_y = 4
n_el_z = 10

g.point([0, 0, 0], ID=0)
#g.point([w/2.0, 0.0, 0.0], 1)
g.point([w, 0, 0], 1)
g.point([w, l, 0], 2, marker=50)
g.point([0, l, 0], 3)
g.point([0, 0, h], 4)
g.point([w, 0, h], 5)
g.point([w, l, h], 6)
g.point([0, l, h], 7)

# Add splines

g.spline([0, 1], 0, el_on_curve = n_el_x)
g.spline([1, 2], 1, el_on_curve = n_el_z)
g.spline([2, 3], 2, el_on_curve = n_el_x)
g.spline([3, 0], 3, el_on_curve = n_el_z)
g.spline([0, 4], 4, el_on_curve = n_el_y)
g.spline([1, 5], 5, el_on_curve = n_el_y)
g.spline([2, 6], 6, el_on_curve = n_el_y)
g.spline([3, 7], 7, el_on_curve = n_el_y)
g.spline([4, 5], 8, el_on_curve = n_el_x)
g.spline([5, 6], 9, el_on_curve = n_el_z)
g.spline([6, 7], 10, el_on_curve = n_el_x)
g.spline([7, 4], 11, el_on_curve = n_el_z)
'''
g.spline([0, 1], 0)
g.spline([1, 2], 1)
g.spline([2, 3], 2)
g.spline([3, 0], 3)
g.spline([0, 4], 4)
g.spline([1, 5], 5)
g.spline([2, 6], 6)
g.spline([3, 7], 7)
g.spline([4, 5], 8)
g.spline([5, 6], 9)
g.spline([6, 7], 10)
g.spline([7, 4], 11)
'''
# Add surfaces

marker_bottom = 40
marker_top = 41
marker_fixed_left = 42
marker_back = 43
marker_fixed_right = 44
marker_front = 45

#g.surface([0, 1, 2, 3], 0, marker=marker_bottom)
#g.Surface([8, 9, 10, 11], 1, marker=marker_top)
#g.Surface([0, 4, 8, 5], 2, marker=marker_fixed_left)
#g.Surface([1, 5, 9, 6], 3, marker=marker_back)
#g.Surface([2, 6, 10, 7], 4, marker=marker_fixed_right)
#g.Surface([3, 4, 11, 7], 5, marker=marker_front)

g.structuredSurface([0, 1, 2, 3], 0, marker=marker_bottom)
g.structuredSurface([8, 9, 10, 11], 1, marker=marker_top)
g.structuredSurface([0, 4, 8, 5], 2, marker=marker_fixed_left)
g.structuredSurface([1, 5, 9, 6], 3, marker=marker_back)
g.structuredSurface([2, 6, 10, 7], 4, marker=marker_fixed_right)
g.structuredSurface([3, 4, 11, 7], 5, marker=marker_front)

g.structuredVolume([0,1,2,3,4,5], 0, marker=90)






el_type = 5 
dofs_per_node = 1
elSizeFactor = 0.01

# Create mesh

#coord, edof, dof, bdof, elementmarkers = cfm.GmshMeshGenerator.create(g, el_type=5, dofs_per_node=1, el_size_factor=0.1)
#coord, edof, dof, bdof, elementmarkers = cfm.mesh(g, el_type, 0.1, dofs_per_node)
#mesh.create()
"""
mesh = cfm.GmshMesh(g)
mesh.elType = 5             # Type of mesh
mesh.dofsPerNode = 1        # Factor that changes element sizes
mesh.elSizeFactor = 0.01    # Factor that changes element sizes
coord, edof, dof, bdof, elementmarkers = mesh.create()
"""
coord, edof, dof, bdof, elementmarkers = cfm.mesh(g, el_type, elSizeFactor, dofs_per_node)
#print(edof[0])
print('coord',coord)
#print(dof)
#print(bdof)
ex, ey, ez = cfc.coordxtr(edof, coord, dof)







nnode = np.size(coord, axis = 0)
ndof = np.size(dof, axis = 0)*np.size(dof, axis = 1)
nel = np.size(edof, axis = 0)

Lambda = 1.4 # Thermal conductivity for concrete
D = np.identity(3)*Lambda

#print(D)
#print(ex)
#bcPrescr = np.array([1,2,3,4,5,51,52,53,54,55,56,57,58,59,60,111,112,113,114,115,166,167,168,169,170,221,222,223,224,225])
#bcPrescr = np.array([1,2,3,4,5,56,57,58,59,60,111,112,113,114,115,147,166,167,168,169,170,221,222,223,224,225])
#bcPrescr = np.array([[1],[2],[3],[4],[5],[56],[57],[58],[59],[60],[111],[112],[113],[114],[115],[166],[167],[168],[169],[170],[221],[222],[223],[224],[225]])
#bc = np.array([147, 30])
#bc = np.zeros((1,26))
#bc[:,15] = 30
#print(bcPrescr)
#print(bc)

#K = np.int_(np.zeros((ndof,ndof)))
K = np.zeros((ndof,ndof))

f = np.zeros([ndof,1])
eq = np.zeros([nel,1])

print('eq',eq)

#eq_el = 140

eq_els = np.array([[70],[89]])

#eq_el1 = 66
#eq_el2 = 93

eq[eq_els[0]] = 30000
eq[eq_els[1]] = -30000

#print(edof[0,:])
#Ke = cfc.flw3i8e(ex, ey, ez, ep, D)
#K = cfc.assem(edof,K,Ke)
for eltopo, elx, ely, elz, eqs in zip(edof, ex, ey, ez, eq):
    #Ke = cfc.flw3i8e(elx, ely, elz, [2], D)
    Ke,fe = cfc.flw3i8e(elx, ely, elz, [2], D, eqs)
    #print(Ke)
    #cfc.assem(eltopo, K, Ke)
    cfc.assem(eltopo, K, Ke, f, fe)
"""
for i in range(nel):
    #Ke, fe = cfc.flw3i8e(ex[i], ey[i], ez[i], ep, D, eq[i])
    Ke = cfc.flw3i8e(ex[i], ey[i], ez[i], ep, D)
    #print(Ke)
    #print(edof[i,:])
    #print(K)
    #print(i)
    #K,f = cfc.assem(edof[i,:],K,Ke,f,fe)
    K = cfc.assem(edof[i,:],K,Ke)
    #Ke = cfc.flw3i8e(ex[i,:], ey[i,:], ez[i,:], ep, D)
    #print(Ke)
    #K = cfc.assem(edof[i,:],K,Ke)
"""
#f = np.zeros([ndof,1])
#f[7,0] = -3000
bc = np.array([],'i')
bcVal = np.array([],'i')
"""
bc = np.array([],'i')
bcVal = np.array([],'i')
bcs = np.array([
    [2],[20], [38], [158]
    ])
bc_50 = np.array([
    [50],[50], [50], [50]
    ])
np.append(bc, bcs, axis=0)
np.append(bcVal, bc_50, axis=0)
"""
#bc = np.array([2,20,38,158])
#bc = np.array([1,2,4,5,6,  9,10,11,   23,42,165])
#bcVal = np.array([0,0,50,0,0,  0,0,0,    50,50,50])
#bc = np.array([1,2,5,6,  9,10,11])
#bcVal = np.array([0,0,0,0,  0,0,0])

#sys.exit()
# Nearly all surfaces are 20 deg
#bc, bcVal = cfu.apply_bc_3d(bdof, bc, bcVal, marker_bottom, 20.0)
#bc, bcVal = cfu.apply_bc_3d(bdof, bc, bcVal, marker_top, 20.0)
bc, bcVal = cfu.apply_bc_3d(bdof, bc, bcVal, marker_fixed_left, 20.0) # bottom is 0
#bc, bcVal = cfu.apply_bc_3d(bdof, bc, bcVal, marker_back, 20.0)
bc, bcVal = cfu.apply_bc_3d(bdof, bc, bcVal, marker_fixed_right, 20.0) # top is 20
#bc, bcVal = cfu.apply_bc_3d(bdof, bc, bcVal, marker_front, 20.0)
#bc, bcVal = cfu.apply_bc_3d(bdof, bc, bcVal, 50, 50.0)

#f = np.zeros([ndof,1])
#cfu.apply_force_total_3d(bdof, f, marker_top, value = 30, dimension=1)

# Visually accuires 4 nodes for a top corner, these are to be 30 deg
# They are found to be 2, 20, 38, 158
#for i in range(len(bc)):
#    if bc[i] == 2 or bc[i] == 20 or bc[i] == 38 or bc[i] == 158:
#        bcVal[i] = 50
#    print('DOF: ',bc[i],' value: ',bcVal[i])


print('DOF: ',bc,' value: ',bcVal)

#bcPrescr = np.array([1,2,3,4,5,6,13,14,15,16,17,18])
#print(bc)
#print(bcVal)
#print(K[0,:])
#print(f)
#print(type(bc))
#print(bc)
#print(bcVal)
T,r = cfc.solveq(K, f, bc, bcVal)

#ed = cfc.extractEldisp(edof,T)
#print(edof[0])
#print(T)
#print(np.size(T))
ed = cfc.extract_eldisp(edof,T)
#print(ed)

#print(ed)
es = np.zeros((8,3,nel))
edi = np.zeros((8,3,nel))
eci = np.zeros((8,3,nel))

#print('ex: ',ex[0])
#print('ey: ',ey[0])
#print('ez: ',ez[0])
#print('ed: ',ed[0])

for i in range(nel):
    es[0:8,:,i], edi[0:8,:,i], eci[0:8,:,i] = cfc.flw3i8s(ex[i],ey[i],ez[i],[2],D,ed[i])
#es, et, eci = cfc.flw3i8s(ex, ey, ez, [2], D, ed)
"""
# --- Gauss points & shape functions from soli8e/soli8s ---
# This is used to get stresses at nodes & modal analysis later

g1=0.577350269189626;
gp = np.mat([
    [-1,-1,-1],
    [ 1,-1,-1],
    [ 1, 1,-1],
    [-1, 1,-1],
    [-1,-1, 1],
    [ 1,-1, 1],
    [ 1, 1, 1],
    [-1, 1, 1]
])*g1

xsi = gp[:,0]
eta = gp[:,1]
zet = gp[:,2]

N = np.multiply(np.multiply((1-xsi),(1-eta)),(1-zet))/8.
N = np.append(N,np.multiply(np.multiply((1+xsi),(1-eta)),(1-zet))/8.,axis=1)
N = np.append(N,np.multiply(np.multiply((1+xsi),(1+eta)),(1-zet))/8.,axis=1)
N = np.append(N,np.multiply(np.multiply((1-xsi),(1+eta)),(1-zet))/8.,axis=1)
N = np.append(N,np.multiply(np.multiply((1-xsi),(1-eta)),(1+zet))/8.,axis=1)
N = np.append(N,np.multiply(np.multiply((1+xsi),(1-eta)),(1+zet))/8.,axis=1)
N = np.append(N,np.multiply(np.multiply((1+xsi),(1+eta)),(1+zet))/8.,axis=1)
N = np.append(N,np.multiply(np.multiply((1-xsi),(1+eta)),(1+zet))/8.,axis=1)

#print(N)

calc = ((np.transpose(N) * N) / np.transpose(N)) # saving for quicker calculation
"""


#ns = np.zeros((nel,8));
#print(ns[:,0,0])
#print(calc)
#print(es[:,0,0])
vectors = np.zeros((nel,3));
flux = np.zeros((nel,3));
flux_tot = np.zeros((nel,1));
for i in range(nel):
    #print('Flux at element in x-dir.',[es[:,0,i]])
    #print('Flux at element in y-dir.',[es[:,1,i]])
    #print('Flux at element in z-dir.',[es[:,2,i]])
    #flux[i] = np.sqrt(np.average([es[:,0,i]])**2 + np.average([es[:,1,i]])**2 + np.average([es[:,2,i]])**2)
    #print('Flux at element',flux[i])
    
    flux[i,:] = [np.average([es[:,0,i]]), np.average([es[:,1,i]]), np.average([es[:,2,i]])]
    #flow_x = calc*np.transpose([es[:,0,i]])
    #print('Flux at element in x-dir.',flux_x)
    #flux_y = np.average([es[:,1,i]])
    #flow_y = calc*np.transpose([es[:,1,i]])
    #print('Flux at element in y-dir.',flux_y)
    #flux_z = np.average([es[:,2,i]])
    #flow_z = calc*np.transpose([es[:,2,i]])
    #print('Flux at element in z-dir.',flux_z)

    flux_tot[i] = np.sqrt(flux[i,0]**2 + flux[i,1]**2 + flux[i,2]**2)
    #print('Flux at element',flux[i])

    




    #vectors[i,0] = np.average(flow_x)
    #vectors[i,1] = np.average(flow_y)
    #vectors[i,2] = np.average(flow_z)

    #ns[i,0] = np.sqrt(flow_x[0]**2 + flow_y[0]**2 + flow_z[0]**2)
    #ns[i,1] = np.sqrt(flow_x[1]**2 + flow_y[1]**2 + flow_z[1]**2)
    #ns[i,2] = np.sqrt(flow_x[2]**2 + flow_y[2]**2 + flow_z[2]**2)
    #ns[i,3] = np.sqrt(flow_x[3]**2 + flow_y[3]**2 + flow_z[3]**2)
    #ns[i,4] = np.sqrt(flow_x[4]**2 + flow_y[4]**2 + flow_z[4]**2)
    #ns[i,5] = np.sqrt(flow_x[5]**2 + flow_y[5]**2 + flow_z[5]**2)
    #ns[i,6] = np.sqrt(flow_x[6]**2 + flow_y[6]**2 + flow_z[6]**2)
    #ns[i,7] = np.sqrt(flow_x[7]**2 + flow_y[7]**2 + flow_z[7]**2)

#ns = np.transpose(ns)
#print(ns[0])
#print('vectors: ',vectors)










'''

edof,coord,dof,T,ed,bc,eq,es,et,eci = cfvv.import_mat('exv3',['edof','coord','dof','a','ed','bc','eq','es','et','eci'])

nel = np.size(edof,0)
nnode = np.size(coord,0)

#print(bc)

bcVal = bc[:,1]
bc = bc[:,0]
eq_el = np.array([[10],[149]])


flux_tot = np.zeros((nel,1));
for i in range(nel):
    print('Flux at element in x-dir.',es[i,0])
    #flux_y = np.average([es[:,1,i]])
    #flow_y = calc*np.transpose([es[:,1,i]])
    print('Flux at element in y-dir.',es[i,1])
    #flux_z = np.average([es[:,2,i]])
    #flow_z = calc*np.transpose([es[:,2,i]])
    print('Flux at element in z-dir.',es[i,2])
    #flux[i,:] = [np.average([es[:,0,i]]), np.average([es[:,1,i]]), np.average([es[:,2,i]])]
    flux_tot[i] = np.sqrt(es[i,0]**2 + es[i,1]**2 + es[i,2]**2)

ex,ey,ez = cfc.coordxtr(edof,coord,dof)

'''


#print(bc,bcVal)

# For simple point coordinates of geometry
#points = g.getPointCoords()
points = g.points

# For lines with point connectivity
lines = g.curves

# For surfaces with line connectivity
surfaces = g.surfaces

cfvv.draw_geometry(points,lines,surfaces)
#cfvv.draw_mesh(edof,coord,dof,3,scale=0.002)
#cfvv.show_and_wait()
#sys.exit()

print('eq',eq)

cfvv.figure(2)
cfvv.draw_mesh(edof,coord,dof,3,alpha=1,scale=0.005,bcPrescr=bc, bc=bcVal, eq_els=eq_els, eq=eq[eq_els])
#cfvv.add_axes()
#cfvv.draw_mesh(edof,coord,dof,3,alpha=1,scale=0.005,bcPrescr=bc, bc=bcVal, eq_els=[eq_el], eq=eq[eq_el,:])
#cfvv.draw_mesh(edof,coord,dof,3,scale=0.002)
#cfvv.show_and_wait()



disp = np.zeros((nnode,1))

cfvv.figure(3)
#print('flux',flux,flux.shape[0],flux.shape[1])
cfvv.draw_displaced_mesh(edof,coord,dof,3,scalars=flux_tot,colormap='coolwarm',wireframe=True)
#cfvv.add_axes(xrange=1,yrange=1,zrange=1)
cfvv.add_axes(xrange=[-0.1,0.5], yrange=[-0.1,1.1], zrange=[-0.1,0.5])
#cfvv.add_axes()

#cfvv.eliso(mesh)
#cfvv.elcont(mesh)
cfvv.add_scalar_bar('Heat flux [W/m^2]')
cfvv.elflux(ex,ey,ez,flux,colormap='coolwarm')
#cfvv.draw_displaced_mesh(edof,coord,dof,3,disp,ns,colormap='coolwarm')
#cfvv.add_scalar_bar('Heat flow rate [W]')
#cfvv.draw_mesh(edof,coord,dof,3,scale=0.002)
#cfvv.show_and_wait()

cfvv.figure(4)
#print('T',T,T.shape[0],T.shape[1])
mesh = cfvv.draw_displaced_mesh(edof,coord,dof,3,scalars=ed,colormap='coolwarm',colors=5,scalar_title='Temp. [C]')
cfvv.add_mesh_axes(mesh)
cfvv.add_scalar_bar('Temp. [C]')
#cfvv.draw_mesh(edof,coord,dof,3,scale=0.002)
#print(T[153],T[154])
#cfvv.show_and_wait()



# Animation

cfvv.figure(5)

cfvv.animation(edof,coord,dof,3,scalars=ed,dt=250,steps=20,colormap='coolwarm',colors=5)

# For exporting animation
#cfvv.animation(edof,coord,dof,3,scalars=ed,dt=250,steps=20,colormap='coolwarm',colors=5,export=True,file='export/exv3/anim/exv3',scalar_title='Temp. [C]')

cfvv.add_scalar_bar('Temp. [C]')

#Start Calfem-vedo visualization
cfvv.show_and_wait()



# For exporting mesh
#cfvv.export_vtk('export/exv3/exv3', mesh)



#cfvv.figure(5)
#cfvv.eliso(mesh)
#cfvv.show_and_wait()

#cfvv.figure(5)
#cfvv.draw_displaced_mesh(edof,coord,dof,3,disp,T,alpha=0.5,scale=0.01)
#cfvv.add_vectors(edof,coord,dof,flux,3)
#cfvv.add_scalar_bar('Temp. [C]')
#cfvv.draw_mesh(edof,coord,dof,3,scale=0.002)
#cfvv.show_and_wait()


    #ns[:,:,i] = calc*np.transpose([es[:,0,i]])
    #ns[:,:,i] = calc*np.transpose([es[:,1,i]])#calc*es[:,1,i]
    #ns[:,:,i] = calc*np.transpose([es[:,2,i]])#calc*es[:,2,i]

#print(ns)

#ns_test

"""
ngp=ep[1]*ep[1]*ep[1];

g1=0.577350269189626; w1=1;
gp[:,1]=[-1, 1, 1,-1,-1, 1, 1,-1]*g1; w[:,1]=[ 1, 1, 1, 1, 1, 1, 1, 1]*w1;
gp[:,2]=[-1,-1, 1, 1,-1,-1, 1, 1]*g1; w[:,2]=[ 1, 1, 1, 1, 1, 1, 1, 1]*w1;
gp[:,3]=[-1,-1,-1,-1, 1, 1, 1, 1]*g1; w[:,3]=[ 1, 1, 1, 1, 1, 1, 1, 1]*w1;

wp=w[:,1]*w[:,2]*w[:,3];
xsi=gp[:,1];  eta=gp[:,2]; zet=gp[:,3];  r2=ngp*3;

N[:,1]=(1-xsi)*(1-eta)*(1-zet)/8;  N[:,5]=(1-xsi)*(1-eta)*(1+zet)/8;
N[:,2]=(1+xsi)*(1-eta)*(1-zet)/8;  N[:,6]=(1+xsi)*(1-eta)*(1+zet)/8;
N[:,3]=(1+xsi)*(1+eta)*(1-zet)/8;  N[:,7]=(1+xsi)*(1+eta)*(1+zet)/8;
N[:,4]=(1-xsi)*(1+eta)*(1-zet)/8;  N[:,8]=(1-xsi)*(1+eta)*(1+zet)/8;

calc = np.inv( (np.transpose(N) * N) ) * np.transpose(N);

conductivity = np.zeros([nel,8])

for i in range(nel):
#for i = (1:nel):
    #ns[:,0,i] = calc*es[:,0,i];
    #ns[:,1,i] = calc*es[:,1,i];
    #ns[:,2,i] = calc*es[:,2,i];
    conductivity[i,0] = np.sqrt(calc*es[0,0,i]**2 + calc*es[0,1,i]**2 + calc*es[0,2,i]**2)
    conductivity[i,1] = np.sqrt(calc*es[1,0,i]**2 + calc*es[1,1,i]**2 + calc*es[1,2,i]**2)
    conductivity[i,0] = np.sqrt(calc*es[2,0,i]**2 + calc*es[2,1,i]**2 + calc*es[2,2,i]**2)
    conductivity[i,1] = np.sqrt(calc*es[3,0,i]**2 + calc*es[3,1,i]**2 + calc*es[3,2,i]**2)
    conductivity[i,0] = np.sqrt(calc*es[4,0,i]**2 + calc*es[4,1,i]**2 + calc*es[4,2,i]**2)
    conductivity[i,1] = np.sqrt(calc*es[5,0,i]**2 + calc*es[5,1,i]**2 + calc*es[5,2,i]**2)
    conductivity[i,0] = np.sqrt(calc*es[6,0,i]**2 + calc*es[6,1,i]**2 + calc*es[6,2,i]**2)
    conductivity[i,1] = np.sqrt(calc*es[7,0,i]**2 + calc*es[7,1,i]**2 + calc*es[7,2,i]**2)
"""


#cfvv.draw_geometry(edof,coord,dof,ed,3,scale=0.002)


#cfvv.beam3d.draw_displaced_geometry(edof,coord,dof,a,def_scale=5)

"""
# Send data of undeformed geometry
cfvv.beam3d.geometry(edof,coord,dof,0.02,0.2)

eo = np.array([0, 0, 1])

E = 210000000
v = 0.3
G = E/(2*(1+v))

#HEA300-beams
A = 11250*0.000001      # m^2
A_web = 2227*0.000001   # m^2
Iy = 63.1*0.000001      # m^4
hy = 0.29*0.5           # m
Iz = 182.6*0.000001     # m^4
hz = 0.3*0.5            # m
Kv = 0.856*0.000001     # m^4

ep = [E, G, A, Iy, Iz, Kv]

#eq = np.zeros([nel,4])
#eq[23] = [0,-2000,0,0]
#eq[24] = [0,-2000,0,0]
#eq[25] = [0,-2000,0,0]

K = np.zeros([ndof,ndof])
f = np.zeros([ndof,1])

for i in range(nel):
    #Ke,fe = cfc.beam3e(ex[i], ey[i], ez[i], eo, ep, eq[i])
    Ke = cfc.beam3e(ex[i], ey[i], ez[i], eo, ep)
    
    K = cfc.assem(edof[i],K,Ke)
    #f = cfc.assem(edof[i],f,fe)




f[7,0] = -3000
f[13,0] = -3000
f[19,0] = -3000
f[49,0] = -3000
f[55,0] = -3000
f[61,0] = -3000

#f[19,0] = -3000
#f[37,0] = -3000
#f[20,0] = 3000 # Punktlast i z-led
#f[38,0] = 500 # Punktlast i z-led
#f[61,0] = -3000
#f[79,0] = -3000

bcPrescr = np.array([1, 2, 3, 4, 5, 6, 19, 22, 23, 24, 37, 40, 41, 42, 43, 44, 45, 46, 47, 48, 61, 64, 65, 66, 79, 82, 83, 84])
a,r = cfc.solveq(K, f, bcPrescr)

ed = cfc.extractEldisp(edof,a)

# Number of points along the beam
nseg=13 # 13 points in a 3m long beam = 250mm long segments


es = np.zeros((nel*nseg,6))
edi = np.zeros((nel*nseg,4))
eci = np.zeros((nel*nseg,1))

for i in range(29):
    #es[nseg*i:nseg*i+nseg,:], edi[nseg*i:nseg*i+nseg,:], eci[nseg*i:nseg*i+nseg,:] = cfc.beam3s(ex[i],ey[i],ez[i],eo,ep,ed[i],eq[i],nseg)
    es[nseg*i:nseg*i+nseg,:], edi[nseg*i:nseg*i+nseg,:], eci[nseg*i:nseg*i+nseg,:] = cfc.beam3s(ex[i],ey[i],ez[i],eo,ep,ed[i],[0,0,0,0],nseg)
    #es[nseg*i:nseg*i+nseg,:] = cfc.beam3s(ex[i],ey[i],ez[i],eo,ep,ed[i],[0,0,0,0],nseg)

N = es[:,0]
Vy = es[:,1]
Vz = es[:,2]
T = es[:,3]
My = es[:,4]
Mz = es[:,5]

normal_stresses = np.zeros((nel*nseg,1))
shear_stresses_y = np.zeros((nel*nseg,1))
shear_stresses_z = np.zeros((nel*nseg,1))

# Stress calculation based on element forces
for i in range(nel*nseg):
    # Calculate least favorable normal stress using Navier's formula
    if N[i] < 0:
        normal_stresses[i] = N[i]/A - np.absolute(My[i]/Iy*hz) - np.absolute(Mz[i]/Iz*hy)
        #normal_stresses[i] = N[i]/A - My[i]/Iy*hz - Mz[i]/Iz*hy
    else:
        normal_stresses[i] = N[i]/A + np.absolute(My[i]/Iy*hz) + np.absolute(Mz[i]/Iz*hy)
        #normal_stresses[i] = N[i]/A + My[i]/Iy*hz + Mz[i]/Iz*hy

    # Calculate shear stress in y-direction (Assuming only web taking shear stresses)
    shear_stresses_y[i] = Vy[i]/A_web

    # Calculate shear stress in y-direction (Assuming only flanges taking shear stresses)
    shear_stresses_z[i] = Vz[i]/(A-A_web)

    
# Below the data for the undeformed mesh is sent, along with element values.
# Normal stresses are sent by default, but comment it out and uncomment 
# shear_stresses_y/shear_stresses_z to visualize them

# Send data of deformed geometry & normal stresses as element values
cfvv.beam3d.def_geometry(edof,coord,dof,a,normal_stresses,'Max normal stress',def_scale=5,nseg=nseg)

# Send data of deformed geometry & normal stresses as element values
#cfvv.beam3d.def_geometry(edof,coord,dof,a,shear_stresses_y,'Shear stress y',def_scale=5,nseg=nseg)

# Send data of deformed geometry & normal stresses as element values
#cfvv.beam3d.def_geometry(edof,coord,dof,a,shear_stresses_z,'Shear stress z',def_scale=5,nseg=nseg)

"""

