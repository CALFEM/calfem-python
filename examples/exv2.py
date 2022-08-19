# -*- coding: utf-8 -*-
"""
3D example using Vedo, bar & beam elements

@author: Andreas Åmand
"""

import numpy as np
import calfem.core as cfc
import calfem.vedo as cfv

coord = np.array([
    [0, 0, 0],
    [3, 0, 0],
    [6, 0, 0],
    [9, 0, 0],
    [3, 3, 0],
    [6, 3, 0],
    [9, 3, 0],
    [0, 0, 3],
    [3, 0, 3],
    [6, 0, 3],
    [9, 0, 3],
    [3, 3, 3],
    [6, 3, 3],
    [9, 3, 3]
])

dof = np.array([
    [1, 2, 3, 4, 5, 6],
    [7, 8, 9, 10, 11, 12],
    [13, 14, 15, 16, 17, 18],
    [19, 20, 21, 22, 23, 24],
    [25, 26, 27, 28, 29, 30],
    [31, 32, 33, 34, 35, 36],
    [37, 38, 39, 40, 41, 42],
    [43, 44, 45, 46, 47, 48],
    [49, 50, 51, 52, 53, 54],
    [55, 56, 57, 58, 59, 60],
    [61, 62, 63, 64, 65, 66],
    [67, 68, 69, 70, 71, 72],
    [73, 74, 75, 76, 77, 78],
    [79, 80, 81, 82, 83, 84]
])

edof_beams = np.array([
    '''Left side'''
    [1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12],
    [7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18],
    [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24],
    [1,  2,  3,  4,  5,  6,  25, 26, 27, 28, 29, 30],
    [7,  8,  9,  10, 11, 12, 25, 26, 27, 28, 29, 30],
    [13, 14, 15, 16, 17, 18, 31, 32, 33, 34, 35, 36],
    [19, 20, 21, 22, 23, 24, 37, 38, 39, 40, 41, 42],
    [25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36],
    [31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
    '''Right side'''
    [43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54],
    [49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60],
    [55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66],
    [43, 44, 45, 46, 47, 48, 67, 68, 69, 70, 71, 72],
    [49, 50, 51, 52, 53, 54, 67, 68, 69, 70, 71, 72],
    [55, 56, 57, 58, 59, 60, 73, 74, 75, 76, 77, 78],
    [61, 62, 63, 64, 65, 66, 79, 80, 81, 82, 83, 84],
    [67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78],
    [73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84],
    '''Bottom'''
    [1,  2,  3,  4,  5,  6,  43, 44, 45, 46, 47, 48],
    [7,  8,  9,  10, 11, 12, 49, 50, 51, 52, 53, 54],
    [13, 14, 15, 16, 17, 18, 55, 56, 57, 58, 59, 60],
    [19, 20, 21, 22, 23, 24, 61, 62, 63, 64, 65, 66]
])

nnode = np.size(coord, axis = 0)
ndof = np.size(dof, axis = 0)*np.size(dof, axis = 1)
nel_beams = np.size(edof_beams, axis = 0)

ex_beams,ey_beams,ez_beams = cfc.coordxtr(edof_beams,coord,dof)

eo = np.array([
    [0, 0, 1],
    [0, 0, 1],
    [0, 0, 1],
    [0, 0, 1],
    [0, 0, 1],
    [0, 0, 1],
    [0, 0, 1],
    [0, 0, 1],
    [0, 0, 1],

    [0, 0, 1],
    [0, 0, 1],
    [0, 0, 1],
    [0, 0, 1],
    [0, 0, 1],
    [0, 0, 1],
    [0, 0, 1],
    [0, 0, 1],
    [0, 0, 1],

    [-1, 0, 0],
    [-1, 0, 0],
    [-1, 0, 0],
    [-1, 0, 0]
    ])

E = 210000000               # Pa
v = 0.3
G = E/(2*(1+v))             # Pa

#HEA300-beams
A_beams = 11250*0.000001    # m^2
A_web_beams = 2227*0.000001 # m^2
Iy = 63.1*0.000001          # m^4
hy = 0.29*0.5               # m
Iz = 182.6*0.000001         # m^4
hz = 0.3*0.5                # m
Kv = 0.856*0.000001         # m^4
ep_beams = [E, G, A_beams, Iy, Iz, Kv]

K = np.zeros([ndof,ndof])
f = np.zeros([ndof,1])

# 1 kN point load in z-direction at symmetry (top)
f[38] = 1000 # N
# Transverse beams at y=0
# Dead load of 3x3 m concrete, t = 0.2 m, 25 kN/m^2
q1 = 25*1000*0.2*3 #N/m
eq1 = [0,-q1,0,0]



for i in range(nel_beams):
    if i < 18:
        Ke = cfc.beam3e(ex_beams[i], ey_beams[i], ez_beams[i], eo[i], ep_beams)
        K = cfc.assem(edof_beams[i],K,Ke)
    else:
        Ke, fe = cfc.beam3e(ex_beams[i], ey_beams[i], ez_beams[i], eo[i], ep_beams,eq1)
        K, f = cfc.assem(edof_beams[i],K,Ke,f,fe)

edof_bars = np.array([
    [13, 14, 15, 25, 26, 27],
    [19, 20, 21, 31, 32, 33],
    [55, 56, 57, 67, 68, 69],
    [61, 62, 63, 73, 74, 75],
    [25, 26, 27, 67, 68, 69],
    [31, 32, 33, 73, 74, 75],
    [37, 38, 39, 79, 80, 81]
])

nel_bars = np.size(edof_bars, axis = 0)

ex_bars = np.array([
    [6,3],
    [9,6],
    [6,3],
    [9,6],
    [3,3],
    [6,6],
    [9,9],
])

ey_bars = np.array([
    [0,3],
    [0,3],
    [0,3],
    [0,3],
    [3,3],
    [3,3],
    [3,3],
])

ez_bars = np.array([
    [0,0],
    [0,0],
    [3,3],
    [3,3],
    [0,3],
    [0,3],
    [0,3],
])

#Solid 100mm circular bars
A_bars = 0.05*np.pi*np.pi   # m^2

ep_bars = [E, A_bars]

for i in range(nel_bars):
    Ke = cfc.bar3e(ex_bars[i], ey_bars[i], ez_bars[i], ep_bars)
    K = cfc.assem(edof_bars[i],K,Ke)

bcPrescr = np.array([1, 2, 3, 4, 5, 6, 19, 22, 23, 24, 37, 40, 41, 42, 43, 44, 45, 46, 47, 48, 61, 64, 65, 66, 79, 82, 83, 84])
bcVal = np.zeros([1,28])
a,r = cfc.solveq(K, f, bcPrescr)

ed_beams = cfc.extractEldisp(edof_beams,a)

# Number of points along the beam
nseg=13 # 13 points in a 3m long beam = 250mm long segments

es_beams = np.zeros((nel_beams*nseg,6))
edi_beams = np.zeros((nel_beams*nseg,4))
eci_beams = np.zeros((nel_beams*nseg,1))

for i in range(nel_beams):
    if i < 18:
        es_beams[nseg*i:nseg*i+nseg,:], edi_beams[nseg*i:nseg*i+nseg,:], eci_beams[nseg*i:nseg*i+nseg,:] = cfc.beam3s(ex_beams[i],ey_beams[i],ez_beams[i],eo[i],ep_beams,ed_beams[i],[0,0,0,0],nseg)
    else:
        es_beams[nseg*i:nseg*i+nseg,:], edi_beams[nseg*i:nseg*i+nseg,:], eci_beams[nseg*i:nseg*i+nseg,:] = cfc.beam3s(ex_beams[i],ey_beams[i],ez_beams[i],eo[i],ep_beams,ed_beams[i],eq1,nseg)

N_beams = es_beams[:,0]
Vy = es_beams[:,1]
Vz = es_beams[:,2]
T = es_beams[:,3]
My = es_beams[:,4]
Mz = es_beams[:,5]

ed_bars = cfc.extractEldisp(edof_bars,a)

es_bars = np.zeros((nel_bars,1))

for i in range(nel_bars):
    es_bars[i,:] = cfc.bar3s(ex_bars[i],ey_bars[i],ez_bars[i],ep_bars,ed_bars[i])

N_bars = es_bars

normal_stresses_beams = np.zeros(nel_beams*nseg)

# Stress calculation based on element forces
for i in range(nel_beams*nseg):
    if N_beams[i] < 0:
        normal_stresses_beams[i] = N_beams[i]/A_beams - np.absolute(My[i]/Iy*hz) - np.absolute(Mz[i]/Iz*hy)
    else:
        normal_stresses_beams[i] = N_beams[i]/A_beams + np.absolute(My[i]/Iy*hz) + np.absolute(Mz[i]/Iz*hy)

normal_stresses_bars = np.zeros(nel_bars)

for i in range(nel_bars):
    # Calculate least favorable normal stress
    normal_stresses_bars[i] = N_bars[i]/A_bars

eq_els = np.array([[18],[19],[20],[21]])

eq = np.zeros([nel_beams,4])

# For color-coding elements with force allpied
eq[eq_els[0,:]] = eq1
eq[eq_els[1,:]] = eq1
eq[eq_els[2,:]] = eq1
eq[eq_els[3,:]] = eq1

bcPrescr = np.transpose(bcPrescr)

cfv.draw_mesh(edof_beams,coord,dof,5,nseg=nseg,alpha=0.2, eq_els=eq_els, eq=eq[eq_els])
beams = cfv.draw_displaced_mesh(edof_beams,coord,dof,5,a,normal_stresses_beams/1000000,nseg=nseg,scalar_title='Max normal stress [MPa]')

### IMPORTANT: Update to 12 dofs here, otherwise visualization breaks
### In reality, 3D-bars only have 6 dofs
### Alternative solution is to have separate dof-matrices (dof_beams & dof_bars)
edof_bars = np.array([
    [13, 14, 15, 16, 17, 18, 25, 26, 27, 28, 29, 30],#[13, 14, 15, 25, 26, 27],
    [19, 20, 21, 22, 23, 24, 31, 32, 33, 34, 35, 36],#[19, 20, 21, 31, 32, 33],
    [55, 56, 57, 58, 59, 60, 67, 68, 69, 70, 71, 72],#[55, 56, 57, 67, 68, 69],
    [61, 62, 63, 64, 65, 66, 73, 74, 75, 76, 77, 78],#[55, 56, 57, 73, 74, 75],
    [25, 26, 27, 28, 29, 30, 67, 68, 69, 70, 71, 72],#[25, 26, 27, 67, 68, 69],
    [31, 32, 33, 34, 35, 36, 73, 74, 75, 76, 77, 78],#[31, 32, 33, 73, 74, 75],
    [37, 38, 39, 40, 41, 42, 79, 80, 81, 82, 83, 84]#[37, 38, 39, 79, 80, 81]
])

cfv.draw_mesh(edof_bars,coord,dof,2,alpha=0.2)
vmin, vmax = np.min(normal_stresses_beams), np.max(normal_stresses_beams)
bars = cfv.draw_displaced_mesh(edof_bars,coord,dof,2,a,normal_stresses_bars/1000000,vmin=vmin,vmax=vmax,scalar_title='Max normal stress [MPa]')

Mz_beams = np.zeros((4,nseg))
eci_beams_upd = np.zeros((4,nseg))
for i in range(Mz_beams.shape[0]):
    Mz_beams[i] = Mz[(18+i)*nseg:(18+i)*nseg+nseg]
    eci_beams_upd[i] = np.transpose(eci_beams[(18+i)*nseg:(18+i)*nseg+nseg])

cfv.eldia(ex_beams[18:22],ey_beams[18:22],ez_beams[18:22],Mz_beams/1000,eci_beams_upd,label='M_x [kNm]')
cfv.add_rulers()

cfv.figure(2)
steps = 20
cfv.add_text(f'Looping from undef. to def. state w/ {steps} steps',pos='top-middle')
cfv.animation(edof_beams,coord,dof,5,a,normal_stresses_beams/1000000,nseg=nseg,dt=125,steps=20)
cfv.animation(edof_bars,coord,dof,2,a,normal_stresses_bars/1000000,nseg=nseg,dt=125,steps=20,vmax=vmax,vmin=vmin)
cfv.add_scalar_bar('Max normal stress [MPa]')

#Start Calfem-vedo visualization
cfv.show_and_wait()

# For exporting meshes
cfv.export_vtk('export/exv2/exv2_beams', [beams,bars])

