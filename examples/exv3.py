# -*- coding: utf-8 -*-
"""
3D example using Vedo, flow elements

@author: Andreas Åmand
"""

import numpy as np
import calfem.core as cfc
import calfem.vedo as cfvv
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.utils as cfu

l = 1
h = 0.4
w = 0.4

n_el_x = 4
n_el_y = 4
n_el_z = 10

g = cfg.Geometry()

g.point([0, 0, 0], 0)
g.point([w, 0, 0], 1)
g.point([w, l, 0], 2)
g.point([0, l, 0], 3)
g.point([0, 0, h], 4)
g.point([w, 0, h], 5)
g.point([w, l, h], 6)
g.point([0, l, h], 7)

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

marker_bottom = 40
marker_top = 41
marker_fixed_left = 42
marker_back = 43
marker_fixed_right = 44
marker_front = 45

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

coord, edof, dof, bdof, elementmarkers = cfm.mesh(g, el_type, elSizeFactor, dofs_per_node)
ex, ey, ez = cfc.coordxtr(edof, coord, dof)

nnode = np.size(coord, axis = 0)
ndof = np.size(dof, axis = 0)*np.size(dof, axis = 1)
nel = np.size(edof, axis = 0)

Lambda = 1.4 # Thermal conductivity for concrete
D = np.identity(3)*Lambda

K = np.zeros((ndof,ndof))

f = np.zeros([ndof,1])
eq = np.zeros([nel,1])

eq_els = np.array([[70],[89]])

eq[eq_els[0]] = 30000
eq[eq_els[1]] = -30000

for eltopo, elx, ely, elz, eqs in zip(edof, ex, ey, ez, eq):
    Ke,fe = cfc.flw3i8e(elx, ely, elz, [2], D, eqs)
    cfc.assem(eltopo, K, Ke, f, fe)

bc = np.array([],'i')
bcVal = np.array([],'i')

bc, bcVal = cfu.apply_bc_3d(bdof, bc, bcVal, marker_fixed_left, 20.0) # bottom is 0
bc, bcVal = cfu.apply_bc_3d(bdof, bc, bcVal, marker_fixed_right, 20.0) # top is 20

T,r = cfc.solveq(K, f, bc, bcVal)

ed = cfc.extract_eldisp(edof,T)

es = np.zeros((8,3,nel))
edi = np.zeros((8,3,nel))
eci = np.zeros((8,3,nel))

for i in range(nel):
    es[0:8,:,i], edi[0:8,:,i], eci[0:8,:,i] = cfc.flw3i8s(ex[i],ey[i],ez[i],[2],D,ed[i])

vectors = np.zeros((nel,3));
flux = np.zeros((nel,3));
flux_tot = np.zeros((nel,1));
for i in range(nel):
    flux[i,:] = [np.average([es[:,0,i]]), np.average([es[:,1,i]]), np.average([es[:,2,i]])]
    flux_tot[i] = np.sqrt(flux[i,0]**2 + flux[i,1]**2 + flux[i,2]**2)

# For simple point coordinates of geometry
points = g.points

# For lines with point connectivity
lines = g.curves

# For surfaces with line connectivity
surfaces = g.surfaces

cfv.draw_geometry(points,lines,surfaces)

cfv.figure(2)
cfv.draw_mesh(edof,coord,dof,3,alpha=1,scale=0.005,bcPrescr=bc, bc=bcVal, eq_els=eq_els, eq=eq[eq_els])

disp = np.zeros((nnode,1))

cfv.figure(3)
cfv.draw_displaced_mesh(edof,coord,dof,3,scalars=flux_tot,colormap='coolwarm',wireframe=True)
cfv.add_axes(xrange=[-0.1,0.5], yrange=[-0.1,1.1], zrange=[-0.1,0.5])

cfv.add_scalar_bar('Heat flux [W/m^2]')
cfv.elflux(ex,ey,ez,flux,colormap='coolwarm')

cfv.figure(4)
mesh = cfv.draw_displaced_mesh(edof,coord,dof,3,scalars=ed,colormap='coolwarm',colors=5,scalar_title='Temp. [C]')
cfv.add_mesh_axes(mesh)
cfv.add_scalar_bar('Temp. [C]')


# Animation

cfv.figure(5)
cfv.animation(edof,coord,dof,3,scalars=ed,dt=250,steps=20,colormap='coolwarm',colors=5,export=True,file='export/exv3/anim/exv3',scalar_title='Temp. [C]')
cfv.add_scalar_bar('Temp. [C]')

#Start Calfem-vedo visualization
cfv.show_and_wait()

# Export the mesh
cfv.export_vtk('export/exv3/exv3', mesh)

# For not exporting animation
#cfv.animation(edof,coord,dof,3,scalars=ed,dt=250,steps=20,colormap='coolwarm',colors=5)

