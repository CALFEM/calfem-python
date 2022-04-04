# -*- coding: utf-8 -*-
"""
3D example using Vedo, solid elements

@author: Andreas Åmand
"""

import numpy as np
import calfem.core as cfc
import calfem.vedo as cfv
import calfem.vedo_utils as cfvu

edof,coord,dof,a,ed,bc,f_dofs,Stress_tensors,vM_el,vM_n,lamb,eig = cfvv.import_mat('exv4',['edof','coord','dof','a','ed','bc','force_dofs','Stress_tensors','vM_el','vM_n','lambda','eig'])

ex,ey,ez = cfc.coordxtr(edof,coord,dof)

eigenmode = 0 # Choose what eigenmode to display in figure 5/6

ndof = np.size(dof, axis = 0)*np.size(dof, axis = 1)
ncoord = np.size(coord, axis = 0)
nel = np.size(edof, axis = 0)

mode_a = np.zeros((nel, 1))
tot_deform = np.zeros(8)
for i in range(nel):
	coords = cfvu.get_coord_from_edof(edof[i,:],dof,4)
	for j in range(8):
		deform = cfvu.get_a_from_coord(coords[j],3,eig[:,eigenmode])
		tot_deform[j] = np.sqrt(deform[0]**2 + deform[1]**2 + deform[2]**2)

	mode_a[i,:] = np.average(tot_deform)

Freq=np.sqrt(lamb[eigenmode]/(2*np.pi))

''' Principal stresses '''

ps_val = np.zeros((nel,3))
ps_vec = np.zeros((nel,3,3))
for i in range(nel):
    ps_val[i,:], ps_vec[i,:,:] = np.linalg.eig(Stress_tensors[:,:,i])



upd_ed = np.zeros((nel,8))
for i in range(nel):
    upd_ed[i,0] = np.sqrt( ed[i,0]**2 + ed[i,1]**2 + ed[i,2]**2 )
    upd_ed[i,1] = np.sqrt( ed[i,3]**2 + ed[i,4]**2 + ed[i,5]**2 )

    upd_ed[i,2] = np.sqrt( ed[i,6]**2 + ed[i,7]**2 + ed[i,8]**2 )
    upd_ed[i,3] = np.sqrt( ed[i,9]**2 + ed[i,10]**2 + ed[i,11]**2 )

    upd_ed[i,4] = np.sqrt( ed[i,12]**2 + ed[i,13]**2 + ed[i,14]**2 )
    upd_ed[i,5] = np.sqrt( ed[i,15]**2 + ed[i,16]**2 + ed[i,17]**2 )

    upd_ed[i,6] = np.sqrt( ed[i,18]**2 + ed[i,19]**2 + ed[i,20]**2 )
    upd_ed[i,7] = np.sqrt( ed[i,21]**2 + ed[i,22]**2 + ed[i,23]**2 )

bcPrescr = bc
bc = np.zeros((np.size(bc[:,0]),1))
f = -5000*np.ones((np.size(f_dofs[:,0]),1))

# First plot, undeformed mesh
cfv.figure(1)

cfv.draw_mesh(edof,coord,dof,4,scale=0.005,bcPrescr=bcPrescr[:,0],bc=bc[:,0],fPrescr=f_dofs[:,0],f=f[:,0])
cfv.add_text('Undeformed mesh + Forces & BCs for static analysis')
cfv.show_and_wait()

# Second plot, deformed mesh with element stresses
cfv.figure(2)

scalefact = 3 #deformation scale factor
static = cfv.draw_displaced_mesh(edof,coord,dof,4,a,vM_el/1000000,def_scale=scalefact,scalar_title='von Mises [MPa]')

cfv.add_text('Static analysis: self-weight & ecc. vertical load', pos='top-left')
cfv.add_text(f'Deformation scalefactor: {scalefact}',pos='top-right')
cfv.add_scalar_bar('von Mises [MPa]')
cfv.show_and_wait()

# Third plot, animation of figure 2
cfv.figure(3)
scalefact = 3 #deformation scale factor

cfv.animation(edof,coord,dof,4,a,vM_el/1000000,def_scale=scalefact)

cfv.add_text('Static analysis: self-weight & ecc. vertical load', pos='top-left')
cfv.add_text(f'Deformation scalefactor: {scalefact}',pos='top-right')

cfv.add_scalar_bar('von Mises [MPa]')

#Start Calfem-vedo visualization
cfv.show_and_wait()

# Fourth plot, principal stresses for static analysis
cfv.figure(4)

# Return the mesh for export
cfv.draw_displaced_mesh(edof,coord,dof,4,a,upd_ed*1000,wireframe=True)
cfv.elprinc(ex,ey,ez,ps_val/1000000,ps_vec,ed,colormap='coolwarm',unit='MPa')
cfv.add_scalar_bar('Deformation [mm]')
cfv.add_text('Static analysis',pos='top-left')
cfv.add_text('Deformation scalefactor: 1',pos='top-right')
cfv.add_text('Princ. stress vectors',pos='top-middle')

#Start Calfem-vedo visualization
cfv.show_and_wait()

# Fifth plot, first mode from eigenvalue analysis
cfv.figure(5)

scalefact = 100 #deformation scale factor
modal = cfv.draw_displaced_mesh(edof,coord,dof,4,eig[:,eigenmode],mode_a*1000,def_scale=scalefact,lines=True,scalar_title='Tot. el. displacement [mm]')
cfv.add_text(f'Modal analysis: {eigenmode+1}st mode',pos='top-left')
cfv.add_text(f'Frequency: {round(Freq[0],2)} Hz')
cfv.add_text(f'Deformation scalefactor: {scalefact}',pos='top-right')
cfv.add_scalar_bar('Tot. el. displacement [mm]')
cfv.add_projection(plane='xz',rulers=True)
cfv.show_and_wait()

# Sixth plot, animation of figure 5
cfv.figure(6)

cfv.add_text(f'Modal analysis: {eigenmode+1}st mode',pos='top-left')
cfv.add_text(f'Frequency: {round(Freq[0],2)} Hz')
cfv.add_text(f'Deformation scalefactor: {scalefact}',pos='top-right')

cfv.animation(edof,coord,dof,4,eig[:,eigenmode],mode_a*1000,def_scale=scalefact,negative=True,scalar_title='Tot. el. displacement [mm]',export=True,file='export/exv4/anim/exv4_modal')

cfv.add_scalar_bar('Tot. el. displacement [mm]')

#Start Calfem-vedo visualization
cfv.show_and_wait()

# Export the two meshes
cfv.export_vtk('export/exv4/exv4_static', static)
cfv.export_vtk('export/exv4/exv4_modal', modal)

# For not exporting animation
#cfv.animation(edof,coord,dof,4,eig[:,eigenmode],mode_a*1000,def_scale=scalefact,negative=True)

