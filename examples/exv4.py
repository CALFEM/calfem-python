# -*- coding: utf-8 -*-
"""
3D example using Vedo, solid elements

@author: Andreas Åmand
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
import vedo_utils as cfvu
#from scipy.io import loadmat

#edof,coord,dof,a,es,ns,lamb,eig = cfvv.import_mat('exv4',['edof','coord','dof','a','es','ns','lambda','eig'])
edof,coord,dof,a,ed,bc,f_dofs,Stress_tensors,vM_el,vM_n,lamb,eig = cfvv.import_mat('exv4',['edof','coord','dof','a','ed','bc','force_dofs','Stress_tensors','vM_el','vM_n','lambda','eig'])

ex,ey,ez = cfc.coordxtr(edof,coord,dof)

eigenmode = 0 # Choose what eigenmode to display in figure 4/5

ndof = np.size(dof, axis = 0)*np.size(dof, axis = 1)
ncoord = np.size(coord, axis = 0)
print('Number of CALFEM coordinates: ',ncoord)
nel = np.size(edof, axis = 0)

''' Principal stresses '''

ps_val = np.zeros((nel,3))
ps_vec = np.zeros((nel,3,3))
for i in range(nel):
    ps_val[i,:], ps_vec[i,:,:] = np.linalg.eig(Stress_tensors[:,:,i])



#print('Principal stresses el. 1')
#print(PS_val,PS_vec)
#print('---')









#print(data)

#solid_data = loadmat('exv4.mat')

#edof = solid_data['edof']
#edof = np.delete(edof,0,1)
#coord = solid_data['coord']
#dof = solid_data['dof']
#a = solid_data['a']
#ed = solid_data['ed']
#es = solid_data['es']
#et = solid_data['et']
#eci = solid_data['eci']
#ns = solid_data['ns']
#L = solid_data['L']
#X = solid_data['X']






#inter = np.zeros((ncoord, 1))

#for i in range(ncoord):
	#inter[i,:] = X[3*i,0]*X[3*i,0] + X[3*i,1]*X[3*i,1] + X[3*i,2]*X[3*i,2]

mode_a = np.zeros((nel, 1))
tot_deform = np.zeros(8)
for i in range(nel):
	coords = cfvu.get_coord_from_edof(edof[i,:],dof,4)
	#print(coords)
	#eig[coords,0]
	for j in range(8):
		deform = cfvu.get_a_from_coord(coords[j],3,eig[:,eigenmode])
		tot_deform[j] = np.sqrt(deform[0]**2 + deform[1]**2 + deform[2]**2)

	mode_a[i,:] = np.average(tot_deform)
	#print(x1)
	"""
	mode_a[i,:] = np.sqrt(
		X[coords[0],0]*X[coords[0],0] + 
		X[coords[1],0]*X[coords[1],0] + 
		X[coords[2],0]*X[coords[2],0] + 
		X[coords[3],0]*X[coords[3],0] + 
		X[coords[4],0]*X[coords[4],0] + 
		X[coords[5],0]*X[coords[5],0] + 
		X[coords[6],0]*X[coords[6],0] + 
		X[coords[7],0]*X[coords[7],0]
		)
	"""
#print(np.size(mode_a, axis = 0))
#print(np.size(mode_a, axis = 1))
#print(np.size(X, axis = 0))
#print(np.size(X, axis = 1))
#mode_a[:,0] = X[:,0]
#print(mode_a)

Freq=np.sqrt(lamb[eigenmode]/(2*np.pi))
#print(L)
#print(np.size(Freq, axis = 1))

#print(np.size(X[0,:], axis = 0))
#print(np.size(X[0,:], axis = 1))
#print(X[:,0])

#print(np.size(a, axis = 0))
#print(np.size(a, axis = 1))
#print(a)

upd_ed = np.zeros((nel,8))
#ed = np.asarray(ed)
for i in range(nel):
    print(ed[i,0]**2 + ed[i,1]**2, ed[i,2]**2)
    print(np.sqrt( ed[i,0]**2 + ed[i,1]**2 + ed[i,2]**2 ))
    upd_ed[i,0] = np.sqrt( ed[i,0]**2 + ed[i,1]**2 + ed[i,2]**2 )
    upd_ed[i,1] = np.sqrt( ed[i,3]**2 + ed[i,4]**2 + ed[i,5]**2 )

    upd_ed[i,2] = np.sqrt( ed[i,6]**2 + ed[i,7]**2 + ed[i,8]**2 )
    upd_ed[i,3] = np.sqrt( ed[i,9]**2 + ed[i,10]**2 + ed[i,11]**2 )

    upd_ed[i,4] = np.sqrt( ed[i,12]**2 + ed[i,13]**2 + ed[i,14]**2 )
    upd_ed[i,5] = np.sqrt( ed[i,15]**2 + ed[i,16]**2 + ed[i,17]**2 )

    upd_ed[i,6] = np.sqrt( ed[i,18]**2 + ed[i,19]**2 + ed[i,20]**2 )
    upd_ed[i,7] = np.sqrt( ed[i,21]**2 + ed[i,22]**2 + ed[i,23]**2 )

"""
stresses = np.zeros([nel,6])

#print(nel)

#print(np.average(es[0,:,0]))
#print(np.average(es[:,0,8]))

for i in range(0, nel):
	#print(i)
	stresses[i,0] = np.average(es[:,0,i])
	stresses[i,1] = np.average(es[:,1,i])
	stresses[i,2] = np.average(es[:,2,i])
	stresses[i,3] = np.average(es[:,3,i])
	stresses[i,4] = np.average(es[:,4,i])
	stresses[i,5] = np.average(es[:,5,i])

#print(stresses)

von_mises_elements = np.zeros([nel,1])

for i in range(0, nel):
	#print(i)
	von_mises_elements[i] = np.sqrt( 0.5 * ( np.square(stresses[i,0]-stresses[i,1]) + np.square(stresses[i,1]-stresses[i,2]) + np.square(stresses[i,2]-stresses[i,0]) ) + 3 * (np.square(stresses[i,3]) + np.square(stresses[i,4]) + np.square(stresses[i,5])) )

#print(von_mises)

von_mises_nodes = np.zeros([nel,8])



for i in range(0, nel):
	#print(i)
	von_mises_nodes[i,0] = np.sqrt( 0.5 * ( np.square(ns[0,0,i]-ns[0,1,i]) + np.square(ns[0,1,i]-ns[0,2,i]) + np.square(ns[0,2,i]-ns[0,0,i]) ) + 3 * (np.square(ns[0,3,i]) + np.square(ns[0,4,i]) + np.square(ns[0,5,i])) )
	von_mises_nodes[i,1] = np.sqrt( 0.5 * ( np.square(ns[1,0,i]-ns[1,1,i]) + np.square(ns[1,1,i]-ns[1,2,i]) + np.square(ns[1,2,i]-ns[1,0,i]) ) + 3 * (np.square(ns[1,3,i]) + np.square(ns[1,4,i]) + np.square(ns[1,5,i])) )
	von_mises_nodes[i,2] = np.sqrt( 0.5 * ( np.square(ns[2,0,i]-ns[2,1,i]) + np.square(ns[2,1,i]-ns[2,2,i]) + np.square(ns[2,2,i]-ns[2,0,i]) ) + 3 * (np.square(ns[2,3,i]) + np.square(ns[2,4,i]) + np.square(ns[2,5,i])) )
	von_mises_nodes[i,3] = np.sqrt( 0.5 * ( np.square(ns[3,0,i]-ns[3,1,i]) + np.square(ns[3,1,i]-ns[3,2,i]) + np.square(ns[3,2,i]-ns[3,0,i]) ) + 3 * (np.square(ns[3,3,i]) + np.square(ns[3,4,i]) + np.square(ns[3,5,i])) )
	von_mises_nodes[i,4] = np.sqrt( 0.5 * ( np.square(ns[4,0,i]-ns[4,1,i]) + np.square(ns[4,1,i]-ns[4,2,i]) + np.square(ns[4,2,i]-ns[4,0,i]) ) + 3 * (np.square(ns[4,3,i]) + np.square(ns[4,4,i]) + np.square(ns[4,5,i])) )
	von_mises_nodes[i,5] = np.sqrt( 0.5 * ( np.square(ns[5,0,i]-ns[5,1,i]) + np.square(ns[5,1,i]-ns[5,2,i]) + np.square(ns[5,2,i]-ns[5,0,i]) ) + 3 * (np.square(ns[5,3,i]) + np.square(ns[5,4,i]) + np.square(ns[5,5,i])) )
	von_mises_nodes[i,6] = np.sqrt( 0.5 * ( np.square(ns[6,0,i]-ns[6,1,i]) + np.square(ns[6,1,i]-ns[6,2,i]) + np.square(ns[6,2,i]-ns[6,0,i]) ) + 3 * (np.square(ns[6,3,i]) + np.square(ns[6,4,i]) + np.square(ns[6,5,i])) )
	von_mises_nodes[i,7] = np.sqrt( 0.5 * ( np.square(ns[7,0,i]-ns[7,1,i]) + np.square(ns[7,1,i]-ns[7,2,i]) + np.square(ns[7,2,i]-ns[7,0,i]) ) + 3 * (np.square(ns[7,3,i]) + np.square(ns[7,4,i]) + np.square(ns[7,5,i])) )
"""




#print(von_mises_elements)
#print(von_mises_nodes)


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

# First plot, undeformed mesh
cfvv.figure(1)

bcPrescr = bc
bc = np.zeros((np.size(bc[:,0]),1))
print('bc',bc[:,0])
f = -5000*np.ones((np.size(f_dofs[:,0]),1))
print('f',f[:,0])

print('bcP',bcPrescr[:,0])
print('bc',bc[:,0])
#sys.exit()

cfvv.draw_mesh(edof,coord,dof,4,scale=0.005,bcPrescr=bcPrescr[:,0],bc=bc[:,0],fPrescr=f_dofs[:,0],f=f[:,0])
cfvv.add_text('Undeformed mesh + Forces & BCs for static analysis')
#cfvv.add_rulers()
cfvv.show_and_wait()
#sys.exit()

# Third plot, deformed mesh with element stresses
cfvv.figure(2)

print('Number of von Mises stresses in elements: ',np.size(vM_el))

print(a)

scalefact = 3 #deformation scale factor
#mesh1 = cfvv.test(edof,ex,ey,ez,a,von_mises_elements/1000000,def_scale=scalefact)
#mesh1 = cfvv.draw_displaced_mesh(edof,coord,dof,4,a,von_mises_elements/1000000,def_scale=scalefact)
static = cfvv.draw_displaced_mesh(edof,coord,dof,4,a,vM_el/1000000,def_scale=scalefact,scalar_title='von Mises [MPa]')
#cfvv.add_vectors(ex,ey,ez)
#cfvv.tensors(ex,ey,ez,ps)

cfvv.add_text('Static analysis: self-weight & ecc. vertical load', pos='top-left')
cfvv.add_text(f'Deformation scalefactor: {scalefact}',pos='top-right')
cfvv.add_scalar_bar('von Mises [MPa]')
cfvv.show_and_wait()







cfvv.figure(3)
# Second plot, first mode from eigenvalue analysis
#cfvv.add_text('Eigenvalue analysis: first mode',window=1)
scalefact = 3 #deformation scale factor
#mode_mesh = cfvv.draw_displaced_geometry(edof,coord,dof,4,X[:,0],mode_a,def_scale=scalefact,scale=0.002,render_nodes=False,merge=True)
#cfvv.add_text(f'Deformation scalefactor: {scalefact}',pos='top-left')
#cfvv.add_scalar_bar(mode_mesh,'Tot. el. displacement')

#cfvv.animation(edof,coord,dof,4,eig[:,eigenmode],10,mode_a*1000,def_scale=scalefact,export=True,file='anim/exv4b')



cfvv.animation(edof,coord,dof,4,a,vM_el/1000000,def_scale=scalefact)
#cfvv.animation(edof,coord,dof,4,a,def_scale=scalefact,export=True,file='export/exv4/anim/exv4_static')

# For exporting animation
#cfvv.animation(edof,coord,dof,4,a,vM_el/1000000,def_scale=scalefact,export=True,file='export/exv4/anim/exv4_static',scalar_title='von Mises [MPa]')


#cfvv.animate(edof,coord,dof,4,eig[:,0],10,def_scale=scalefact,export=True)
cfvv.add_text('Static analysis: self-weight & ecc. vertical load', pos='top-left')
cfvv.add_text(f'Deformation scalefactor: {scalefact}',pos='top-right')
#cfvv.add_text('Note: next fig. might take some time to render', pos='bottom-middle')

cfvv.add_scalar_bar('von Mises [MPa]')

#Start Calfem-vedo visualization
cfvv.show_and_wait()









#cfvv.figure(3)
#mesh2 = cfvv.draw_displaced_mesh(edof,coord,dof,4,a,vM_n/1000000,def_scale=scalefact,lines=True)
#cfvv.add_scalar_bar('von Mises [MPa]')
#cfvv.add_text('Static analysis: self-weight & ecc. vertical load',pos='top-left')
#cfvv.export_vtk('exv4a', mesh2)
#cfvv.show_and_wait()


# Fourth plot, deformed mesh with nodal stresses
cfvv.figure(4)

print('Number of von Mises stresses at nodes: ',np.size(vM_n))
print(vM_n[0])
print('Number of element displacements: ',np.size(upd_ed))
print(upd_ed[0])
# Return the mesh for export
#mesh2 = cfvv.test(edof,ex,ey,ez,a,von_mises_nodes/1000000,def_scale=scalefact,merge=True)
#mesh2 = cfvv.draw_displaced_mesh(edof,coord,dof,4,a,von_mises_nodes/1000000,def_scale=scalefact,merge=True)
#mesh2 = cfvv.draw_displaced_mesh(edof,coord,dof,4,a,vM_n/1000000,def_scale=scalefact,merge=True)
cfvv.draw_displaced_mesh(edof,coord,dof,4,a,upd_ed*1000,wireframe=True)
cfvv.elprinc(ex,ey,ez,ps_val/1000000,ps_vec,ed,colormap='coolwarm',unit='MPa')
cfvv.add_scalar_bar('Deformation [mm]')
#cfvv.add_scalar_bar('Stress [MPa]',pos=[0.8,0.65],text_pos='top-right',on='vectors')
cfvv.add_text('Static analysis',pos='top-left')
cfvv.add_text('Deformation scalefactor: 1',pos='top-right')
cfvv.add_text('Princ. stress vectors',pos='top-middle')
#cfvv.add_text(f'Deformation scalefactor: {scalefact}',pos='top-right')
# Export the mesh to 'exv4a.vtk'





#Start Calfem-vedo visualization
cfvv.show_and_wait()
#sys.exit()

# Second plot, first mode from eigenvalue analysis
cfvv.figure(5)

#print(eig[:,0])
#disp = np.zeros((ncoord*3,1))
#disp[:,0] = eig[:,0]

scalefact = 100 #deformation scale factor
#scalefact = 1 #deformation scale factor
#cfvv.test(edof,ex,ey,ez,eigen[:,0],mode_a*1000,def_scale=scalefact)
modal = cfvv.draw_displaced_mesh(edof,coord,dof,4,eig[:,eigenmode],mode_a*1000,def_scale=scalefact,lines=True,scalar_title='Tot. el. displacement [mm]')
cfvv.add_text(f'Modal analysis: {eigenmode+1}st mode',pos='top-left')
cfvv.add_text(f'Frequency: {round(Freq[0],2)} Hz')
cfvv.add_text(f'Deformation scalefactor: {scalefact}',pos='top-right')
cfvv.add_scalar_bar('Tot. el. displacement [mm]')
cfvv.add_projection(plane='xz',rulers=True)
cfvv.show_and_wait()



### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###


cfvv.figure(6)
# Second plot, first mode from eigenvalue analysis
#cfvv.add_text('Eigenvalue analysis: first mode',window=1)
scalefact = 100 #deformation scale factor
#mode_mesh = cfvv.draw_displaced_geometry(edof,coord,dof,4,X[:,0],mode_a,def_scale=scalefact,scale=0.002,render_nodes=False,merge=True)
cfvv.add_text(f'Modal analysis: {eigenmode+1}st mode',pos='top-left')
cfvv.add_text(f'Frequency: {round(Freq[0],2)} Hz')
cfvv.add_text(f'Deformation scalefactor: {scalefact}',pos='top-right')
#cfvv.add_text(f'Deformation scalefactor: {scalefact}',pos='top-left')
#cfvv.add_scalar_bar(mode_mesh,'Tot. el. displacement')

#cfvv.animation(edof,coord,dof,4,eig[:,eigenmode],10,mode_a*1000,def_scale=scalefact,export=True,file='anim/exv4b')
cfvv.animation(edof,coord,dof,4,eig[:,eigenmode],mode_a*1000,def_scale=scalefact,negative=True)
#cfvv.animate(edof,coord,dof,4,eig[:,0],10,def_scale=scalefact,export=True)

# For exporting animation
#cfvv.animation(edof,coord,dof,4,eig[:,eigenmode],mode_a*1000,def_scale=scalefact,negative=True,scalar_title='Tot. el. displacement [mm]',export=True,file='export/exv4/anim/exv4_modal')


cfvv.add_scalar_bar('Tot. el. displacement [mm]')

#Start Calfem-vedo visualization
cfvv.show_and_wait()

#For exporting meshes
#cfvv.export_vtk('export/exv4/exv4_static', static)
#cfvv.export_vtk('export/exv4/exv4_modal', modal)







'''
        #print(a)

        #ncoord = np.size(coord, axis = 0)
        #nnode = np.size(coord, axis = 0)

        ex,ey,ez = cfc.coordxtr(edof,coord,dof)

        ed = cfc.extractEldisp(edof,a)
        if element_type == 3:
            if val != 'nodal_values_by_el':
                coord2, topo, node_dofs, a_node, node_scalars = vdu.convert_to_node_topo(edof,ex,ey,ez,ed,ignore_first=False,dofs_per_node=1)
            else:
                coord2, topo, node_dofs, a_node, node_scalars = vdu.convert_to_node_topo(edof,ex,ey,ez,ed,values,ignore_first=False,dofs_per_node=1)
        elif element_type == 4:
            if val != 'nodal_values_by_el':
                coord2, topo, node_dofs, a_node, node_scalars = vdu.convert_to_node_topo(edof,ex,ey,ez,ed,ignore_first=False)
            else:
                coord2, topo, node_dofs, a_node, node_scalars = vdu.convert_to_node_topo(edof,ex,ey,ez,ed,values,ignore_first=False)
        #a_node = vdu.convert_a(coord,coord2,a,3)

        #def_coord = np.zeros([nnode,3])

        def_coord = coord2 + a_node*def_scale

        #print(a_node)
        
        #for i in range(np.size(def_coord, axis = 0)):
            #def_coord[i] = coord2[i] + a_node[i]
        #    def_coord[i,0] = coord2[i,0] + a[i*3]
        #    def_coord[i,1] = coord2[i,1] + a[i*3+1]
        #    def_coord[i,2] = coord2[i,2] + a[i*3+2]

        


        #print(def_coord)
        
        """
        for i in range(nnode):
        #a_dx, a_dy, a_dz = get_a_from_coord(i,3,a,def_scale)
        #x = coord[i][0]+a_dx
        #y = coord[i][1]+a_dy
        #z = coord[i][2]+a_dz
        #def_coord[i] = [x,y,z]
            def_coord[i,0] = a[i*3]
            def_coord[i,1] = a[i*3+1]
            def_coord[i,2] = a[i*3+2]
        """

        #print(a)
        #print(np.size(a, axis = 0))
        #print(np.size(a, axis = 1))
        """
        for i in range(0, ncoord):
            #if a.any() == None:
            #    x = coord[i,0]
            #    y = coord[i,1]
            #    z = coord[i,2]
            #else:
            #a_dx, a_dy, a_dz = get_a_from_coord(i,3,a,def_scale)

            #x = coord[i,0]+a_dx
            #y = coord[i,1]+a_dy
            #z = coord[i,2]+a_dz

            x = coord[i][0]+a[i][0]*scale
            y = coord[i][1]+a[i][1]*scale
            z = coord[i][2]+a[i][2]*scale

            def_coord[i] = [x,y,z]

            #def_nodes.append(v.Sphere(c='white').scale(1.5*scale).pos([x,y,z]).alpha(alpha))
        """
        #meshes = []
        #nel = np.size(edof, axis = 0)
        
        

        #print(topo)

        #mesh = v.Mesh([def_coord, topo]).lw(1)
        




        ct = vtk.VTK_HEXAHEDRON

        celltypes = [ct] * nel

        ug=v.UGrid([def_coord, topo, celltypes])
        ug.points(def_coord)
        
        mesh = ug.tomesh().lw(1).alpha(alpha)

        #v.settings.useDepthPeeling = True

        #print(val)

        #print('Cell connectivity: ',mesh.faces())

        #elif val and val == 'nodal_values':
        if val and val == 'el_values':
            #print(val)
            #vmin, vmax = np.min(values), np.max(values)
            
            el_values = vdu.convert_el_values(edof,values)
            mesh.celldata["val"] = el_values

            mesh.cmap(colormap, "val", on="cells")
'''




