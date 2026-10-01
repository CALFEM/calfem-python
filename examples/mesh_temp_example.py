import sys
#!pip install calfem_python
#!apt-get install libglu1-mesa -y

# -*- coding: utf-8 -*-
"Läs in programbibliotek som används"
import numpy as np
import calfem.core as cfc
import calfem.vis_mpl as cfv
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.utils as cfu
import matplotlib.pyplot as plt


g = cfg.Geometry()
g.point([0.0, 0.0]) # point 0
g.point([12.0, 0.0]) # point 1
g.point([17.0, 0.0]) # point 2
g.point([17.0, 7.0]) # point 3
g.point([12.0, 7.0]) # point 4
g.point([12.0, 10.0]) # point 6
g.point([0.0, 10.0]) # point 7
g.point([0.0, 7.0]) # point 8

elmult=2
nrofel1=9*elmult
nrofel2=4*elmult
nrofel3=7*elmult
nrofel4=3*elmult

g.spline([0, 1],el_on_curve=nrofel1) # line 0
g.spline([1, 2],el_on_curve=nrofel2) # line 1
g.spline([2, 3],marker=20,el_on_curve=nrofel3) # line 2
g.spline([3, 4],el_on_curve=nrofel2) # line 3
g.spline([4, 5],el_on_curve=nrofel4) # line 4
g.spline([5, 6],marker=10,el_on_curve=nrofel1) # line 5
g.spline([6, 7],marker=10,el_on_curve=nrofel4) # line 6
g.spline([7, 0],marker=10,el_on_curve=nrofel3) # line 7
g.spline([1, 4],el_on_curve=nrofel3) # line 8
g.spline([7, 4],el_on_curve=nrofel1) # line 9

g.struct_surface([0, 8, 9, 7])
g.struct_surface([1, 2, 3, 8])
g.struct_surface([9, 4, 5, 6])

cfv.draw_geometry(g)
cfv.showAndWait()


mesh = cfm.GmshMesh(g)
mesh.el_type = 3          # Element type is quadrilateral
mesh.dofs_per_node = 1     # Degrees of freedom per node

coords, edof, dofs, bdofs, elementmarkers = mesh.create()

# Draw the mesh.
cfv.figure()
cfv.drawMesh(coords=coords,edof=edof,dofs_per_node=mesh.dofsPerNode,el_type=mesh.elType)

n_dofs = np.size(dofs)
n_el = np.size(edof,0)

ex, ey = cfc.coordxtr(edof, coords, dofs)
K = np.zeros([n_dofs,n_dofs])
ep=[1]
D=np.array([[0.1, 0.0],[0.0, 0.1],])

for eltopo, elx, ely in zip(edof, ex, ey):
    Ke = cfc.flw2qe(elx, ely, ep, D)
    cfc.assem(eltopo, K, Ke)

f = np.zeros([n_dofs,1])

bc = np.array([],'i')
bcVal = np.array([],'f')

bc, bcVal = cfu.applybc(bdofs, bc, bcVal, 20, 20.0, 0)
bc, bcVal = cfu.applybc(bdofs, bc, bcVal, 10, 0.0, 0)

a,r = cfc.solveq(K,f,bc,bcVal)

ed = cfc.extract_ed(edof, a)

es = np.zeros([n_el,2])
et = np.zeros([n_el,2])
for i in range(0, n_el): #elx, ely, eld in zip(ex, ey, ed):
  [es[i,:], et[i,:]] = cfc.flw2qs(ex[i,:], ey[i,:], ep, D, ed[i,:])
  #[es, et] = cfc.flw2qs(elx, ely, ep, D, eld)


cfv.figure()
cfv.drawMesh(coords=coords,edof=edof,dofs_per_node=mesh.dofsPerNode,el_type=mesh.elType)
cfv.draw_nodal_values_contour(a, coords, edof,levels=20)

cfv.figure()
cfv.drawMesh(coords=coords,edof=edof,dofs_per_node=mesh.dofsPerNode,el_type=mesh.elType)
cfv.draw_nodal_values_contourf(a, coords, edof,levels=20)

cfv.figure()
cfv.drawMesh(coords=coords,edof=edof,dofs_per_node=mesh.dofsPerNode,el_type=mesh.elType)
cfv.draw_nodal_values_shaded(a, coords, edof)


fig, ax = plt.subplots()
#cfv.figure()
cfv.drawMesh(coords=coords,edof=edof,dofs_per_node=mesh.dofsPerNode,el_type=mesh.elType)

sfac = cfv.elflux2(
    ex,
    ey,
    es,
    plotcolor=[4],   #  red
    ax=ax
)

ax.set_aspect("equal")
plt.show()
