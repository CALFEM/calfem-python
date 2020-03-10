# -*- coding: utf-8 -*-

'''Example 04.

Structured 3D meshing. Adding texts and labels to figures. Altering axis properties.
''' 

import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_mpl as cfv

# ---- Define geometry ------------------------------------------------------

g = cfg.geometry()

# Add Points

g.point([0, 0, 0], ID=0)
g.point([0.5, -0.3, -0.3], 1)
g.point([1, 0, 0], 2)
g.point([1, 1, 0], 3)
g.point([0, 1, 0], 4, marker = 11) # Set some markers no reason.
g.point([0, 0, 1], 5, marker = 11) # (markers can be given to points as well
                                      # as curves and surfaces)
g.point([1, 0, 1], 6, marker = 11)
g.point([1, 1, 1], 7)
g.point([0, 1, 1], 8)

# Add splines

g.spline([0, 1, 2], 0, marker = 33, el_on_curve = 5)
g.spline([2, 3], 1, marker = 23, el_on_curve = 5)
g.spline([3, 4], 2, marker = 23, el_on_curve = 5)
g.spline([4, 0], 3, el_on_curve = 5)
g.spline([0, 5], 4, el_on_curve = 5)
g.spline([2, 6], 5, el_on_curve = 5)
g.spline([3, 7], 6, el_on_curve = 5)
g.spline([4, 8], 7, el_on_curve = 5)
g.spline([5, 6], 8, el_on_curve = 5)
g.spline([6, 7], 9, el_on_curve = 5)
g.spline([7, 8], 10, el_on_curve = 5)
g.spline([8, 5], 11, el_on_curve = 5)

# Add surfaces

g.structuredSurface([0, 1, 2, 3], 0, marker=45)
g.structuredSurface([8, 9, 10, 11], 1)
g.structuredSurface([0, 4, 8, 5], 2, marker=55)
g.structuredSurface([1, 5, 9, 6], 3, marker=55)
g.structuredSurface([2, 6, 10, 7], 4)
g.structuredSurface([3, 4, 11, 7], 5)

# Add Volume:
#  addStructuredVolume() takes three args. The first is a list of surface IDs 
#  (structured surfaces). The surfaces should make a hexahedron 
#  (i.e. 6 surfaces). Other kinds of structured volumes than hexahedra will
#  not work for hexahedral elements, which is the only type of 3D element that 
#  CALFEM handles. The two optional parameters are the volume ID and 
#  volume marker.

g.structuredVolume([0,1,2,3,4,5], 0, marker=90)

# ---- Create mesh ----------------------------------------------------------

# Element type 5 is hexahedron. (See user manual for more element types)

el_type = 5 

# Degrees of freedom per node.

dofs_per_node = 1 

# Create mesh

coords, edof, dofs, bdofs, _ = cfm.mesh(g, el_type, dofs_per_node)

#coords, edof, dofs, bdofs, _ = cfm.mesh(
#        g, el_type, dofs_per_node, gmsh_exec_path="D:\\vsmn20-software\\gmsh\gmsh.exe")

# ---- Visualise mesh -------------------------------------------------------

# Hold Left Mouse button to rotate.
# Hold right mouse button to zoom.
# Hold SHIFT and left mouse button to pan.
# Hold SHIFT and right mouse button to change the field of view.
# Hold Ctrl and left mouse button to roll the camera. 

# Draw geometry

cfv.draw_geometry(g, draw_points=False)

# Draw mesh

cfv.figure()
cfv.draw_mesh(coords=coords, edof=edof, dofs_per_node=dofs_per_node, el_type=el_type, filled=True)

# Add a text in world space

#cfv.add_text("This is a Text", pos=(1, 0.5, 0.5), angle=45)  

# Add a label in the screen space

#our_label = cfv.add_label("This is a Label", pos=(20,30), angle=-45) 

# We can change the attributes of labels and texts, such as color and position.

#our_label.text = "Label, changed." 

# Make it red. (1,0,0) would also have worked.

#our_label.textColor = 'r'  

# Matlab style axes (three axes in the background instead of a cube)

#cfv.gca().axis.showBox = 0 

#Change the limits of the axes.

#cfv.gca().SetLimits(rangeX=(0,2), rangeY=(-1,1.5), rangeZ=(-0.5,2), margin=0.02) 

# Enter main loop

cfv.show_and_wait()