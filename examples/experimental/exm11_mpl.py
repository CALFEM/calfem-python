# -*- coding: utf-8 -*-

'''Example 10

The use case from the user manual. 
The example does not contain anything that is not covered in the previous examples.
'''

import calfem.core as cfc
import calfem.vis_mpl as cfv
import calfem.utils as cfu
import calfem.shapes as cfs
import calfem.solver as cfslv

cfu.enableLogging()
  
# ---- General parameters ---------------------------------------------------

# Define marker constants instead of using numbers in the code

cfu.info("Creating rectangle")

rect = cfs.Rectangle(5.0, 1.0, element_type=3, dofs_per_node=2, max_area=0.08)
rect.t = 0.2
rect.v = 0.35
rect.E = 2e9
rect.ptype = 1
rect.ep = [rect.ptype, rect.t]
rect.D = cfc.hooke(rect.ptype, rect.E, rect.v)

cfu.info("Creating mesh...")

mesh = cfs.ShapeMesh(rect)    

# ---- Solve problem --------------------------------------------------------

solver = cfslv.Plan2DSolver(mesh)

solver.addBC(rect.left_id, 0.0)
solver.addForceTotal(rect.top_id, -10e5, dimension=2)

results = solver.execute()       
       
# ---- Visualise results ----------------------------------------------------

cfu.info("Drawing results...")

cfv.figure() 
cfv.draw_geometry(rect.geometry(), title="Geometry")

cfv.figure() 
cfv.draw_mesh(mesh.coords, mesh.edof, rect.dofs_per_node, rect.element_type, 
             filled=True, title="Mesh") #Draws the mesh.

cfv.figure()
cfv.draw_displacements(results.a, mesh.coords, mesh.edof, rect.dofs_per_node, rect.element_type, 
                      draw_undisplaced_mesh=False, title="Displacements", 
                      magnfac=1)

cfv.figure()
cfv.draw_element_values(results.el_forces, mesh.coords, mesh.edof, rect.dofs_per_node, rect.element_type, results.a, 
                      draw_elements=True, draw_undisplaced_mesh=False, 
                      title="Effective Stress", magnfac=1)
                      
#cfv.colorBar().SetLabel("Effective stress")

cfu.info("Done drawing...")

cfv.show_and_wait()