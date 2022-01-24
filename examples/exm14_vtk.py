# -*- coding: utf-8 -*-

'''Example 05

This example shows how to make an unstructured 3D mesh (tetrahedron elements, which calfem cant actually use).
It also demonstrates how to do subplots and create two axes that are viewed from the same camera.
''' 

import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_vtk as cfv

import prim3d as p3

import vtk

# ---- Define geometry ------------------------------------------------------

#box1 = p3.Box(0.0, 0.0, 0.0)
#box2 = p3.Box(1.0, 0.0, 0.0)

#model = p3.GeometryModel()
#model.add_volume(box1.volume)
#model.add_volume(box2.volume)

g = cfg.Geometry()

# ---- Create mesh ----------------------------------------------------------

# Element type 4 is tetrahedron. (See user manual for more element types).

el_type = 5

# Degrees of freedom per node.

dofs_per_node = 1 

# Create mesh

mesh_generator = cfm.GmshMesh("test.geo", el_type, 0.05, dofs_per_node)
mesh_generator.mesh_dir = "mesh"

mesh_generator.gmsh_options = {"Mesh.Algorithm":6, "Mesh.Algorithm3D":1, "Mesh.RecombinationAlgorithm":1,
"Mesh.RecombineAll":1, "Mesh.Tetrahedra":1}

coords, edof, dofs, bdofs, elementmarkers = mesh_generator.create(is3D=True)

#coords, edof, dofs, bdofs, elementmarkers = cfm.mesh(model.geometry, el_type, 0.3, dofs_per_node)

print(coords)
print(edof)
#print(bdofs)
#print(elementmarkers)

#print(bdofs[100])
#print(edof[0])

#coords, edof, dofs, bdofs, _ = cfm.mesh(
#        g, el_type, 0.3, dofs_per_node, gmsh_exec_path="D:\\vsmn20-software\\gmsh\gmsh.exe")


# ---- Visualise mesh -------------------------------------------------------

cfv.draw_mesh(coords, edof, el_type)