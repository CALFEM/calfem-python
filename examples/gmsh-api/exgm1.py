import gmsh
import sys

lc = 0

gmsh.initialize()
gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1, 1)
gmsh.model.occ.synchronize()

#gmsh.model.add_physical_group(1, [1, 2, 4], 5)
#ps = gmsh.model.add_physical_group(2, [1])
#gmsh.model.set_physical_name(2, ps, "My surface")

#gmsh.option.setNumber("Mesh.MeshSizeFactor", 0.2)
#gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay for 2D meshes
#gmsh.option.setNumber("Mesh.Algorithm", 8) # Frontal-Delaunay for quads
#gmsh.model.mesh.setAlgorithm(2, 33, 1)
#gmsh.option.setNumber("Mesh.Smoothing", 100)
#gmsh.option.setNumber("Mesh.RecombineAll", 1)
#gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2) # or 3

gmsh.model.mesh.generate(3)

entities = gmsh.model.get_entities()

tag = 1
dim = 3

node_tags, node_coords, node_params = gmsh.model.mesh.get_nodes(dim, tag)
elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.get_elements(dim, tag)
type = gmsh.model.get_type(dim, tag)
name = gmsh.model.getEntityName(dim, tag)

print("type =", type)
print("tag  =", tag)

print("node_tags:", len(node_tags))
print(node_tags)
print("node_coords:", int(node_coords.shape[0]/3))
print(node_coords.reshape((int(node_coords.shape[0]/3), 3)))
print("node_params:", len(node_params))
print(node_params)
print("elem_types:", len(elem_types))
print(elem_types)
print("elem_tags:", len(elem_tags))
print(elem_tags[0])
print("elem_node_tags:", len(elem_node_tags))
print(elem_node_tags[0])

# for e in entities:
#     dim = e[0]
#     tag = e[1]
#     print(dim, tag)
#     node_tags, node_coords, node_params = gmsh.model.mesh.get_nodes(dim, tag)
#     elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.get_elements(dim, tag)

#     type = gmsh.model.getType(e[0], e[1])
#     name = gmsh.model.getEntityName(e[0], e[1])
#     #if len(name): name += ' '
#     #print("Entity " + name + str(e) + " of type " + type)

#     num_elem = sum(len(i) for i in elem_tags)
#     #print(" - Mesh has " + str(len(node_tags)) + " nodes and " + str(num_elem) +
#     #      " elements")

#     for i in elem_node_tags:
#         print(i)

    

    # # * Upward and downward adjacencies:
    # up, down = gmsh.model.get_adjacencies(e[0], e[1])
    # if len(up):
    #     print(" - Upward adjacencies: " + str(up))
    # if len(down):
    #     print(" - Downward adjacencies: " + str(down))

    # physicalTags = gmsh.model.getPhysicalGroupsForEntity(dim, tag)
    # if len(physicalTags):
    #     s = ''
    #     for p in physicalTags:
    #         n = gmsh.model.getPhysicalName(dim, p)
    #         if n: n += ' '
    #         s += n + '(' + str(dim) + ', ' + str(p) + ') '
    #     print(" - Physical groups: " + s)

    # partitions = gmsh.model.getPartitions(e[0], e[1])
    # if len(partitions):
    #     print(" - Partition tags: " + str(partitions) + " - parent entity " +
    #           str(gmsh.model.getParent(e[0], e[1])))

    # # * List all types of elements making up the mesh of the entity:
    # for t in elem_types:
    #     name, dim, order, numv, parv, _ = gmsh.model.mesh.getElementProperties(t)
    #     print(" - Element type: " + name + ", order " + str(order) + " (" +
    #           str(numv) + " nodes in param coord: " + str(parv) + ")")

gmsh.write("exgmsh1.msh")

#if '-nopopup' not in sys.argv:
#gmsh.fltk.run()

gmsh.finalize()