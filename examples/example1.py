import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_mpl as cfv

g = cfg.Geometry() 

# Add points

g.point([0, 0]) # 0
g.point([1, 0]) # 1
g.point([1, 1]) # 2
g.point([0, 1]) # 3

# Add points for circle

r = 0.20

g.point([0.5, 0.5])   # 4
g.point([0.5, 0.5+r]) # 5
g.point([0.5, 0.5-r]) # 6
g.point([0.5+r, 0.5]) # 7
g.point([0.5-r, 0.5]) # 8

# Add lines

g.spline([0, 1]) # 0
g.spline([1, 2]) # 1
g.spline([2, 3]) # 2
g.spline([3, 0]) # 3

# TODO: Add circles HERE

g.circle([5, 4, 7]) # 4
g.circle([7, 4, 6]) # 5
g.circle([6, 4, 8]) # 6
g.circle([8, 4, 5]) # 7

# Add surface and hole UPDATE

g.surface([0, 1, 2, 3], [[7, 6, 5, 4]])

#cfv.figure(fig_size=(10.0,10.0))
#cfv.draw_geometry(g)

mesh = cfm.GmshMesh(g)

# Mesh properties

mesh.el_type = 3
mesh.dofs_per_node = 1  # Degrees of freedom per node.
mesh.el_size_factor = 0.040  # Factor that changes element sizes.

coords, edof, dofs, bdofs, elementmarkers = mesh.create()


cfv.figure(fig_size=(10,10))
cfv.draw_mesh(
    coords=coords,
    edof=edof,
    dofs_per_node=mesh.dofs_per_node,
    el_type=mesh.el_type,
    filled=True,
    title="Example 01"
)

cfv.show_and_wait_mpl()
