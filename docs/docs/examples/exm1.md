<div class="cell markdown">

# Mesh example: Simple mesh geometry

</div>

<div class="cell markdown">

This example is from Calfem for Python Mesh Manual (Mesh_Ex_01.py)

Shows how to create simple geometry from splines and ellipse arcs, and
how to mesh a quad mesh in GmshMesher. Also demonstrates drawGeometry(),
drawMesh, and drawing texts and labels in a figure.

</div>

<div class="cell code" execution_count="1">

``` python
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_mpl as cfv
```

</div>

<div class="cell markdown">

## Define geometry

</div>

<div class="cell markdown">

Create geometry object.

</div>

<div class="cell code" execution_count="2">

``` python
g = cfg.Geometry()
```

</div>

<div class="cell markdown">

### Add points to geometry object.

The first parameter is the coordinates. These can be in 2D or 3D. The
other parameters are not defined in this example. These parameters are
ID, marker, and elSize. Since we do not specify an ID the points are
automatically assigned IDs, starting from 0.

</div>

<div class="cell code" execution_count="3">

``` python
g.point([0, 0])
g.point([2, 0])
g.point([2, 1])
g.point([0, 1])
g.point([0.5, 0.3])
g.point([0.3, 0.7])
g.point([0.7, 0.7])
g.point([0.8, 0.5])
g.point([1.7, 0.5])
g.point([1.5, 0.5])
g.point([1.7, 0.7])
```

</div>

<div class="cell markdown">

### Add curves

</div>

<div class="cell markdown">

There are four types of curves. In this example we create an ellipse arc
and some splines. The first parameter is a list of point IDs that define
the curve. Curves can have have IDs and markers. In this example the IDs
are undefined so the curves are automatically assigned IDs. The markers
can be used for identifying regions/boundaries in the model.

</div>

<div class="cell code" execution_count="4">

``` python
g.ellipse([7, 8, 9, 10], marker=50)

g.spline([0, 1], marker=80)           # 1 - A spline. Splines pass through the
                                      #     points in the first parameter.
g.spline([2, 1])                      # 2
g.spline([3, 2])                      # 3
g.spline([0, 3])                      # 4
g.spline([7, 9], marker=50)           # 5
g.spline([10, 9])                     # 6
g.spline([4, 5, 6, 4])                # 7 - This is a closed spline.
                                      #     The start and end points are the same
```

</div>

<div class="cell markdown">

### Add a surface

Surfaces are defined by its curve boundaries. The first parameter is a
list of curve IDs that specify the outer boundary of the surface. The
second parameter is a list of lists of curve IDs that specify holes in
the surface. In this example there are two holes. The boundaries and
holes must be closed paths. We can see that \[7\] is closed because
curve 7 is a closed spline. addSurface creates a flat surface, so all
curves must lie on the same plane.

</div>

<div class="cell code" execution_count="5">

``` python
g.surface([4, 3, 2, 1], [[7], [5, 6, 0]])
```

</div>

<div class="cell markdown">

## Generate mesh

### Define mesh properties

</div>

<div class="cell code" execution_count="6">

``` python
mesh = cfm.GmshMesh(g)
```

</div>

<div class="cell markdown">

Element type 3 is quad. 2 is triangle. See user manual for more element
types)

</div>

<div class="cell code" execution_count="7">

``` python
mesh.el_type = 3
mesh.dofs_per_node = 1  # Degrees of freedom per node.
mesh.el_size_factor = 0.05  # Factor that changes element sizes.
```

</div>

<div class="cell markdown">

### Create mesh

The first four return values are the same as those that trimesh2d()
returns. coords is as list of node coordinates. edof is the element
topology (element degrees of freedom). dofs is a lists of all degrees of
freedom bdofs is a dictionary of boundary dofs (dofs of geometric
entities with markers). elementmarkers is a list of markers, and is used
for finding the marker of a given element (index).

</div>

<div class="cell code" execution_count="8">

``` python
coords, edof, dofs, bdofs, elementmarkers = mesh.create()
```

<div class="output stream stdout">

    Info    : GMSH -> Python-module

</div>

</div>

<div class="cell markdown">

## Visualise mesh

</div>

<div class="cell markdown">

### Draw geometry

</div>

<div class="cell code" execution_count="9">

``` python
#%matplotlib notebook
%matplotlib inline
# %matplotlib widget
```

</div>

<div class="cell markdown">

### Draw mesh

</div>

<div class="cell code" execution_count="10">

``` python
cfv.figure(fig_size=(10,10))
cfv.draw_mesh(
    coords=coords,
    edof=edof,
    dofs_per_node=mesh.dofs_per_node,
    el_type=mesh.el_type,
    filled=True,
    title="Example 01"
)
```

<div class="output display_data">

![](11bfb1fb0d4722027f27fd47d4f3c58c32cf0285.png)

</div>

</div>

<div class="cell code">

``` python
```

</div>
