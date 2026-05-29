<div class="cell markdown">

# Mesh example: More geometry

This example is from from Calfem for Pyhon Mesh Manual (Mesh_Ex_02.py)

Creating geometry from B-Splines and circle arcs. Also shows how to set
ID numbers for geometry entities and how to specify element density.

</div>

<div class="cell code" execution_count="2">

``` python
import numpy as np

import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_mpl as cfv
```

</div>

<div class="cell markdown">

**Define geometry**

All mesh generation starts by creating a Geometry object.

</div>

<div class="cell code" execution_count="3">

``` python
g = cfg.Geometry()
```

</div>

<div class="cell markdown">

**Add points**

In this example we set the IDs manually.

</div>

<div class="cell code" execution_count="4">

``` python
g.point([-2,  0], ID=0)
```

</div>

<div class="cell markdown">

elSize determines the size of the elements near this point.

</div>

<div class="cell code" execution_count="5">

``` python
g.point([0,  1], ID=1, el_size=5)
```

</div>

<div class="cell markdown">

elSize is 1 by default. Larger number means less dense mesh. Size means
the length of the sides of the elements.

</div>

<div class="cell code" execution_count="6">

``` python
g.point([1,  0], 2, el_size=5)
g.point([0, -2], 3)  
g.point([0,  0], 4, el_size=5)
g.point([.5, .2], 5)
g.point([-.5, .5], 6)
g.point([-.7, -.5], 7)
```

</div>

<div class="cell markdown">

**Add curves**

The 3 points that define the circle arc are \[start, center, end\]. The
arc must be smaller than Pi.

</div>

<div class="cell code" execution_count="7">

``` python
g.circle([1, 4, 2], 2)
```

</div>

<div class="cell markdown">

BSplines are similar to Splines, but do not necessarily pass through the
control points.

</div>

<div class="cell code" execution_count="8">

``` python
g.bspline([5, 6, 7, 5], 5)
g.bspline([1, 0, 3, 2], 4)
```

</div>

<div class="cell markdown">

**Add surfaces**

</div>

<div class="cell code" execution_count="9">

``` python
g.surface([4, 2], [[5]])
```

</div>

<div class="cell markdown">

Markers do not have to be set when the curve is created. It can be done
afterwards. Set marker=80 for curves 2 and 4:

</div>

<div class="cell code" execution_count="10">

``` python
for curveID in [2, 4]:
    g.curve_marker(curveID, 80)
```

</div>

<div class="cell markdown">

**Generate mesh**

Create a mesh object and set its attributes.

</div>

<div class="cell code" execution_count="11">

``` python
mesh = cfm.GmshMesh(g)

mesh.el_type = 3
mesh.dofs_per_node = 2
mesh.el_size_factor = 0.05
```

</div>

<div class="cell markdown">

Create the finite element mesh using the create() method of the mesh
object.

</div>

<div class="cell code" execution_count="12">

``` python
coords, edof, dofs, bdofs, elementmarkers = mesh.create()
```

<div class="output stream stdout">

    Info    : GMSH -> Python-module

</div>

</div>

<div class="cell markdown">

**Visualise mesh**

***Draw geometry***

</div>

<div class="cell code" execution_count="13">

``` python
%matplotlib inline
# %matplotlib widget
```

</div>

<div class="cell code" execution_count="14">

``` python
cfv.figure(fig_size=(8,8))
cfv.draw_geometry(
    g,
    label_curves=True,
    title="Example 2 - Geometry"
)
```

<div class="output display_data">

![](c682f8fcc4226d9c22a34745f1ef814b45cfd74f.png)

</div>

</div>

<div class="cell markdown">

**Draw mesh**

</div>

<div class="cell code" execution_count="15">

``` python
cfv.figure(fig_size=(8,8))
cfv.draw_mesh(
    coords=coords,
    edof=edof,
    dofs_per_node=mesh.dofs_per_node,
    el_type=mesh.el_type,
    filled=True,
    title="Example 2 - Mesh"
)
```

<div class="output display_data">

![](7bb56c15b5125072c76300c64bff8b5cbe4a52cb.png)

</div>

</div>
