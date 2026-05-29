<div class="cell markdown">

# Mesh example: Meshing in 2D

Shows structured meshing in 2D.

</div>

<div class="cell code" execution_count="1">

``` python
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_mpl as cfv
```

</div>

<div class="cell markdown">

**Define geometry**

</div>

<div class="cell code" execution_count="4">

``` python
g = cfg.Geometry()
```

</div>

<div class="cell markdown">

**Define points**

</div>

<div class="cell code" execution_count="5">

``` python
g.point([0,0])
g.point([1.2, 0])
g.point([1, 1.3])
g.point([0, 1])
g.point([2, 0.5])
```

</div>

<div class="cell markdown">

**Add splines**

The first four curves are structured curves, i.e the number of nodes
along the curves is pre-determined. Parameter elOnCurve states how many
elements are placed along the curve. Parameters elDistribType and
elDistribVal are optional parameters that specify how elements are
distributed.

- "bump" means elements are bunched up at the ends or the middle of the
  curve. In this case elDistribVal is smaller than 1, so elements crowd
  at the edges.
- "progression" means each element along the curve is larger/smaller
  than the previous one. A larger elDistribVal makes the elements larger
  at the end of the curves.

</div>

<div class="cell code" execution_count="6">

``` python
g.spline([0,1], el_on_curve=10, el_distrib_type="bump", el_distrib_val=0.2)
g.spline([1,2], el_on_curve=20, el_distrib_type="progression", el_distrib_val=1.1)
g.spline([2,3], el_on_curve=10, el_distrib_type="bump", el_distrib_val=0.2)
g.spline([0,3], el_on_curve=20, el_distrib_type="progression", el_distrib_val=1.1) #Change order of points to reverse progression distribution
g.spline([2, 4, 1])
```

</div>

<div class="cell markdown">

**Add surfaces**

A structured surface must contain 4 curves that have the parameter
'elOnCurve' defined. The number of elements on two opposite curves must
be the same (In this case, curves 0 & 2 and 1 & 3).

</div>

<div class="cell code" execution_count="8">

``` python
g.structured_surface([0,1,2,3]) 
g.surface([4,1])
```

</div>

<div class="cell markdown">

**Create mesh**

</div>

<div class="cell code" execution_count="9">

``` python
mesh = cfm.GmshMesh(g)

mesh.el_type = 3       # quad
mesh.dofs_per_node = 1 
mesh.el_size_factor = 0.01

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

<div class="cell code" execution_count="10">

``` python
%matplotlib inline
```

</div>

<div class="cell code" execution_count="11">

``` python
cfv.figure(fig_size=(10,10))
cfv.draw_geometry(g)
```

<div class="output display_data">

![](5a73c58fdda5f90f8d912441af4b027df34f95b0.png)

</div>

</div>

<div class="cell markdown">

**Draw mesh**

</div>

<div class="cell code" execution_count="12">

``` python
cfv.figure(fig_size=(10,10))
cfv.draw_mesh(coords=coords, edof=edof, dofs_per_node=mesh.dofs_per_node, el_type=mesh.el_type, filled=True)
```

<div class="output display_data">

![](bf6034c69941fc26c6ee79e92b0bb30e58c7e19b.png)

</div>

</div>
