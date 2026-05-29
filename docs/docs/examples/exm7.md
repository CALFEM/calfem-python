<div class="cell markdown">

# Mesh example: Temperature distribution in a plate

Meshing 8-node-isoparametric elements (second order incomplete quads).
Shows use of surfacemarkers/elementmarkers to apply different properties
to elements in different regions.

</div>

<div class="cell code" execution_count="1">

``` python
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_mpl as cfv
import calfem.utils as cfu
import calfem.core as cfc

import numpy as np
```

</div>

<div class="cell code" execution_count="2">

``` python
# %matplotlib notebook
%matplotlib inline
```

</div>

<div class="cell markdown">

## Problem variables

</div>

<div class="cell code" execution_count="3">

``` python
kx1 = 100
ky1 = 100
kx2 = 10
ky2 = 10
t = 1.0

# Gauss points or integration points 

n = 2 
ep = [t, n]

D1 = np.matrix([
    [kx1, 0.],
    [0., ky1]
])
D2 = np.matrix([
    [kx2, 0.],
    [0., ky2]
])

# markers 10 & 11 will be used to specify different regions with different 
# conductivity.

Ddict = {10 : D1, 11 : D2} 
```

</div>

<div class="cell markdown">

## Define geometry

</div>

<div class="cell code" execution_count="4">

``` python
g = cfg.geometry()
```

</div>

<div class="cell markdown">

### Add points

</div>

<div class="cell code" execution_count="5">

``` python
points = [
    [0,0], 
    [0,100], 
    [0,150], 
    [100,0], 
    [150,0], 
    [100,-100], 
    [150,-100]
]

for p in points:
    g.point(p)
```

</div>

<div class="cell markdown">

### Add splines

</div>

<div class="cell code" execution_count="6">

``` python
g.spline([1,2], marker=2, el_on_curve=4)
g.spline([3,4], el_on_curve=4)
g.circle([1,0,3], el_on_curve = 10)
g.circle([2,0,4], el_on_curve = 10)
g.spline([3,5], el_on_curve = 6)
g.spline([5,6], marker=3, el_on_curve = 4)
g.spline([6,4], el_on_curve = 6)
```

</div>

<div class="cell markdown">

### Add surfaces

</div>

<div class="cell code" execution_count="7">

``` python
g.structured_surface([0,2,1,3], marker = 10)
g.structured_surface([1,4,5,6], marker = 11)
```

</div>

<div class="cell markdown">

## Generate mesh

</div>

<div class="cell code" execution_count="8">

``` python
el_type = 16 
dofs_per_node = 1 

mesh = cfm.GmshMesh(g, el_type, dofs_per_node) 

coords, edof, dofs, bdofs, elementmarkers = mesh.create()
```

<div class="output stream stdout">

    Info    : GMSH -> Python-module

</div>

</div>

<div class="cell markdown">

## Solve problem

</div>

<div class="cell markdown">

### Assemble system matrix

</div>

<div class="cell code" execution_count="9">

``` python
n_dofs = np.size(dofs)
ex, ey = cfc.coordxtr(edof, coords, dofs)

K = np.zeros([n_dofs,n_dofs])

for eltopo, elx, ely, elMarker in zip(edof, ex, ey, elementmarkers):

    # Calc element stiffness matrix: Conductivity matrix D is taken 
    # from Ddict and depends on which region (which marker) the element is in.

    Ke = cfc.flw2i8e(elx, ely, ep, Ddict[elMarker]) 
    cfc.assem(eltopo, K, Ke)
```

</div>

<div class="cell markdown">

### Solving equation system

</div>

<div class="cell code" execution_count="10">

``` python
f = np.zeros([n_dofs,1])

bc = np.array([],'i')
bc_val = np.array([],'i')

bc, bc_val = cfu.applybc(bdofs,bc,bc_val,2,30.0)
bc, bc_val = cfu.applybc(bdofs,bc,bc_val,3,0.0)

a,r = cfc.solveq(K,f,bc,bc_val)
```

</div>

<div class="cell markdown">

### Compute element forces

</div>

<div class="cell code" execution_count="11">

``` python
ed = cfc.extract_eldisp(edof,a)

for i in range(np.shape(ex)[0]):
    es, et, eci = cfc.flw2i8s(ex[i,:], ey[i,:], ep, Ddict[elementmarkers[i]], ed[i,:])
```

</div>

<div class="cell markdown">

## Visualise results

</div>

<div class="cell code" execution_count="12">

``` python
cfv.figure(fig_size=(10, 10))
cfv.draw_geometry(g, title="Geometry")
```

<div class="output display_data">

![](62d10b408e9333b3a1aa8ddd3bf6140f52180684.png)

</div>

</div>

<div class="cell code" execution_count="13">

``` python
cfv.figure(fig_size=(10, 10))
cfv.draw_mesh(coords, edof, dofs_per_node, el_type, filled=False)
```

<div class="output display_data">

![](875bf425be6805c7d0f32ece66959c4291008716.png)

</div>

</div>

<div class="cell code" execution_count="17">

``` python
cfv.figure(fig_size=(10, 10))
cfv.draw_nodal_values_shaded(a, coords, edof, title="Temperature", dofs_per_node=mesh.dofs_per_node, el_type=mesh.el_type, draw_elements=True)
cbar = cfv.colorbar(orientation="vertical")
cbar.set_label("Temperature")

cfv.text("The bend has high conductivity", (77,135))
cfv.text("This part has low conductivity", (20,-90))
```

<div class="output execute_result" execution_count="17">

    Text(20, -90, 'This part has low conductivity')

</div>

<div class="output display_data">

![](a49b9cc8a5186172d392f3e2a1138344a7ae9900.png)

</div>

</div>

<div class="cell code">

``` python
```

</div>
