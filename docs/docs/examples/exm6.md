<div class="cell markdown">

# Mesh example: Plane stress problem

This example is from Calfem for Python Mesh Manual (Mesh_Ex_06.py)

Solves a plane stress 2D problem using a structured mesh. Shows how to
draw von Mises effective stress as an element value with
drawElementValues(). Shows use of GmshMesher attribute 'nodesOnCurve'
(dictionary that says which nodes are on a given geometry curve)

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

<div class="cell markdown">

**Define problem variables**

</div>

<div class="cell code" execution_count="19">

``` python
t = 0.2
v = 0.35
E = 2.1e9
ptype = 1
ep = [ptype,t]
D = cfc.hooke(ptype, E, v)
```

</div>

<div class="cell markdown">

**Define geometry**

</div>

<div class="cell code" execution_count="20">

``` python
g = cfg.geometry()
```

</div>

<div class="cell markdown">

**Create points**

Just a shorthand. We use this to make the circle arcs.

</div>

<div class="cell code" execution_count="21">

``` python
s2 = 1/np.sqrt(2) 
```

</div>

<div class="cell markdown">

Using this we then define our points:

</div>

<div class="cell code" execution_count="22">

``` python
points = [[0, 3], [2.5, 3], [3, 3], [4-s2, 3-s2], [4, 2],     #0-4
          [4+s2, 3-s2], [5, 3], [5.5, 3], [8,3], [0, 1.5],    #5-9
          [2.5, 1.5], [4, 1.5], [5.5, 1.5], [8, 1.5], [0, 0], #10-14
          [2.5, 0], [3, 0], [4-s2, s2], [4, 1], [4+s2, s2],   #15-19
          [5, 0], [5.5, 0], [8,0], [4,3], [4,0]]              #20-24
          
for xp, yp in points:
    g.point([xp*0.1, yp*0.1])
```

</div>

<div class="cell markdown">

Now we create our curves:

</div>

<div class="cell code" execution_count="23">

``` python
splines = [[0,1], [1,2], [6,7], [7,8], [8,13],          #0-4
           [13,22], [22,21], [21,20], [16,15], [15,14], #5-9
           [14,9], [9,0], [9,10], [10,1], [10, 15],     #10-14
           [10,11], [11,4], [11,18], [11,12], [12,7],   #15-19
           [12,21], [12,13], [3,10], [5,12], [10,17],   #20-24
           [12,19]]                                     #25
           
for s in splines:
    g.spline(s, el_on_curve=10)
```

</div>

<div class="cell markdown">

In this case we use special functions to assign IDs to our curves.

</div>

<div class="cell code" execution_count="24">

``` python
g.curve_marker(ID=4,  marker=7) #Assign marker 7 to the splines on the right.
g.curve_marker(ID=5,  marker=7) # We will apply a force on nodes with marker 7.
g.curve_marker(ID=10, marker=5) #Assign marker 5 to the splines on the left.
g.curve_marker(ID=11, marker=5) # The nodes with marker 5 will be locked in place.
```

</div>

<div class="cell markdown">

Next we create our circle arcs.

</div>

<div class="cell code" execution_count="25">

``` python
# Points in circle arcs are [start, center, end]

circle_arcs = [[2, 23, 3], [3, 23, 4], [4, 23, 5], [5, 23, 6],           #26-29
              [16, 24, 17], [17, 24, 18], [18, 24, 19], [19, 24, 20]]   #30-33
              
for c in circle_arcs:
    g.circle(c, el_on_curve=10)
```

</div>

<div class="cell markdown">

Finally we create our structured surfaces:

</div>

<div class="cell code" execution_count="26">

``` python
g.struct_surf([11,12,13,0]) #0
g.struct_surf([14, 12, 10, 9])
g.struct_surf([8, 30, 24, 14])
g.struct_surf([24, 31, 17, 15])
g.struct_surf([15, 16, 27, 22]) #4
g.struct_surf([22, 26, 1, 13]) 
g.struct_surf([16, 18, 23, 28])
g.struct_surf([19, 2, 29, 23])
g.struct_surf([19, 21, 4, 3]) #8
g.struct_surf([20, 6, 5, 21])
g.struct_surf([25, 20, 7, 33])
g.struct_surf([32, 17, 18, 25]) #11
```

</div>

<div class="cell markdown">

**Create mesh**

</div>

<div class="cell code" execution_count="27">

``` python
mesh = cfm.GmshMesh(g)

mesh.el_type = 3
mesh.dofs_per_node = 2

coords, edof, dofs, bdofs, elementmarkers = mesh.create()
```

<div class="output stream stdout">

    Info    : GMSH -> Python-module

</div>

</div>

<div class="cell markdown">

**Solve problem**

***Assemble system matrix***

</div>

<div class="cell code" execution_count="28">

``` python
n_dofs = np.size(dofs)
ex, ey = cfc.coordxtr(edof, coords, dofs)
K = np.zeros([n_dofs,n_dofs])

for eltopo, elx, ely in zip(edof, ex, ey):
    Ke = cfc.planqe(elx, ely, ep, D)
    cfc.assem(eltopo, K, Ke)
```

</div>

<div class="cell markdown">

**Solve equation system**

</div>

<div class="cell code" execution_count="29">

``` python
f = np.zeros([n_dofs,1])

bc = np.array([],'i')
bcVal = np.array([],'f')

bc, bcVal = cfu.applybc(bdofs, bc, bcVal, 5, 0.0, 0)

cfu.applyforce(bdofs, f, 7, 10e5, 1)

a,r = cfc.solveq(K,f,bc,bcVal)
```

</div>

<div class="cell markdown">

**Compute element forces**

</div>

<div class="cell code" execution_count="30">

``` python
ed = cfc.extract_eldisp(edof,a)

von_mises = []

# For each element:
for i in range(edof.shape[0]): 

    # Determine element stresses and strains in the element.
    es, et = cfc.planqs(ex[i,:], ey[i,:], ep, D, ed[i,:]) 

    # Calc and append effective stress to list.    
    von_mises.append(np.sqrt(np.power(es[0],2) - es[0]*es[1] + np.power(es[1],2) + 3*es[2] ) ) 

    ## es: [sigx sigy tauxy]
```

</div>

<div class="cell markdown">

**Visualise results**

</div>

<div class="cell code" execution_count="31">

``` python
%matplotlib inline
```

</div>

<div class="cell code" execution_count="32">

``` python
cfv.figure(fig_size=(10,10))
cfv.draw_geometry(g, draw_points=True, label_curves=True, label_points=True)
```

<div class="output display_data">

![](ad3352d8c4dfa1dadeecc3e36df1c7d11c60bef5.png)

</div>

</div>

<div class="cell code" execution_count="33">

``` python
cfv.figure(fig_size=(10,10))
cfv.draw_mesh(coords, edof, dofs_per_node=mesh.dofs_per_node, el_type=mesh.el_type)
```

<div class="output display_data">

![](8eb93e5d47d041dba14bfb4420fbc207e1da6e8c.png)

</div>

</div>

<div class="cell code" execution_count="34">

``` python
cfv.figure(fig_size=(10,10))
cfv.draw_element_values(von_mises, coords, edof, mesh.dofs_per_node, mesh.el_type, None, draw_elements=False, draw_undisplaced_mesh=False, title="Example 6 - Effective stress")
```

<div class="output display_data">

![](dbc6bee47c99874f9582a406707270429b399ca5.png)

</div>

</div>

<div class="cell code" execution_count="35" scrolled="true">

``` python
cfv.figure(fig_size=(10,10))
cfv.draw_displacements(a, coords, edof, mesh.dofs_per_node, mesh.el_type, draw_undisplaced_mesh=True, title="Example 06 - Displacements")
```

<div class="output display_data">

![](6bf983a88f85384c21fc1f0d50a01ff130efb2f6.png)

</div>

</div>

<div class="cell code">

``` python
```

</div>
