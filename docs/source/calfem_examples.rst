Using CALFEM for Python
=======================

What is CALFEM?
---------------

CALFEM is an interactive computer program for teaching the finite element method (FEM). The name CALFEM is an abbreviation of "Computer Aided Learning of the Finite Element Method". The program can be used for different types of structural mechanics problems and field problems.

CALFEM, the program and its built-in philosophy, have been developed at the Division of Structural Mechanics, Lund University, starting in the late 70's. Many coworkers, former and present, have been engaged in the development at different stages.

What is CALFEM for Python?
--------------------------
* Subset of CALFEM routines implemented in Python
* Using NumPy for matrices
* Additional mesh generation routines supporting Triangle and GMSH
* Plotting with Matplotlib and visvis

CALFEM Python modules
---------------------

* **calfem.core**  

  * Element routines
  * System routines
* **calfem.utils**  

  * I/O routines
  * Misc. routines
* **calfem.geometry**  

  * Routines for defining problem geometry used for input in mesh generation
* **calfem.mesh**  

  * Mesh generation routines  
* **calfem.vis/calfem.vis_mpl**  

  * Routines for visualising geometry, meshes and results.

Examples
--------
The example codes show what CALFEM can do for you. The examples are divided into two:    

- Numerical examples

- Mesh examples

The next is tutorial on using Calfem for Python for numerical finite element, i.e., solving FEM equation to obtain nodal displacements given loading forces and 
stiffness matrix. The example can be found in `examples` directories both on 
calfem-python root directory (for .py files) and docs directory (for .ipynb files).
