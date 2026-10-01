# Using CALFEM for Python

For the main introduction to CALFEM, the finite element workflow, and
the standard function groups, see the
[CALFEM manual](https://calfem-python-manual.readthedocs.io/en/latest/introduction.html).
This page only summarizes the Python package layout and points to the
local Python-specific documentation.

## What is CALFEM for Python?

- Subset of CALFEM routines implemented in Python
- Using NumPy for matrices
- Additional mesh generation routines supporting Triangle and GMSH
- Plotting with Matplotlib, plotly (interactive) and vedo (3D)

## CALFEM Python modules

- **calfem.core**
  - Element routines
  - System routines
- **calfem.utils**
  - I/O routines
  - Misc. routines
- **calfem.geometry**
  - Routines for defining problem geometry used for input in mesh
    generation
- **calfem.mesh**
  - Mesh generation routines
- **calfem.vis/calfem.vis_mpl**
  - Routines for visualising geometry, meshes and results with
    Matplotlib. calfem.vis is an alias for calfem.vis_mpl.
- **calfem.vis_plotly**
  - The same routines as calfem.vis_mpl with interactive plotly
    figures, e.g. for Jupyter notebooks.
- **calfem.vis_vedo**
  - Visualisation of 3D models with vedo.
- **calfem.io**
  - Export of meshes and results to ParaView and other tools
    (via meshio), saving and loading of meshes and geometries.
- **calfem.vis_visvis**
  - The former visvis based calfem.vis. Deprecated, visvis is no
    longer maintained.

## Examples

The CALFEM manual contains the
[canonical worked examples](https://calfem-python-manual.readthedocs.io/en/latest/examples.html).
This documentation keeps only the Python-specific
[mesh examples](mesh_examples.md) locally.

For background on the standard element, system, graphics, utility, and
matrix functions, use the manual together with the local
[function reference](calfem_reference.md).
