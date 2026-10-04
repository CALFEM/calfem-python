# CALFEM for Python

## Documentation

[https://calfem-for-python.readthedocs.io/en/latest/](https://calfem-for-python.readthedocs.io/en/latest/)

## Manuals

Original manual: [manual.pdf](https://github.com/CALFEM/calfem-python/tree/master/reports/manual.pdf)

Manual for with improved mesh: [manual-mesh-module.pdf](https://github.com/CALFEM/calfem-python/tree/master/reports/manual-mesh-module.pdf)

## Background

The computer program CALFEM is written for the software MATLAB and is an interactive tool for learning the finite element method. CALFEM is an abbreviation
of ”Computer Aided Learning of the Finite Element Method” and been developed by the Division of Structural Mechanics at Lund University since the late 70’s.

## Why CALFEM for Python?

While both the MATLAB and Python versions of CALFEM are open-source (MIT Licensed), the key difference lies in the environments they operate in. MATLAB is not open-source and requires expensive licenses (for commercial use), which can be a barrier for many users. In contrast, Python is a free, open-source programming language, making CALFEM for Python more accessible and cost-effective for a broader audience, including those in academic, personal, or commercial settings. 

CALFEM for Python is released under the MIT license, which enables its use in open-source as well as commercial projects.

## Installation

Install CALFEM for Python using

```bash
pip install calfem-python
```

This also installs the required dependencies NumPy, SciPy, Matplotlib, tabulate and the Gmsh Python module used for mesh generation. No separate Gmsh installation is needed.

## Visualisation

`calfem.vis` (an alias for `calfem.vis_mpl`) draws geometries, meshes and results using Matplotlib. `calfem.vis_plotly` provides the same functions as interactive plotly figures. The former visvis based module is available as `calfem.vis_visvis`, but is deprecated.

## Optional dependencies

Some modules need additional packages, which can be installed as extras:

| Extra | Installs | Needed for |
| --- | --- | --- |
| `plotly` | plotly | Interactive plots with `calfem.vis_plotly` |
| `io` | meshio | Mesh import/export with `calfem.io` |
| `vedo` | vedo | 3D visualisation with `calfem.vis_vedo` |
| `vtk` | vtk | Visualisation with `calfem.vis_vtk` |
| `pyvtk` | pyvtk | VTK export with `calfem.utils.export_vtk_stress` |
| `qtpy` | qtpy | The geometry editor `calfem.editor`, using an already installed Qt binding |
| `pyside6` | qtpy, PySide6 | The geometry editor `calfem.editor` with PySide6 |
| `pyqt6` | qtpy, PyQt6 | The geometry editor `calfem.editor` with PyQt6 |
| `visvis` | visvis | The deprecated `calfem.vis_visvis` |

For example:

```bash
pip install calfem-python[plotly,io]
```

When using conda, install the Qt binding from conda-forge instead (e.g. `conda install -c conda-forge pyside6`) together with the `qtpy` extra. Mixing a pip installed Qt binding with Qt libraries from conda can fail with DLL load errors.

## References

* Forsman, K, 2017. VisCon: Ett visualiseringsverktyg för tvådimensionell konsolidering i undervisningssammanhang - http://www.byggmek.lth.se/fileadmin/byggnadsmekanik/publications/tvsm5000/web5225.pdf 

* Edholm, A., 2013. Meshing and visualisation routines in the Python version of CALFEM.  - http://www.byggmek.lth.se/fileadmin/byggnadsmekanik/publications/tvsm5000/web5187.pdf 

* Ottosson, A., 2010. Implementation of CALFEM for Python - http://www.byggmek.lth.se/fileadmin/byggnadsmekanik/publications/tvsm5000/web5167.pdf 

* Eriksson, K, 2021. CALFEM Geometry Editor - An interactive geometry editor for CALFEM

* Åmand, A, 2022. Development of visualisation functions for CALFEM for Python


