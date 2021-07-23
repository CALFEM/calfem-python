# CALFEM for Python

## Credits

CALFEM for Python is a collaborative effort. The main concept is of course the MATLAB version of CALFEM. The idea of a Python version of CALFEM started from the need to reduce the complexity of a course teaching the basics of how to implement a larger finite element code with user interface, solver and visualisation. We had already implemented the user interface and visualisation parts in Python, but the solver was implemented as a Fortran module. A lot of time was spent in the integration of the Fortran-parts instead of being able to focus on the actual design of the finite element implementation.

The first implementation of CALFEM for Python started in a master thesis work by Andreas Ottosson in 2010. In this thesis the ideas on how to implement the MATLAB version of CALFEM in Python using NumPy was studied. It was shown that it was quite possible to implement a version of CALFEM that used Python instead. About 80 percent of the core CALFEM was converted in this work.

During 2012-13 mesh generation was added using the GMSH-tool in a master thesis work by Andreas Edholm. This is now the default mesh generator in CALFEM. 

In 2017 CALFEM for Python was used in the visualisation tool VisCon for studying 2D earth consolidation. The Software was developed by Karin Forsman. This work lead to many improvements in the packaging of CALFEM for Python. 

The development effort is coordinated by Jonas Lindemann, which has also contributed a lot of modules and overall structure of the project. 

## Manuals

Original manual: [manual.pdf](https://github.com/CALFEM/calfem-python/raw/master/manual.pdf)

Manual for with improved mesh: [manual-mesh-module.pdf](https://github.com/CALFEM/calfem-python/raw/master/manual-mesh-module.pdf)

## Background

The computer program CALFEM is written for the software MATLAB and is an
interactive tool for learning the finite element method. CALFEM is an abbreviation
of ”Computer Aided Learning of the Finite Element Method” and been developed
by the Division of Structural Mechanics at Lund University since the late 70’s.

## Why CALFEM for Python?

Unlike MATLAB, which have expensive licenses, Python is free to use and distribute
both for personal and commercial use. An implementation to Py

## References

* Forsman, K, 2017. VisCon: Ett visualiseringsverktyg för tvådimensionell konsolidering i undervisningssammanhang - http://www.byggmek.lth.se/fileadmin/byggnadsmekanik/publications/tvsm5000/web5225.pdf 

* Edholm, A., 2013. Meshing and visualisation routines in the Python version of CALFEM.  - http://www.byggmek.lth.se/fileadmin/byggnadsmekanik/publications/tvsm5000/web5187.pdf 

* Ottosson, A., 2010. Implementation of CALFEM for Python - http://www.byggmek.lth.se/fileadmin/byggnadsmekanik/publications/tvsm5000/web5167.pdf 

