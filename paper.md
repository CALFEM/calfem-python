---
title: 'CALFEM for Python: An Educational Package for Teaching the Finite Element Method'
tags:
  - Python
  - finite element method
  - structural mechanics
  - computational mechanics
  - engineering education
  - numerical methods
authors:
  - name: Jonas Lindemann
    orcid: 0009-0001-9799-4640  # Replace with your actual ORCID
    affiliation: 1
  - name: Per-Erik Austrell
    orcid: 0000-0002-1234-5678  # Replace with your actual ORCID
    affiliation: 1
  - name: Håkan Carlsson
    orcid: 0000-0003-8765-4321  # Replace with your actual ORCID
    affiliation: 1
  - name: Ola Dahlblom
    orcid: 0000-0002-2345-6789  # Replace with your actual ORCID
    affiliation: 1
  - name: Susanne Heyeden
    orcid: 0000-0001-3456-7890  # Replace with your actual ORCID
    affiliation: 1
  - name: Anders Olsson
    orcid: 0000-0002-4567-8901  # Replace with your actual ORCID
    affiliation: 1
  - name: Karl-Gunnar Olsson
    orcid: 0000-0003-5678-9012  # Replace with your actual ORCID
    affiliation: 1
  # Add other key contributors here
  # - name: Contributor Name
  #   orcid: 0000-0000-0000-0000
  #   affiliation: 1
affiliations:
  - name: Division of Structural Mechanics, Lund University, Sweden
    index: 1
date: 18 October 2025
bibliography: paper.bib
---

# Summary

CALFEM for Python is an open-source software package designed for teaching and learning the finite element method (FEM) in structural mechanics and solid mechanics courses. Originally developed as a MATLAB toolkit at Lund University since the late 1970s, CALFEM (Computer Aided Learning of the Finite Element Method) has been ported to Python to provide students and educators with a freely accessible, transparent implementation of fundamental FEM concepts. The package includes core finite element routines, mesh generation capabilities via GMSH and Triangle, and comprehensive visualization tools using matplotlib and visvis. CALFEM for Python bridges the gap between theoretical FEM concepts and practical implementation, enabling students to understand, modify, and extend finite element algorithms without the cost barriers of commercial software or proprietary development environments.

# Statement of Need

The finite element method is a cornerstone of modern engineering analysis, yet teaching FEM effectively presents several challenges. Commercial FEM software packages, while powerful, often act as "black boxes" that obscure the underlying algorithms from students. MATLAB-based teaching tools like the original CALFEM provide transparency but require expensive licenses that may be inaccessible to students and institutions, particularly in developing countries. Furthermore, as universities increasingly adopt Python for engineering education due to its open-source nature and widespread use in industry and research, there is a critical need for FEM teaching tools that align with this pedagogical shift.

CALFEM for Python addresses these needs by providing:

- **Accessibility**: Free and open-source (MIT licensed) software that eliminates cost barriers for students worldwide
- **Pedagogical Clarity**: Clean, readable Python implementations that allow students to examine and understand FEM algorithms at the code level
- **Practical Integration**: Seamless integration with the scientific Python ecosystem (NumPy, SciPy, matplotlib), preparing students for modern computational workflows
- **Comprehensive Functionality**: Coverage of essential FEM concepts including element formulation, assembly procedures, constraint handling, and post-processing
- **Active Learning**: Extensible codebase that encourages students to modify and experiment with FEM implementations

The package has been specifically designed for graduate and postgraduate courses in structural mechanics and has been developed in direct response to Lund University's transition to Python-based computational education. Unlike research-focused FEM packages such as FEniCS or commercial packages like ABAQUS, CALFEM for Python prioritizes educational transparency and ease of understanding over computational efficiency or advanced features.

# Learning Objectives and Content

CALFEM for Python supports the following learning objectives in finite element courses:

1. **Understanding FEM Fundamentals**: Students learn element formulation, stiffness matrix assembly, boundary condition application, and solution procedures through transparent code implementations
2. **Hands-on Implementation**: Students can modify element routines, implement new element types, and experiment with different solution strategies
3. **Mesh Generation**: Integration with GMSH and Triangle teaches students about discretization quality and its impact on solution accuracy
4. **Visualization and Interpretation**: Built-in plotting capabilities help students visualize meshes, deformations, stresses, and other field quantities
5. **Research Preparation**: Graduate students gain skills in algorithm development and FEM programming that transfer to research applications

## Package Components

The package is organized into several modules:

- **`calfem.core`**: Core FEM routines including element stiffness matrices (bar, beam, plane stress/strain, solid elements), assembly functions, equation solvers, and boundary condition handlers
- **`calfem.geometry`**: Geometry definition tools for creating computational domains
- **`calfem.mesh`**: Mesh generation interfaces to GMSH and Triangle for structured and unstructured meshes
- **`calfem.vis_mpl`**: Matplotlib-based visualization for 2D and 3D results
- **`calfem.utils`**: Utility functions for data handling and result extraction

## Educational Materials

CALFEM for Python includes:

- Comprehensive API documentation hosted on ReadTheDocs
- Multiple worked examples covering various problem types (structural frames, heat transfer, solid mechanics)
- Jupyter notebook tutorials demonstrating typical analysis workflows
- PDF manuals detailing theory and implementation (including the original CALFEM manual and expanded mesh generation guide)

# Instructional Design and Experience

## Pedagogical Approach

CALFEM for Python follows a constructivist learning approach where students build understanding by actively working with FEM implementations. The package design emphasizes:

- **Code Readability**: Functions use descriptive names and clear structure that mirrors mathematical formulations
- **Gradual Complexity**: Examples progress from simple 1D bar elements to complex 3D solid mechanics problems
- **Exploratory Learning**: Students are encouraged to modify examples, change parameters, and observe results
- **Error Transparency**: Clear error messages and debugging capabilities help students understand when and why solutions fail

## Teaching Experience

CALFEM for Python has been used at Lund University for teaching finite element methods to graduate and postgraduate students in civil engineering, mechanical engineering, and engineering mechanics since 2010. The package has proven particularly effective for:

- **Bridging Theory and Practice**: Students connect theoretical lectures on FEM with hands-on coding exercises
- **Thesis Projects**: Multiple master's thesis projects have extended the package's capabilities (mesh generation, visualization, specialized applications)
- **Research Training**: Graduate students use the package as a foundation for developing custom FEM codes for their research

## Adoption by Other Instructors

The package is designed for straightforward adoption by other instructors:

- **Easy Installation**: Available via pip (`pip install calfem-python`)
- **Minimal Dependencies**: Relies primarily on standard scientific Python libraries
- **Flexible Integration**: Can be used for homework assignments, laboratory exercises, or project work
- **Example Repository**: Extensive examples provide starting points for assignments and projects

Instructors can adapt the existing examples to their specific course needs or guide students in developing custom applications. The open-source nature allows customization for specialized topics (e.g., non-linear analysis, dynamics) without license restrictions.

# Comparison to Similar Educational Software

Several other FEM packages exist for education and research:

- **FEniCS**: A powerful automated FEM framework that excels at complex PDEs but has a steeper learning curve and abstracts away many implementation details
- **SfePy**: A research-oriented FEM package with extensive capabilities but less focus on pedagogical clarity
- **GetFEM++**: Primarily C++-based with Python bindings, targeting research more than education
- **deal.II**: Excellent C++ FEM library with good documentation but requires C++ proficiency
- **PyFEM**: Similar educational focus but less mature and less comprehensive in scope

CALFEM for Python fills a specific niche by providing undergraduate/graduate-level transparency with comprehensive documentation specifically designed for structural mechanics education. Its MATLAB heritage means it aligns well with classical FEM textbooks, while its Python implementation prepares students for modern computational practices.

# Acknowledgements

The development of CALFEM for Python builds upon decades of work on the original MATLAB version by the Division of Structural Mechanics at Lund University. Significant contributions to the Python port have been made through master's thesis projects, including work on mesh generation (Edholm 2013), initial Python implementation (Ottosson 2010), visualization (Åmand 2022), and geometry editors (Eriksson 2021). We acknowledge the continued support of Lund University and the contribution of students and faculty who have used and improved the package.

# References