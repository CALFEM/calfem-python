#!/bin/env python
# -*- coding: iso-8859-15 -*-
"""
This is a utility module for the CALFEM Python library. It contains various utility functions that is used throughout the library. It includes functions for reading and writing files, applying boundary conditions, displaying messages, and exporting data in different formats.
"""

import os
import sys
import pickle
import scipy.io

import numpy as np
import calfem.core as cfc
import logging as cflog
import tabulate as tab

have_pyvtk = True

try:
    import pyvtk as vtk
except: 
    have_pyvtk = False

haveMatplotLib = True
haveMlab = True

have_mlab = haveMlab
have_matplotlib = haveMatplotLib

have_ipython = True

try:
    from IPython.core.display import display, HTML
except:
    have_ipython = False



def error(msg):
    cflog.error(" "+msg)

def info(msg):
    cflog.info(" "+msg)

def set_debug_level(level):
    cflog.getLogger().setLevel(level)

def type_of_script():
    try:
        ipy_str = str(type(get_ipython()))
        if 'zmqshell' in ipy_str:
            return 'jupyter'
        if 'terminal' in ipy_str:
            return 'ipython'
    except:
        return 'terminal'

def disp(msg):
    if type_of_script() == 'jupyter':
        display(HTML(f"{msg}"))
    else:
        print(msg)

def str_disp(msg):
    return f"{msg}\n"

def disp_par(msg):
    if type_of_script() == 'jupyter':
        display(HTML(f"<p>{msg}</p>"))
    else:
        print(f"\nmsg\n")

def str_disp_par(msg):
    return f"\n{msg}\n"

def disp_bold(msg):
    if type_of_script() == 'jupyter':
        display(HTML(f"<b>{msg}</b>"))
    else:
        print(f"**{msg}**")

def str_disp_bold(msg):
    return f"**{msg}**"
    
def disp_bold_par(msg):
    if type_of_script() == 'jupyter':
        display(HTML(f"<p><b>{msg}</b></p>"))
    else:
        print(f"\n**{msg}**\n")

def str_disp_bold_par(msg):
    return f"\n**{msg}**\n"

def disp_h1(msg):
    if type_of_script() == 'jupyter':
        display(HTML(f"<h1>{msg}</h1>"))
    else:
        print(f"\n# {msg}\n")

def str_disp_h1(msg):    
    return f"\n# {msg}\n"

def disp_h2(msg):
    if type_of_script() == 'jupyter':
        display(HTML(f"<h2>{msg}</h2>"))
    else:
        print(f"\n## {msg}\n")

def str_disp_h2(msg):
    return f"\n## {msg}\n"

def disp_h3(msg):
    if type_of_script() == 'jupyter':
        display(HTML(f"<h3>{msg}</h3>"))
    else:
        print(f"\n### {msg}\n")

def str_disp_h3(msg):
    return f"\n### {msg}\n"

def disp_array(a, headers=[], fmt=".4e", tablefmt="psql", showindex=False):
    """
    Print a numpy array in a nice way.
    """
    if type_of_script() == 'jupyter':
        display(tab.tabulate(np.asarray(a), tablefmt="html", floatfmt=".4e", showindex=showindex, headers=headers))
    else:
        print(tab.tabulate(np.asarray(a), tablefmt=tablefmt, floatfmt=fmt, showindex=showindex, headers=headers))

def str_disp_array(a, headers=[], fmt=".4e", tablefmt="psql", showindex=False):
    """
    Return a numpy array in a nice way as a string.
    """
    return tab.tabulate(np.asarray(a), tablefmt=tablefmt, floatfmt=fmt, showindex=showindex, headers=headers)

class ElementProperties(object):
    def __init__(self):
        self.ep = {}
        self.attributes = {}

    def add(self, markerId, ep=[]):
        if not markerId in self.ep:
            self.ep[markerId] = ep

    def addAttribute(self, markerId, name, value):
        if not markerId in self.attributes:
            self.attributes[markerId] = {}
            self.attributes[markerId][name] = value


def enable_logging(def_level=cflog.INFO):
    cflog.basicConfig(
        format='%(asctime)s:%(levelname)s:%(message)s', level=def_level)

enableLogging = enable_logging


def disable_logging():
    cflog.basicConfig(
        format='%(asctime)s:%(levelname)s:%(message)s', level=cflog.NOTSET)

disableLogging = disable_logging


def read_int(f):
    """
    Read a row from file, f, and return a list of integers.
    """
    return list(map(int, f.readline().split()))

readInt = read_int


def read_float(f):
    """
    Read a row from file, f, and return a list of floats.
    """
    return list(map(float, f.readline().split()))

readFloat = read_float


def read_single_int(f):
    """
    Read a single integer from a row in file f. All other values on row are discarded.
    """
    return readInt(f)[0]

readSingleInt = read_single_int


def read_single_float(f):
    """
    Read a single float from a row in file f. All other values on row are discarded.
    """
    return readFloat(f)[0]

readSingleFloat = read_single_float


def write_single_float(f, floatValue):
    f.write("%g\n" % floatValue)

writeSingleFloat = write_single_float


def write_single_int(f, intValue):
    f.write("%d\n" % intValue)

writeSingleInt = write_single_int


def write_float_list(f, floatList):
    for floatValue in floatList:
        f.write("%g " % floatValue)
    f.write("\n")

writeFloatList = write_float_list

def write_int_list(f, intList):
    for intValue in intList:
        f.write("%d " % intValue)
    f.write("\n")

writeIntList = write_int_list

def which(filename):
    """
    Return complete path to executable given by filename.
    """
    if not ('PATH' in os.environ) or os.environ['PATH'] == '':
        p = os.defpath
    else:
        p = os.environ['PATH']

    pathlist = p.split(os.pathsep)
    pathlist.insert(0, ".")
    pathlist.insert(0, "/bin")
    pathlist.insert(0, "/usr/bin")
    pathlist.insert(0, "/opt/local/bin")
    pathlist.insert(0, "/usr/local/bin")
    pathlist.insert(0, "/Applications/Gmsh.app/Contents/MacOS")

    # Add paths from site-packages

    for path in sys.path:
        if "site-packages" in path:
            pathlist.insert(0, path)

    for path in pathlist:
        f = os.path.join(path, filename)

        if os.access(f, os.X_OK):
            return f

    return None


def apply_bc(boundaryDofs, bcPrescr, bcVal, marker, value=0.0, dimension=0):
    """
    Apply boundary condition to bcPresc and bcVal matrices. For 2D problems
    with 2 dofs per node.
    
    Parameters
    ----------
    boundaryDofs : dict
        Dictionary with boundary dofs.
    bcPresc : array_like
        1-dim integer array containing prescribed dofs.
    bcVal : array_like
        1-dim float array containing prescribed values.
    marker : int
        Boundary marker to assign boundary condition.
    value : float, optional
        Value to assign boundary condition.
        If not given 0.0 is assigned.
    dimension : int, optional
        dimension to apply bc. 0 - all, 1 - x, 2 - y

    Returns
    -------
    bcPresc : array_like
        Updated 1-dim integer array containing prescribed dofs.
    bcVal : array_like
        Updated 1-dim float array containing prescribed values.
    """

    if marker in boundaryDofs:
        if (dimension == 0):
            bcAdd = np.array(boundaryDofs[marker])
            bcAddVal = np.ones([np.size(bcAdd)])*value
        elif dimension in [1, 2]:
            bcAdd = boundaryDofs[marker][(dimension-1)::2]
            bcAddVal = np.ones([np.size(bcAdd)])*value
        else:
            print("Error: wrong dimension, ", dimension)

        newBcPrescr, prescrIdx = np.unique(
            np.hstack([bcPrescr, bcAdd]), return_index=True)
        newBcVal = np.hstack([bcVal, bcAddVal])[prescrIdx]

        return newBcPrescr, newBcVal
    else:
        print("Error: Boundary marker", marker, "does not exist.")

applybc = apply_bc


def apply_bc_3d(boundaryDofs, bcPrescr, bcVal, marker, value=0.0, dimension=0):
    """
    Apply boundary condition to bcPresc and bcVal matrices. For 3D problems
    with 3 dofs per node.
    
    Parameters
    ----------
    boundaryDofs : dict
        Dictionary with boundary dofs.
    bcPrescr : array_like
        1-dim integer array containing prescribed dofs.
    bcVal : array_like
        1-dim float array containing prescribed values.
    marker : int
        Boundary marker to assign boundary condition.
    value : float, optional
        Value to assign boundary condition.
        If not given 0.0 is assigned.
    dimension : int, optional
        dimension to apply bc. 0 - all, 1 - x, 2 - y,
        3 - z

    Returns
    -------
    bcPrescr : array_like
        Updated 1-dim integer array containing prescribed dofs.
    bcVal : array_like
        Updated 1-dim float array containing prescribed values.
    """

    if marker in boundaryDofs:
        if (dimension == 0):
            bcAdd = np.array(boundaryDofs[marker])
            bcAddVal = np.ones([np.size(bcAdd)])*value
        elif dimension in [1, 2, 3]:
            bcAdd = boundaryDofs[marker][(dimension-1)::3]
            bcAddVal = np.ones([np.size(bcAdd)])*value
        else:
            print("Error: wrong dimension, ", dimension)

        newBcPrescr, prescrIdx = np.unique(
            np.hstack([bcPrescr, bcAdd]), return_index=True)
        newBcVal = np.hstack([bcVal, bcAddVal])[prescrIdx]

        return newBcPrescr, newBcVal
    else:
        print("Error: Boundary marker", marker, "does not exist.")

applybc3D = apply_bc_3d


def apply_bc_node(nodeIdx, dofs, bcPrescr, bcVal, value=0.0, dimension=0):
    """
    Apply boundary conditions to a specific node.
    This function adds boundary condition prescriptions and values for a given node
    to existing boundary condition arrays.

    Parameters
    ----------
    nodeIdx : int
        Index of the node to apply boundary conditions to.
    dofs : array_like
        Degrees of freedom array. Can be 1D (for single DOF per node) or 2D 
        (for multiple DOFs per node).
    bcPrescr : array_like
        Existing array of prescribed boundary condition DOF indices.
    bcVal : array_like
        Existing array of prescribed boundary condition values.
    value : float, optional
        Value to prescribe for the boundary condition. Default is 0.0.
    dimension : int, optional
        Dimension/direction to apply BC. If 0, applies to all DOFs of the node.
        If 1, 2, or 3, applies to specific dimension (1-indexed). Default is 0.
    
    Returns
    -------
    tuple of numpy.ndarray
        A tuple containing:
        - Updated prescribed DOF indices array (bcPrescr concatenated with new DOFs)
        - Updated prescribed values array (bcVal concatenated with new values)

    Notes
    -----
    When dimension=0, boundary conditions are applied to all degrees of freedom
    for the specified node. When dimension is 1, 2, or 3, the boundary condition
    is applied only to that specific dimension (using 1-based indexing).

    Examples
    --------
    >>> # Apply BC to all DOFs of node 5 with value 0.0
    >>> bc_dofs, bc_vals = apply_bc_node(5, dofs, [], [], 0.0, 0)
    >>> 
    >>> # Apply BC to x-direction (dimension 1) of node 10 with value 5.0
    >>> bc_dofs, bc_vals = apply_bc_node(10, dofs, bc_dofs, bc_vals, 5.0, 1)
    """

    if (dimension == 0):
        bcAdd = np.asarray(dofs[nodeIdx])
        bcAddVal = np.ones([np.size(bcAdd)])*value
    elif dimension in [1, 2, 3]:
        bcAdd = np.asarray(dofs[nodeIdx, dimension-1])
        bcAddVal = np.ones([np.size(bcAdd)])*value
    else:
        print("Error: wrong dimension, ", dimension)

    return np.hstack([bcPrescr, bcAdd]), np.hstack([bcVal, bcAddVal])

applybcnode = apply_bc_node

def apply_force_node(nodeIdx, dofs, f, value=0.0, dimension=0):
    """
    Apply a force to a specific node in the finite element model.
    This function adds a force value to the global force vector at the degrees of freedom
    corresponding to a specified node. The force can be applied to all DOFs of the node
    or to a specific dimension.
    Parameters
    ----------
    nodeIdx : int
        Index of the node where the force is to be applied.
    dofs : array_like
        Degrees of freedom array that maps nodes to their DOF indices in the global system.
        Can be 1D (for single DOF per node) or 2D (for multiple DOFs per node).
    f : array_like
        Global force vector where the force will be added.
    value : float, optional
        Magnitude of the force to be applied. Default is 0.0.
    dimension : int, optional
        Specific dimension/DOF to apply the force to. If 0, applies to all DOFs of the node.
        If 1 or higher, applies to the specified dimension (1-indexed). Default is 0.
    Notes
    -----
    - When dimension=0, the force is applied to all DOFs of the node (assumes 1D dofs array)
    - When dimension>=1, the force is applied to the specific dimension of the node
      (assumes 2D dofs array with shape [node, dimension])
    - The dimension parameter uses 1-based indexing (dimension=1 corresponds to first DOF)
    Examples
    --------
    >>> # Apply force to all DOFs of node 5
    >>> apply_force_node(5, dofs, f, value=100.0, dimension=0)
    >>> # Apply force to x-direction (dimension 1) of node 3
    >>> apply_force_node(3, dofs, f, value=50.0, dimension=1)
    """

    if (dimension == 0):
        f[dofs[nodeIdx]] += value
    elif (dimension == 1):
        f[dofs[nodeIdx, dimension-1]] += value
    else:
        f[dofs[nodeIdx, dimension-1]] += value

applyforcenode = apply_force_node


def apply_force(boundaryDofs, f, marker, value=0.0, dimension=0):
    """
    Apply boundary force to f matrix. The value is
    added to all boundaryDofs defined by marker. Applicable
    to 2D problems with 2 dofs per node.
    
    Parameters:
    
        boundaryDofs        Dictionary with boundary dofs.
        f                   force matrix.
        marker              Boundary marker to assign boundary condition.
        value               Value to assign boundary condition.
                            If not given 0.0 is assigned.
        dimension           dimension to apply force. 0 - all, 1 - x, 2 - y
                            
    """

    if marker in boundaryDofs:
        if dimension == 0:
            f[np.asarray(boundaryDofs[marker])-1] += value
        elif dimension in [1, 2]:
            f[np.asarray(boundaryDofs[marker][(dimension-1)::2])-1] += value
        else:
            print("Error: The dimension, ", dimension, ", is invalid")
    else:
        print("Error: Boundary marker", marker, "does not exist.")

applyforce = apply_force


def apply_traction_linear_element(boundaryElements, coords, dofs, F, marker, q):
    """
    Apply traction on part of boundary with marker.

    q is added to all boundaryDofs defined by marker. Applicable
    to 2D problems with 2 dofs per node. The function works with linear
    line elements. (elm-type 1 in GMSH).

    Parameters
    ----------
    boundaryElements : dict
        Dictionary with boundary elements, the key is a marker and the values are lists of elements.
    coords : array_like
        Coordinates matrix
    dofs : array_like
        Dofs matrix
    F : array_like
        force matrix.
    marker : int
        Boundary marker to assign boundary condition.
    q : array_like
        Value to assign boundary condition.
        shape = [qx qy] in global coordinates

    """
    if marker not in boundaryElements:
        print("Error: Boundary marker", marker, "does not exist.")
        return
    for element in boundaryElements[marker]:
        if element['elm-type'] != 1:
            print("Error: Wrong element type.")
            return

    q = np.matrix(q).T

    # Integration points and weights:
    Xi = [-1/np.sqrt(3), 1/np.sqrt(3)]
    W = [1, 1]

    # Shape functions:
    def N1(x): return 1-(1+x)/2
    def N2(x): return (1+x)/2

    for element in boundaryElements[marker]:
        # Loop through integration points:
        f = np.zeros([4, 1])
        for xi, w in zip(Xi, W):
            N = np.matrix([[N1(xi),      0,  N2(xi),       0],
                           [0, N1(xi),        0, N2(xi)]])
            # The minus one is since the nodes in node-number-list start at 1...
            coord = coords[np.array(element['node-number-list'])-1]
            v1 = coord[0, :]
            v2 = coord[1, :]
            J = np.linalg.norm(v1-v2) / 2
            f += w * N.T * q * J

        # Minus one since dofs start at 1...
        idx = dofs[np.array(element['node-number-list'])-1, :].flatten()-1
        F[idx] += f

applyTractionLinearElement = apply_traction_linear_element


def apply_force_3d(boundaryDofs, f, marker, value=0.0, dimension=0):
    """
    Apply boundary force to f matrix for 3D problems.
    The value is added to all boundaryDofs defined by marker. Applicable
    to 3D problems with 3 degrees of freedom per node.
    Parameters
    ----------
    boundaryDofs : dict
        Dictionary with boundary degrees of freedom.
    f : numpy.ndarray
        Force matrix to be modified.
    marker : int or str
        Boundary marker to identify which boundary condition to apply.
    value : float, optional
        Value to add to the force matrix at specified boundary DOFs.
        Default is 0.0.
    dimension : int, optional
        Dimension to apply force:
        * 0 - all dimensions (default)
        * 1 - x-direction only
        * 2 - y-direction only  
        * 3 - z-direction only
    Notes
    -----
    If the specified marker does not exist in boundaryDofs, an error message
    is printed. If an invalid dimension is specified (not 0, 1, 2, or 3),
    an error message is printed.
    Examples
    --------
    >>> boundaryDofs = {1: [1, 2, 3, 4, 5, 6]}
    >>> f = np.zeros(6)
    >>> apply_force_3d(boundaryDofs, f, 1, value=100.0, dimension=1)
    # Applies force of 100.0 in x-direction to DOFs 1, 4
    """

    if marker in boundaryDofs:
        if dimension == 0:
            f[np.asarray(boundaryDofs[marker])-1] += value
        elif dimension in [1, 2, 3]:
            f[np.asarray(boundaryDofs[marker][(dimension-1)::3])-1] += value
        else:
            print("Error: The dimension, ", dimension, ", is invalid")
    else:
        print("Error: Boundary marker", marker, "does not exist.")

applyforce3D = apply_force_3d


def apply_force_total(boundaryDofs, f, marker, value=0.0, dimension=0):
    """
    Apply boundary force to f matrix. Total force, value, is
    distributed over all boundaryDofs defined by marker. Applicable
    to 2D problems with 2 dofs per node.
    
    Parameters
    ----------
    boundaryDofs : dict
        Dictionary with boundary dofs.
    f : array_like
        Force matrix.
    marker : int
        Boundary marker to assign boundary condition.
    value : float, optional
        Total force value to assign boundary condition.
        If not given 0.0 is assigned.
    dimension : int, optional
        Dimension to apply force. 0 - all, 1 - x, 2 - y
                            
    """

    if marker in boundaryDofs:
        if dimension == 0:
            nDofs = len(boundaryDofs[marker])
            valuePerDof = value / nDofs
            f[np.asarray(boundaryDofs[marker])-1] += valuePerDof
        elif dimension in [1, 2]:
            nDofs = len(boundaryDofs[marker][(dimension-1)::2])
            valuePerDof = value / nDofs
            f[np.asarray(boundaryDofs[marker][(dimension-1)::2]) -
              1] += valuePerDof
        else:
            print("Error: The dimension, ", dimension, ", is invalid")
    else:
        print("Error: Boundary marker", marker, "does not exist.")

applyforcetotal = apply_force_total


def apply_force_total_3d(boundaryDofs, f, marker, value=0.0, dimension=0):
    """
    Apply boundary force to f matrix. Total force, value, is
    distributed over all boundaryDofs defined by marker. Applicable
    to 3D problems with 3 dofs per node.
    
    Parameters
    ----------
    boundaryDofs : dict
        Dictionary with boundary dofs.
    f : array_like
        Force matrix.
    marker : int
        Boundary marker to assign boundary condition.
    value : float, optional
        Total force value to assign boundary condition.
        If not given 0.0 is assigned.
    dimension : int, optional
        Dimension to apply force. 0 - all, 1 - x, 2 - y,
        3 - z
                            
    """

    if marker in boundaryDofs:
        if dimension == 0:
            nDofs = len(boundaryDofs[marker])
            valuePerDof = value / nDofs
            f[np.asarray(boundaryDofs[marker])-1] += valuePerDof
        elif dimension in [1, 2, 3]:
            nDofs = len(boundaryDofs[marker][(dimension-1)::3])
            valuePerDof = value / nDofs
            f[np.asarray(boundaryDofs[marker][(dimension-1)::3]) -
              1] += valuePerDof
        else:
            print("Error: The dimension, ", dimension, ", is invalid")
    else:
        print("Error: Boundary marker", marker, "does not exist.")

applyforcetotal3D = apply_force_total_3d


def export_vtk_stress(filename, coords, topo, a=None, el_scalar=None, el_vec1=None, el_vec2=None):
    """
    Export mesh and results for a 2D stress problem.
    
    Parameters
    ----------
    filename : str
        Filename of vtk-file
    coords : numpy.ndarray
        Element coordinates
    topo : numpy.ndarray
        Element topology (not dof topology). mesh.topo.
    a : numpy.ndarray, optional
        Element displacements 2-dof
    el_scalar : list, optional
        Scalar values for each element
    el_vec1 : list, optional
        Vector value for each element
    el_vec2 : list, optional
        Vector value for each element
    """

    points = np.zeros([coords.shape[0], 3], dtype=np.float64)
    points[:,0:2] = coords
    points = points.tolist()
    polygons = (topo-1).tolist()

    displ = []

    point_data = None
    scalars = None
    vectors1 = None
    vectors2 = None
    cell_data = None

    if a is not None:
        for i in range(0, len(a), 2):
            displ.append([a[i].item(), a[i+1].item(), 0.0])

        point_data = vtk.PointData(vtk.Vectors(displ, name="displacements"))

    if el_scalar is not None:
        print("Adding cell scalars...")
        scalars = vtk.Scalars(el_scalar, name="scalar")
    if el_vec1 is not None:
        print("Adding cell vector 1...")
        vectors1 = vtk.Vectors(el_vec1, name="principal1")
    if el_vec2 is not None:
        print("Adding cell vector 2...")
        vectors2 = vtk.Vectors(el_vec2, name="principal2")

    if el_scalar is not None and el_vec1 is None and el_vec2 is None:
        print("Exporting celldata, el_scalar...")
        cell_data = vtk.CellData(scalars)
    if el_scalar is not None and el_vec1 is None and el_vec2 is not None:
        print("Exporting celldata, el_scalar, el_vec2...")
        cell_data = vtk.CellData(scalars, vectors2)
    if el_scalar is not None and el_vec1 is not None and el_vec2 is None:
        print("Exporting celldata, el_scalar, el_vec1...")
        cell_data = vtk.CellData(scalars, vectors1)
    if el_scalar is not None and el_vec1 is not None and el_vec2 is not None:
        print("Exporting celldata, el_scalar, el_vec1, el_vec2...")
        cell_data = vtk.CellData(scalars, vectors1, vectors2)
    if el_scalar is None and el_vec1 is None and el_vec2 is not None:
        print("Exporting celldata, el_vec2...")
        cell_data = vtk.CellData(vectors2)
    if el_scalar is None and el_vec1 is not None and el_vec2 is None:
        print("Exporting celldata, el_vec1...")
        cell_data = vtk.CellData(vectors1)
    if el_scalar is None and el_vec1 is not None and el_vec2 is None:
        print("Exporting celldata, el_vec1, el_vec2...")
        cell_data = vtk.CellData(vectors1, vectors2)

    structure = vtk.PolyData(points=points, polygons=polygons)

    if cell_data is not None and point_data is not None:
        print("VTK includes cell_data and point_data")
        vtk_data = vtk.VtkData(structure, cell_data, point_data)
    if cell_data is None and point_data is not None:
        print("VTK includes point_data")
        vtk_data = vtk.VtkData(structure, point_data)
    if cell_data is None and point_data is None:
        print("VTK includes only structure")
        vtk_data = vtk.VtkData(structure)

    vtk_data.tofile("exm6.vtk", "ascii")

def scalfact2(ex, ey, ed, rat=0.2):
    """
    Determine scale factor for drawing computational results, such as 
    displacements, section forces or flux.

    Parameters
    ----------
    ex : array_like
        Element node x coordinates
    ey : array_like  
        Element node y coordinates
    ed : array_like
        Element displacement matrix or section force matrix
    rat : float, optional
        Relation between illustrated quantity and element size. 
        Default is 0.2.

    Returns
    -------
    float
        Scale factor for drawing computational results
    """
    # nen:   number of element nodes
    # nel:   number of elements
    nen = -1
    if ex.shape != ey.shape:
        print("ex and ey shapes do not match.")
        return 1.0

    dlmax = 0.
    edmax = 1.

    if np.linalg.matrix_rank(ex) == 1:
        nen = ex.shape[0]
        nel = 1
        dxmax = max(ex.T.max(0)-ex.T.min(0))  # axis 0, return vector
        dymax = max(ey.T.max(0)-ey.T.min(0))
        dlmax = max(dxmax, dymax)
        edmax = abs(ed).max()
    else:
        nen = ex.shape[1]
        nel = ex.shape[0]
        dxmax = max(ex.T.max(0)-ex.T.min(0))
        dymax = max(ey.T.max(0)-ey.T.min(0))
        dlmax = max(dxmax, dymax)
        edmax = abs(ed).max()

    k = rat
    return k*dlmax/edmax


def load_geometry(name):
    """Loads a geometry from a file."""

    with open(name, 'rb') as file:
        test = pickle.load(file)
    return test


def save_geometry(g, name="unnamed_geometry"):
    """Save a geometry to file."""

    if not name.endswith(".cfg"):
        name = name + ".cfg"
    with open(name, 'wb') as file:
        pickle.dump(g, file)


def load_mesh(name):
    """Load a mesh from file."""

    with open(name, 'rb') as file:
        mesh = pickle.load(file)
    return mesh


def save_mesh(mesh, name="Untitled"):
    """Save a mesh to file."""

    if not name.endswith(".cfm"):
        name = name + ".cfm"
    with open(name, 'wb') as file:
        pickle.dump(mesh, file)


def save_arrays(coords, edof, dofs, bdofs, elementmarkers, boundaryElements, markerDict, name="unnamed_arrays"):
    """Save arrays to file."""

    if not name.endswith(".cfma"):
        name = name + ".cfma"
    with open(name, 'wb') as file:
        pickle.dump(coords, file)
        pickle.dump(edof, file)
        pickle.dump(dofs, file)
        #for key in bdofs.items():
        #    print(key, markerDict[key])
        pickle.dump(bdofs, file)
        pickle.dump(elementmarkers, file)
        pickle.dump(boundaryElements, file)
        pickle.dump(markerDict, file)


def load_arrays(name):
    """Load arrays from file."""

    with open(name, 'rb') as file:
        coords = pickle.load(file)
        edof = pickle.load(file)
        dofs = pickle.load(file)
        bdofs = pickle.load(file)
        elementmarkers = pickle.load(file)
        boundaryElements = pickle.load(file)
        markerDict = pickle.load(file)

    return coords, edof, dofs, bdofs, elementmarkers, boundaryElements, markerDict


def save_matlab_arrays(coords, edof, dofs, bdofs, elementmarkers, boundaryElements, markerDict, name="Untitled"):
    """Save arrays as MATLAB .mat files."""

    if not name.endswith(".mat"):
        name = name + ".mat"
    saveDict = {}
    saveDict["coords"] = coords.astype('double')

    # Convert to CALFEM Edof definition with element number as first index

    new_column = np.arange(1, np.size(edof, 0) + 1)[:, np.newaxis]
    edof = np.append(new_column, edof, axis=1)

    saveDict["edof"] = edof.astype('double')
    saveDict["dofs"] = dofs.astype('double')
    # bdofs = {str(k): v for k, v in bdofs.items()} # MATLAB struct needs keys as strings
    print(bdofs)
    print(markerDict)
    newBdof = {}
    for index, bdofs in bdofs.items():
        print(index, bdofs)
        if index == 0:
            newBdof["None"] = bdofs
        else:
            newBdof[markerDict[index]] = bdofs

    saveDict["bdofs"] = newBdof
    elementmarkers = np.asarray(elementmarkers)

    # To avoid problems with one indexing in MATLAB
    
    elementmarkers = elementmarkers + 1
    saveDict["elementmarkers"] = elementmarkers
    scipy.io.savemat(name, saveDict)

def calc_limits(coords):
    """Calculate max an min limits of 3d coordinates"""

    if coords.shape[0]>0:

        max_x = np.max(coords[:,0])
        min_x = np.min(coords[:,0])
        max_y = np.max(coords[:,1])
        min_y = np.min(coords[:,1])
        max_z = np.max(coords[:,2])
        min_z = np.min(coords[:,2])

        return max_x, min_x, max_y, min_y, max_z, min_z
    else:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

def calc_size(coords):
    """Calculate max and min sizes of 3d coordinates."""

    max_x, min_x, max_y, min_y, max_z, min_z = calc_limits(coords)
    lx = max_x - min_x
    ly = max_y - min_y
    lz = max_z - min_z

    return lx, ly, lz

def calc_beam_displ_limits(a, coords, edof, dofs):
    """
    Calculate max and min displacements for beams.

    Parameters
    ----------
    a : array_like
        Global displacement array with 6 dofs / node.
    coords : array_like
        Node coordinates.
    edof : array_like
        Beam topology.
    dofs : array_like
        Node dofs.

    Returns
    -------
    tuple
        Tuple containing (min_displ, max_displ).
    """

    if edof.shape[0]>0:

        ex, ey, ez = cfc.coord_extract(edof, coords, dofs)
        ed = cfc.extract_eldisp(edof, a)

        coords, topo, node_dofs = convert_to_node_topo(edof, ex, ey, ez, n_dofs_per_node=6, ignore_first=False)

        min_displ = 1e300
        max_displ = -1e300

        for el_topo in topo:
            d0 = np.array(a[node_dofs[el_topo[0]][:3]-1]).flatten()
            d1 = np.array(a[node_dofs[el_topo[1]][:3]-1]).flatten()

            l_d0 = np.linalg.norm(d0) 
            l_d1 = np.linalg.norm(d1) 

            if l_d0>max_displ:
                max_displ = l_d0
            if l_d1>max_displ:
                max_displ = l_d1
            if l_d0<min_displ:
                min_displ = l_d0
            if l_d1<min_displ:
                min_displ = l_d1

        return min_displ, max_displ

    else: 

        return 0.0, 0.0
    
def calc_bar_displ_limits(a, coords, edof, dofs):
    """
    Calculate max and min global displacements for bars.

    Parameters
    ----------
    a : array_like
        Global displacement array with 3 dofs / node.
    coords : array_like
        Node coordinates.
    edof : array_like
        Bar topology.
    dofs : array_like
        Node dofs.

    Returns
    -------
    tuple
        Tuple containing (min_displ, max_displ).
    """

    if edof.shape[0]>0:

        ex, ey, ez = cfc.coord_extract(edof, coords, dofs)
        ed = cfc.extract_eldisp(edof, a)

        coords, topo, node_dofs = convert_to_node_topo(edof, ex, ey, ez, n_dofs_per_node=3, ignore_first=False)

        min_displ = 1e300
        max_displ = -1e300

        for el_topo in topo:
            d0 = np.array(a[node_dofs[el_topo[0]][:3]-1]).flatten()
            d1 = np.array(a[node_dofs[el_topo[1]][:3]-1]).flatten()

            l_d0 = np.linalg.norm(d0)
            l_d1 = np.linalg.norm(d1)

            if l_d0>max_displ:
                max_displ = l_d0
            if l_d1>max_displ:
                max_displ = l_d1
            if l_d0<min_displ:
                min_displ = l_d0
            if l_d1<min_displ:
                min_displ = l_d1

        return min_displ, max_displ

    else:
        
        return 0.0, 0.0
    
def convert_to_node_topo(edof, ex, ey, ez, n_dofs_per_node=3, ignore_first=True):
    """
    Routine to convert dof based topology and element coordinates to node based
    topology required for visualisation with VTK and other visualisation frameworks

    Parameters
    ----------
    edof : array_like
        Element topology [nel x (n_dofs_per_node)|(n_dofs_per_node+1)*n_nodes ]
    ex : array_like
        Element x coordinates [nel x n_nodes]
    ey : array_like
        Element y coordinates [nel x n_nodes]
    ez : array_like
        Element z coordinates [nel x n_nodes]
    n_dofs_per_node : int, optional
        Number of dofs per node. Default is 3.
    ignore_first : bool, optional
        Ignore first column of edof. Default is True.
        
    Returns
    -------
    coords : numpy.ndarray
        Array of node coordinates. [n_nodes x 3]
    topo : numpy.ndarray
        Node topology. [nel x n_nodes]
    node_dofs : numpy.ndarray
        Dofs for each node. [n_nodes x n_dofs_per_node]
    """

    node_hash_coords = {}
    node_hash_numbers = {}
    node_hash_dofs = {}
    el_hash_dofs = []

    nel, cols = edof.shape

    if ignore_first:
        tot_dofs = cols-1
    else:
        tot_dofs = cols

    n_nodes = int(tot_dofs / n_dofs_per_node)

    for elx, ely, elz, dofs in zip(ex, ey, ez, edof):

        if ignore_first:
            el_dofs = dofs[1:]
        else:
            el_dofs = dofs

        # 0 1 2  3 4 5  6 7 8  9 12 11 

        el_dof = np.zeros((n_nodes, n_dofs_per_node), dtype=int)
        el_hash_topo = []

        for i in range(n_nodes):
            el_dof[i] = el_dofs[ (i*n_dofs_per_node):((i+1)*n_dofs_per_node) ]
            node_hash_coords[hash(tuple(el_dof[i]))] = [elx[i], ely[i], elz[i]]
            node_hash_numbers[hash(tuple(el_dof[i]))] = -1
            node_hash_dofs[hash(tuple(el_dof[i]))] = el_dof[i]
            el_hash_topo.append(hash(tuple(el_dof[i])))

        el_hash_dofs.append(el_hash_topo)

    coord_count = 0

    coords = []
    node_dofs = []

    for node_hash in node_hash_numbers.keys():
        node_hash_numbers[node_hash] = coord_count
        node_dofs.append(list(node_hash_dofs[node_hash]))
        coord_count +=1

        coords.append(node_hash_coords[node_hash])

    topo = []

    for el_hashes in el_hash_dofs:
        el_topo = []
        for i in range(n_nodes):
            el_topo.append(node_hash_numbers[el_hashes[i]])

        topo.append(el_topo)

    return np.array(coords), np.array(topo), np.array(node_dofs)
