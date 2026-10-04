# -*- coding: utf-8 -*-
"""
CALFEM input/output module

Reading and writing of meshes, results, geometries and arrays.

* write_mesh / read_mesh: exchange meshes and results with other tools
  (ParaView, Gmsh, Abaqus, ...) in any format supported by meshio, e.g.
  VTU (.vtu), XDMF (.xdmf), VTK (.vtk) or Gmsh (.msh). Requires the optional
  dependency meshio (pip install calfem-python[io]).
* save_* / load_*: CALFEM's own pickle based formats and MATLAB export.

Usage::

    import calfem.io as cfio

    cfio.write_mesh("result.vtu", coords, edof, dofs_per_node, el_type,
                    point_data={"a": a}, cell_data={"von_mises": vm})
"""

import pickle

import numpy as np
import scipy.io

# Gmsh element type -> (meshio cell type, nodes per element). Only element
# types with identical node ordering in Gmsh and VTK/meshio are included.
_GMSH_TO_MESHIO = {
    1: ("line", 2),
    2: ("triangle", 3),
    3: ("quad", 4),
    4: ("tetra", 4),
    5: ("hexahedron", 8),
    9: ("triangle6", 6),
    16: ("quad8", 8),
}

_MESHIO_TO_GMSH = {cell_type: el_type
                   for el_type, (cell_type, _) in _GMSH_TO_MESHIO.items()}

_CELL_DIM = {"line": 1, "triangle": 2, "quad": 2, "triangle6": 2,
             "quad8": 2, "tetra": 3, "hexahedron": 3}


def _meshio():
    try:
        import meshio
    except ImportError as e:
        raise ImportError(
            "calfem.io.write_mesh/read_mesh require meshio. "
            "Install it with: pip install calfem-python[io]"
        ) from e
    return meshio


def _element_nodes(edof, dofs_per_node, dofs=None):
    """Zero-based node indices of each element from the dof topology."""
    edof = np.asarray(edof, dtype=int)
    first_dofs = edof[:, 0::dofs_per_node]

    if dofs is None:
        return (first_dofs - 1) // dofs_per_node

    dofs = np.asarray(dofs, dtype=int).reshape(-1, dofs_per_node)
    dof_to_node = np.full(dofs[:, 0].max() + 1, -1, dtype=int)
    dof_to_node[dofs[:, 0]] = np.arange(dofs.shape[0])
    nodes = dof_to_node[first_dofs]
    if (nodes < 0).any():
        raise ValueError("edof contains dofs that are not found in dofs.")
    return nodes


def _as_point_array(values, n_points):
    """Reshape a nodal value array to (n_points,) or (n_points, ncomp)."""
    values = np.asarray(values, dtype=float)
    if values.size % n_points != 0:
        raise ValueError(
            f"Point data of size {values.size} does not match "
            f"{n_points} nodes."
        )
    ncomp = values.size // n_points
    values = values.reshape(n_points, ncomp)
    return _pad_vectors(values)


def _as_cell_array(values, n_cells):
    values = np.asarray(values, dtype=float)
    if values.size % n_cells != 0:
        raise ValueError(
            f"Cell data of size {values.size} does not match "
            f"{n_cells} elements."
        )
    return _pad_vectors(values.reshape(n_cells, -1))


def _pad_vectors(values):
    """1 component -> 1D array, 2 components -> padded to 3D vectors so that
    e.g. ParaView's Warp By Vector and Glyph filters can be used."""
    if values.shape[1] == 1:
        return values[:, 0]
    if values.shape[1] == 2:
        return np.column_stack((values, np.zeros(values.shape[0])))
    return values


def write_mesh(filename, coords, edof, dofs_per_node, el_type,
               point_data=None, cell_data=None, dofs=None, file_format=None):
    """
    Write a mesh and optional results to file using meshio.

    The file format is determined from the file extension, e.g. .vtu or
    .vtk (ParaView), .xdmf, .msh (Gmsh) or .inp (Abaqus), or by file_format.
    XDMF additionally requires h5py.

    Parameters
    ----------
    filename : str or Path
        Name of the file to write.
    coords : array_like
        An N-by-2 or N-by-3 array. Row i contains the x,y,z coordinates of
        node i.
    edof : array_like
        An E-by-L array. Element topology. (E is the number of elements and
        L is the number of dofs per element)
    dofs_per_node : int
        Dofs per node.
    el_type : int
        Element type (Gmsh numbering): 1 line, 2 triangle, 3 quadrangle,
        4 tetrahedron, 5 hexahedron, 9 6-node triangle, 16 8-node
        quadrangle.
    point_data : dict, optional
        Nodal values, {name: values}. values can be a global vector, e.g.
        the displacement vector a with dofs_per_node values per node, or an
        N-by-k array. 2-component vectors are written as 3D vectors with a
        zero z-component.
    cell_data : dict, optional
        Element values, {name: values}. values is an array with one value
        or one row per element, e.g. effective stresses or flux vectors.
    dofs : array_like, optional
        An N-by-dofs_per_node array with the dofs of each node, as returned
        by the mesh generator. Only needed if the dofs are not numbered
        consecutively node by node (node i has dofs i*dofs_per_node+1, ...).
    file_format : str, optional
        meshio file format, e.g. "vtu" or "xdmf". Overrides the extension.
    """
    meshio = _meshio()

    if el_type not in _GMSH_TO_MESHIO:
        raise ValueError(
            f"Element type {el_type} not supported. Supported types: "
            f"{sorted(_GMSH_TO_MESHIO)}"
        )
    cell_type, nodes_per_el = _GMSH_TO_MESHIO[el_type]

    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[1] not in (2, 3):
        raise ValueError("coords must be an N-by-2 or N-by-3 array.")
    n_points = coords.shape[0]
    if coords.shape[1] == 2:
        points = np.column_stack((coords, np.zeros(n_points)))
    else:
        points = coords

    cells = _element_nodes(edof, dofs_per_node, dofs)
    if cells.shape[1] != nodes_per_el:
        raise ValueError(
            f"edof has {cells.shape[1]} nodes per element, element type "
            f"{el_type} ({cell_type}) requires {nodes_per_el}."
        )
    if cells.min() < 0 or cells.max() >= n_points:
        raise ValueError("edof refers to nodes outside coords.")
    n_cells = cells.shape[0]

    mesh_point_data = {
        name: _as_point_array(values, n_points)
        for name, values in (point_data or {}).items()
    }
    mesh_cell_data = {
        name: [_as_cell_array(values, n_cells)]
        for name, values in (cell_data or {}).items()
    }

    mesh = meshio.Mesh(
        points,
        [(cell_type, cells)],
        point_data=mesh_point_data,
        cell_data=mesh_cell_data,
    )
    meshio.write(str(filename), mesh, file_format=file_format)


def read_mesh(filename, dofs_per_node=1, file_format=None, return_data=False):
    """
    Read a mesh from file using meshio.

    Only the cells of the highest dimension are used, e.g. the triangles of
    a 2D mesh that also contains boundary lines. They must all be of one
    supported element type (see write_mesh).

    Parameters
    ----------
    filename : str or Path
        Name of the file to read.
    dofs_per_node : int, optional
        Dofs per node used to create edof and dofs. Default 1.
    file_format : str, optional
        meshio file format. Overrides the file extension.
    return_data : bool, optional
        Also return the point and cell data in the file. Default False.

    Returns
    -------
    coords : ndarray
        Node coordinates, N-by-2 if all z-coordinates are zero, otherwise
        N-by-3.
    edof : ndarray
        Element topology, dofs numbered node by node starting at 1.
    dofs : ndarray
        An N-by-dofs_per_node array with the dofs of each node.
    el_type : int
        Element type (Gmsh numbering).
    point_data : dict
        Nodal values, {name: array}. Only if return_data is True.
    cell_data : dict
        Element values for the returned elements, {name: array}. Only if
        return_data is True.
    """
    meshio = _meshio()
    mesh = meshio.read(str(filename), file_format=file_format)

    blocks = [(i, block) for i, block in enumerate(mesh.cells)
              if block.type in _CELL_DIM]
    if not blocks:
        raise ValueError(
            "No supported elements found. Supported cell types: "
            f"{sorted(_MESHIO_TO_GMSH)}"
        )
    max_dim = max(_CELL_DIM[block.type] for _, block in blocks)
    blocks = [(i, block) for i, block in blocks
              if _CELL_DIM[block.type] == max_dim]
    cell_types = {block.type for _, block in blocks}
    if len(cell_types) > 1:
        raise ValueError(
            f"Mixed element types are not supported: {sorted(cell_types)}"
        )

    cell_type = cell_types.pop()
    el_type = _MESHIO_TO_GMSH[cell_type]
    cells = np.vstack([block.data for _, block in blocks])

    points = np.asarray(mesh.points, dtype=float)
    if points.shape[1] == 3 and np.allclose(points[:, 2], 0.0):
        coords = points[:, :2].copy()
    else:
        coords = points.copy()

    n_points = coords.shape[0]
    dofs = np.arange(1, n_points*dofs_per_node + 1).reshape(n_points,
                                                             dofs_per_node)
    edof = dofs[cells].reshape(cells.shape[0], -1)

    if not return_data:
        return coords, edof, dofs, el_type

    point_data = {name: np.asarray(values)
                  for name, values in mesh.point_data.items()}
    cell_data = {}
    for name, block_values in mesh.cell_data.items():
        cell_data[name] = np.concatenate(
            [np.asarray(block_values[i]) for i, _ in blocks]
        )
    return coords, edof, dofs, el_type, point_data, cell_data


# ------------------------------------------------- CALFEM native formats

def save_geometry(g, name="Untitled"):
    """Save a geometry object to a .cfg file (pickle)."""
    if not name.endswith(".cfg"):
        name = name + ".cfg"
    with open(name, "wb") as file:
        pickle.dump(g, file)


def load_geometry(name):
    """Load a geometry object saved with save_geometry."""
    with open(name, "rb") as file:
        return pickle.load(file)


def save_mesh(mesh, name="Untitled"):
    """Save a mesh object to a .cfm file (pickle)."""
    if not name.endswith(".cfm"):
        name = name + ".cfm"
    with open(name, "wb") as file:
        pickle.dump(mesh, file)


def load_mesh(name):
    """Load a mesh object saved with save_mesh."""
    with open(name, "rb") as file:
        return pickle.load(file)


def save_arrays(coords, edof, dofs, bdofs, elementmarkers, boundary_elements,
                marker_dict, name="Untitled"):
    """Save the arrays returned by the mesh generator to a .cfma file."""
    if not name.endswith(".cfma"):
        name = name + ".cfma"
    with open(name, "wb") as file:
        pickle.dump(coords, file)
        pickle.dump(edof, file)
        pickle.dump(dofs, file)
        pickle.dump(bdofs, file)
        pickle.dump(elementmarkers, file)
        pickle.dump(boundary_elements, file)
        pickle.dump(marker_dict, file)


def load_arrays(name):
    """
    Load arrays saved with save_arrays.

    Returns
    -------
    coords, edof, dofs, bdofs, elementmarkers, boundary_elements, marker_dict
    """
    with open(name, "rb") as file:
        coords = pickle.load(file)
        edof = pickle.load(file)
        dofs = pickle.load(file)
        bdofs = pickle.load(file)
        elementmarkers = pickle.load(file)
        boundary_elements = pickle.load(file)
        marker_dict = pickle.load(file)

    return coords, edof, dofs, bdofs, elementmarkers, boundary_elements, marker_dict


def _matlab_field_name(name):
    """Make a valid MATLAB struct field name."""
    field = "".join(c if c.isalnum() or c == "_" else "_" for c in str(name))
    if not field or not field[0].isalpha():
        field = "marker_" + field
    return field[:63]


def save_matlab_arrays(coords, edof, dofs, bdofs, elementmarkers,
                       boundary_elements=None, marker_dict=None,
                       name="Untitled"):
    """
    Save mesh arrays to a MATLAB .mat file.

    edof is written in the MATLAB CALFEM format with the element number as
    first column. bdofs is written as a struct with one field per marker,
    named after marker_dict[marker] if available, otherwise marker_<marker>.
    elementmarkers are increased by one to suit one-based indexing.
    boundary_elements is accepted for compatibility but not written.
    """
    if not name.endswith(".mat"):
        name = name + ".mat"

    edof = np.asarray(edof)
    element_numbers = np.arange(1, edof.shape[0] + 1)[:, np.newaxis]

    marker_dict = marker_dict or {}
    matlab_bdofs = {}
    for marker, marker_dofs in bdofs.items():
        label = marker_dict.get(marker, marker)
        if isinstance(label, (int, np.integer)):
            label = f"marker_{label}"
        matlab_bdofs[_matlab_field_name(label)] = np.asarray(
            marker_dofs, dtype=float)

    save_dict = {
        "coords": np.asarray(coords, dtype=float),
        "edof": np.hstack((element_numbers, edof)).astype(float),
        "dofs": np.asarray(dofs, dtype=float),
        "bdofs": matlab_bdofs,
        "elementmarkers": np.asarray(elementmarkers) + 1,
    }
    scipy.io.savemat(name, save_dict)


# camelCase aliases for compatibility with calfem._export
saveGeometry = save_geometry
loadGeometry = load_geometry
saveMesh = save_mesh
loadMesh = load_mesh
saveArrays = save_arrays
loadArrays = load_arrays
saveMatlabArrays = save_matlab_arrays
