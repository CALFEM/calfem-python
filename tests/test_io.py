# -*- coding: utf-8 -*-
"""
Tests for the calfem.io module.
"""

import sys
from pathlib import Path

# Add src directory to system path to use local calfem package
script_dir = Path(__file__).parent.parent
src_dir = script_dir / "src"
sys.path.insert(0, str(src_dir))

import numpy as np
import pytest
import scipy.io

import calfem.io as cfio


def quad_mesh(nx=3, ny=2, dofs_per_node=1):
    """Structured quad mesh with dofs numbered node by node."""
    xs, ys = np.meshgrid(np.linspace(0, 3, nx + 1), np.linspace(0, 2, ny + 1))
    coords = np.column_stack((xs.ravel(), ys.ravel()))
    nodes = []
    for j in range(ny):
        for i in range(nx):
            n = j*(nx + 1) + i
            nodes.append([n, n + 1, n + nx + 2, n + nx + 1])
    nodes = np.array(nodes)
    dofs = np.arange(1, coords.shape[0]*dofs_per_node + 1).reshape(
        -1, dofs_per_node)
    edof = dofs[nodes].reshape(nodes.shape[0], -1)
    return coords, edof, dofs, nodes


def hex_mesh():
    coords = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                       [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
                       [2, 0, 0], [2, 1, 0], [2, 0, 1], [2, 1, 1]], float)
    nodes = np.array([[0, 1, 2, 3, 4, 5, 6, 7],
                      [1, 8, 9, 2, 5, 10, 11, 6]])
    dofs = np.arange(1, 3*coords.shape[0] + 1).reshape(-1, 3)
    edof = dofs[nodes].reshape(2, -1)
    return coords, edof, dofs, nodes


# ------------------------------------------------------- write/read mesh

@pytest.fixture
def meshio():
    return pytest.importorskip("meshio")


@pytest.mark.parametrize("ext", ["vtu", "vtk", "xdmf", "msh"])
def test_write_read_roundtrip_2d(meshio, tmp_path, ext):
    if ext == "xdmf":
        pytest.importorskip("h5py")
    coords, edof, dofs, _ = quad_mesh(dofs_per_node=2)
    filename = tmp_path / f"mesh.{ext}"
    cfio.write_mesh(filename, coords, edof, 2, 3)

    coords_r, edof_r, dofs_r, el_type = cfio.read_mesh(filename,
                                                       dofs_per_node=2)
    assert el_type == 3
    assert np.allclose(coords_r, coords)
    assert np.array_equal(edof_r, edof)
    assert np.array_equal(dofs_r, dofs)


def test_write_read_point_and_cell_data(meshio, tmp_path):
    coords, edof, dofs, nodes = quad_mesh(dofs_per_node=2)
    a = np.random.default_rng(0).standard_normal(dofs.size)
    temperature = np.arange(coords.shape[0], dtype=float)
    vm = np.arange(edof.shape[0], dtype=float)
    flux = np.random.default_rng(1).standard_normal((edof.shape[0], 2))
    stress = np.random.default_rng(2).standard_normal((edof.shape[0], 3))

    filename = tmp_path / "result.vtu"
    cfio.write_mesh(filename, coords, edof, 2, 3,
                    point_data={"a": a, "T": temperature},
                    cell_data={"von_mises": vm, "flux": flux,
                               "stress": stress})

    *_, point_data, cell_data = cfio.read_mesh(filename, dofs_per_node=2,
                                               return_data=True)
    # 2-component vectors are padded to 3D
    assert point_data["a"].shape == (coords.shape[0], 3)
    assert np.allclose(point_data["a"][:, :2], a.reshape(-1, 2))
    assert np.allclose(point_data["a"][:, 2], 0)
    assert np.allclose(point_data["T"], temperature)
    assert np.allclose(cell_data["von_mises"], vm)
    assert np.allclose(cell_data["flux"][:, :2], flux)
    assert np.allclose(cell_data["stress"], stress)


def test_write_mesh_3d_hexahedra(meshio, tmp_path):
    coords, edof, dofs, nodes = hex_mesh()
    filename = tmp_path / "hex.vtu"
    cfio.write_mesh(filename, coords, edof, 3, 5,
                    point_data={"u": np.ones(dofs.size)})

    mesh = meshio.read(filename)
    assert mesh.cells[0].type == "hexahedron"
    assert np.array_equal(mesh.cells[0].data, nodes)
    assert np.allclose(mesh.points, coords)
    assert mesh.point_data["u"].shape == (coords.shape[0], 3)

    coords_r, edof_r, dofs_r, el_type = cfio.read_mesh(filename,
                                                       dofs_per_node=3)
    assert el_type == 5
    assert coords_r.shape == (12, 3)
    assert np.array_equal(edof_r, edof)


def test_write_mesh_triangles_with_dofs_array(meshio, tmp_path):
    """Non-consecutive dof numbering is mapped back to nodes via dofs."""
    coords, _, _, quad_nodes = quad_mesh()
    nodes = np.vstack([quad_nodes[:, [0, 1, 2]], quad_nodes[:, [0, 2, 3]]])
    dofs = (np.arange(coords.shape[0])[::-1] + 1).reshape(-1, 1)*10
    edof = dofs[nodes].reshape(nodes.shape[0], -1)

    filename = tmp_path / "tri.vtu"
    cfio.write_mesh(filename, coords, edof, 1, 2, dofs=dofs)
    mesh = meshio.read(filename)
    assert mesh.cells[0].type == "triangle"
    assert np.array_equal(mesh.cells[0].data, nodes)


def test_read_mesh_uses_highest_dimension(meshio, tmp_path):
    coords, _, _, nodes = quad_mesh()
    points = np.column_stack((coords, np.zeros(coords.shape[0])))
    lines = np.array([[0, 1], [1, 2], [2, 3]])
    mesh = meshio.Mesh(points, [("line", lines), ("quad", nodes)])
    filename = tmp_path / "mixed.vtu"
    meshio.write(filename, mesh)

    coords_r, edof, _, el_type = cfio.read_mesh(filename)
    assert el_type == 3
    assert np.array_equal(edof, nodes + 1)


def test_write_mesh_invalid_input(meshio, tmp_path):
    coords, edof, _, _ = quad_mesh()
    filename = tmp_path / "bad.vtu"
    with pytest.raises(ValueError):
        cfio.write_mesh(filename, coords, edof, 1, 99)
    with pytest.raises(ValueError):
        cfio.write_mesh(filename, coords, edof, 1, 2)  # quads as triangles
    with pytest.raises(ValueError):
        cfio.write_mesh(filename, coords, edof, 1, 3,
                        point_data={"x": np.ones(5)})
    with pytest.raises(ValueError):
        cfio.write_mesh(filename, coords, edof, 1, 3,
                        cell_data={"x": np.ones(5)})
    with pytest.raises(ValueError):
        cfio.write_mesh(filename, coords[:4], edof, 1, 3)


def test_missing_meshio_gives_helpful_error(monkeypatch, tmp_path):
    monkeypatch.setitem(sys.modules, "meshio", None)
    coords, edof, _, _ = quad_mesh()
    with pytest.raises(ImportError, match="calfem-python\\[io\\]"):
        cfio.write_mesh(tmp_path / "m.vtu", coords, edof, 1, 3)


# ------------------------------------------------- CALFEM native formats

def test_save_load_arrays(tmp_path):
    coords, edof, dofs, _ = quad_mesh()
    bdofs = {10: [1, 2, 3], 20: [7, 8]}
    markers = [0]*edof.shape[0]
    name = str(tmp_path / "arrays")
    cfio.save_arrays(coords, edof, dofs, bdofs, markers, {}, {10: "left"},
                     name=name)
    out = cfio.load_arrays(name + ".cfma")
    assert np.array_equal(out[0], coords)
    assert np.array_equal(out[1], edof)
    assert out[3] == bdofs
    assert out[6] == {10: "left"}


def test_save_load_mesh_and_geometry(tmp_path):
    obj = {"points": [[0, 0], [1, 0]]}
    cfio.save_mesh(obj, str(tmp_path / "m"))
    assert cfio.load_mesh(str(tmp_path / "m.cfm")) == obj
    cfio.save_geometry(obj, str(tmp_path / "g"))
    assert cfio.load_geometry(str(tmp_path / "g.cfg")) == obj


def test_save_matlab_arrays(tmp_path):
    coords, edof, dofs, _ = quad_mesh()
    bdofs = {10: [1, 2, 3], 20: [7, 8]}
    name = str(tmp_path / "mesh")
    cfio.save_matlab_arrays(coords, edof, dofs, bdofs, [0]*edof.shape[0],
                            {}, {10: "left side"}, name=name)

    data = scipy.io.loadmat(name + ".mat")
    assert np.array_equal(data["edof"][:, 0], np.arange(1, edof.shape[0] + 1))
    assert np.array_equal(data["edof"][:, 1:], edof)
    assert np.allclose(data["coords"], coords)
    assert np.all(data["elementmarkers"] == 1)
    bd = data["bdofs"][0, 0]
    assert set(bd.dtype.names) == {"left_side", "marker_20"}
    assert np.array_equal(bd["left_side"].ravel(), [1, 2, 3])


def test_export_module_backwards_compatible(tmp_path):
    import calfem._export as cfe
    cfe.saveMesh({"a": 1}, str(tmp_path / "m"))
    assert cfe.loadMesh(str(tmp_path / "m.cfm")) == {"a": 1}
    assert cfe.saveMatlabArrays is cfio.save_matlab_arrays
