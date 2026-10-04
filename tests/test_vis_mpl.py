# -*- coding: utf-8 -*-
"""
Tests for calfem.vis_mpl: element flux drawing (elflux2, draw_element_flux)
and displaced meshes (draw_element_values, draw_displacements).
"""

import sys
from pathlib import Path

# Add src directory to system path to use local calfem package
script_dir = Path(__file__).parent.parent
src_dir = script_dir / "src"
sys.path.insert(0, str(src_dir))

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest
import matplotlib.pyplot as plt

import calfem.core as cfc
import calfem.vis_mpl as cfv


def quad_mesh(nx=3, ny=2):
    """Structured quad mesh, one dof per node (dof = node number)."""
    xs, ys = np.meshgrid(np.linspace(0, 3, nx + 1), np.linspace(0, 2, ny + 1))
    coords = np.column_stack((xs.ravel(), ys.ravel()))
    edof = []
    for j in range(ny):
        for i in range(nx):
            n = j*(nx + 1) + i
            edof.append([n, n + 1, n + nx + 2, n + nx + 1])
    edof = np.array(edof) + 1
    dofs = np.arange(1, coords.shape[0] + 1).reshape(-1, 1)
    return coords, edof, dofs


def tri_mesh():
    coords, quad_edof, dofs = quad_mesh()
    edof = np.vstack([quad_edof[:, [0, 1, 2]], quad_edof[:, [0, 2, 3]]])
    return coords, edof, dofs


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


def quiver_data(ax):
    q = [c for c in ax.collections
         if isinstance(c, matplotlib.quiver.Quiver)][-1]
    return q.X, q.Y, q.U, q.V, q


@pytest.mark.parametrize("mesh, el_type", [(quad_mesh, 3), (tri_mesh, 2)])
def test_draw_element_flux_matches_elflux2(mesh, el_type):
    coords, edof, dofs = mesh()
    ex, ey = cfc.coordxtr(edof, coords, dofs)
    flux = np.random.default_rng(0).standard_normal((edof.shape[0], 2))

    _, ax1 = plt.subplots()
    s1 = cfv.elflux2(ex, ey, flux, ax=ax1)
    _, ax2 = plt.subplots()
    s2 = cfv.draw_element_flux(flux, coords, edof, 1, el_type, ax=ax2)

    assert np.isclose(s1, s2)
    for a, b in zip(quiver_data(ax1)[:4], quiver_data(ax2)[:4]):
        assert np.allclose(a, b)


def test_draw_element_flux_centroids_and_scale():
    coords, edof, _ = quad_mesh()
    flux = np.tile([2.0, -1.0], (edof.shape[0], 1))
    _, ax = plt.subplots()
    scale = cfv.draw_element_flux(flux, coords, edof, 1, 3, scale=0.25, ax=ax)
    X, Y, U, V, _ = quiver_data(ax)

    assert scale == 0.25
    assert np.allclose(np.sort(np.unique(X)), [0.5, 1.5, 2.5])
    assert np.allclose(np.sort(np.unique(Y)), [0.5, 1.5])
    assert np.allclose(U, 0.5) and np.allclose(V, -0.25)


def test_draw_element_flux_automatic_scale():
    """Mean arrow length = 0.8 * mean element diagonal."""
    coords, edof, _ = quad_mesh()
    flux = np.tile([3.0, 4.0], (edof.shape[0], 1))
    scale = cfv.draw_element_flux(flux, coords, edof, 1, 3)
    assert np.isclose(scale*5.0, 0.8*np.sqrt(2))


def test_draw_element_flux_zero_flux():
    coords, edof, _ = quad_mesh()
    flux = np.zeros((edof.shape[0], 2))
    assert cfv.draw_element_flux(flux, coords, edof, 1, 3) == 0.0


def test_draw_element_flux_color_by_magnitude_sets_mappable():
    coords, edof, _ = quad_mesh()
    flux = np.random.default_rng(1).standard_normal((edof.shape[0], 2))
    _, ax = plt.subplots()
    cfv.draw_element_flux(flux, coords, edof, 1, 3, color_by_magnitude=True,
                          cmap="viridis", ax=ax)
    q = quiver_data(ax)[4]
    assert np.allclose(q.get_array(), np.linalg.norm(flux, axis=1))
    cfv.colorbar()


def test_draw_element_flux_draw_elements_and_limits():
    coords, edof, _ = quad_mesh()
    flux = np.tile([1.0, 0.0], (edof.shape[0], 1))
    _, ax = plt.subplots()
    cfv.draw_element_flux(flux, coords, edof, 1, 3, draw_elements=True,
                          title="Flux", ax=ax)
    polys = [c for c in ax.collections
             if isinstance(c, matplotlib.collections.PolyCollection)
             and not isinstance(c, matplotlib.quiver.Quiver)]
    assert len(polys) == 1
    assert ax.get_title() == "Flux"
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    assert x0 <= 0 and x1 >= 3 and y0 <= 0 and y1 >= 2


def test_draw_element_flux_multiple_dofs_per_node():
    """edof with 2 dofs per node (e.g. a structural mesh) is handled."""
    coords, edof, _ = quad_mesh()
    edof2 = np.column_stack([np.column_stack((2*edof[:, k] - 1, 2*edof[:, k]))
                             for k in range(4)])
    flux = np.random.default_rng(2).standard_normal((edof.shape[0], 2))
    _, ax1 = plt.subplots()
    cfv.draw_element_flux(flux, coords, edof, 1, 3, ax=ax1)
    _, ax2 = plt.subplots()
    cfv.draw_element_flux(flux, coords, edof2, 2, 3, ax=ax2)
    for a, b in zip(quiver_data(ax1)[:4], quiver_data(ax2)[:4]):
        assert np.allclose(a, b)


def _quad_mesh_2dof():
    coords, edof, _ = quad_mesh()
    edof2 = np.column_stack([np.column_stack((2*edof[:, k] - 1, 2*edof[:, k]))
                             for k in range(4)])
    return coords, edof2


def _polygon_vertices(ax):
    pc = [c for c in ax.collections
          if isinstance(c, matplotlib.collections.PolyCollection)][-1]
    return np.vstack([p.vertices[:4] for p in pc.get_paths()])


@pytest.mark.parametrize("as_array", [False, True])
def test_draw_element_values_displacements(as_array):
    """Displacements are applied both as global vector and N-by-2 array."""
    coords, edof2 = _quad_mesh_2dof()
    a = np.tile([0.5, -0.25], coords.shape[0])
    if as_array:
        a = a.reshape(-1, 2)
    _, ax = plt.subplots()
    cfv.draw_element_values(np.ones(edof2.shape[0]), coords, edof2, 2, 3,
                            displacements=a, magnfac=2.0)
    verts = _polygon_vertices(ax)
    assert np.isclose(verts[:, 0].min(), 1.0)
    assert np.isclose(verts[:, 1].min(), -0.5)


@pytest.mark.parametrize("as_array", [False, True])
def test_draw_displacements_auto_magnfac(as_array):
    """Largest displacement becomes magscale (0.1) times the model size."""
    coords, edof2 = _quad_mesh_2dof()
    a = np.zeros(2*coords.shape[0])
    a[0] = 1e-6                     # node 0 (at the origin) moves in x
    if as_array:
        a = a.reshape(-1, 2)
    _, ax = plt.subplots()
    cfv.draw_displacements(a, coords, edof2, 2, 3)
    verts = _polygon_vertices(ax)
    assert np.isclose(verts[0, 0], 0.1*3.0)


def test_draw_displacements_given_magnfac():
    coords, edof2 = _quad_mesh_2dof()
    a = np.tile([0.0, 0.1], coords.shape[0])
    _, ax = plt.subplots()
    cfv.draw_displacements(a, coords, edof2, 2, 3, magnfac=10.0)
    verts = _polygon_vertices(ax)
    assert np.isclose(verts[:, 1].min(), 1.0)


# ------------------------------------------------------- classic functions

def element_arrays():
    coords, edof, _ = quad_mesh()
    nodes = edof - 1
    return coords[nodes, 0], coords[nodes, 1], nodes, coords


def _all_polygons(ax):
    pc = [c for c in ax.collections
          if isinstance(c, matplotlib.collections.PolyCollection)][0]
    return [p.vertices for p in pc.get_paths()]


def test_eldraw2_element_numbers():
    ex, ey, nodes, _ = element_arrays()
    _, ax = plt.subplots()
    cfv.eldraw2(ex, ey, [1, 2, 1], elnum=range(1, ex.shape[0] + 1))
    texts = [t.get_text() for t in ax.texts]
    assert texts == [str(i) for i in range(1, ex.shape[0] + 1)]
    assert np.isclose(ax.texts[0].get_position()[0], ex[0].mean())


def test_eldraw2_lists_and_wrong_elnum():
    _, ax = plt.subplots()
    cfv.eldraw2([0.0, 1.0], [0.0, 0.0], [1, 1, 0])
    with pytest.raises(ValueError):
        cfv.eldraw2(np.zeros((2, 3)), np.zeros((2, 3)), elnum=[1])


def test_eldraw2_8_node_outline_order():
    ex = np.array([[0, 2, 2, 0, 1, 2, 1, 0]], float)
    ey = np.array([[0, 0, 2, 2, 0, 1, 2, 1]], float)
    _, ax = plt.subplots()
    cfv.eldraw2(ex, ey, [1, 1, 0])
    verts = _all_polygons(ax)[0]
    assert np.allclose(verts[:8, 0], [0, 1, 2, 2, 2, 1, 0, 0])


def test_eldisp2_quads_returns_sfac():
    ex, ey, nodes, coords = element_arrays()
    a = np.column_stack((0.01*coords[:, 0], np.zeros(coords.shape[0])))
    ed = a[nodes].reshape(nodes.shape[0], -1)
    _, ax = plt.subplots()
    sfac = cfv.eldisp2(ex, ey, ed, [1, 2, 1], sfac=10.0)
    assert sfac == 10.0
    polys = _all_polygons(ax)
    assert len(polys) == ex.shape[0]
    assert np.allclose(polys[0][:4, 0], 1.1*ex[0])


def test_eldisp2_auto_scale_and_triangles():
    ex = np.array([[0.0, 2.0, 0.0]])
    ey = np.array([[0.0, 0.0, 1.0]])
    ed = np.array([[0, 0, 4e-3, 0, 0, 0]])
    _, ax = plt.subplots()
    sfac = cfv.eldisp2(ex, ey, ed)
    assert np.isclose(sfac, 0.1*2.0/4e-3)
    assert np.allclose(_all_polygons(ax)[0][1], [2.2, 0.0])


def test_eldisp2_beams():
    ex = np.array([[0.0, 1.0]])
    ey = np.array([[0.0, 0.0]])
    ed = np.array([[0, 0, 0, 0, 0.1, 0]])
    _, ax = plt.subplots()
    cfv.eldisp2(ex, ey, ed, [1, 1, 0], sfac=2.0)
    curve = _all_polygons(ax)[0]
    assert len(curve) > 2
    assert np.allclose(curve[0], [0, 0]) and np.allclose(curve[-1], [1, 0.2])


def test_eldisp2_8_node_edges_pass_through_nodes():
    ex = np.array([[0, 2, 2, 0, 1, 2, 1, 0]], float)
    ey = np.array([[0, 0, 2, 2, 0, 1, 2, 1]], float)
    ed = np.zeros((1, 16))
    ed[0, 9] = 0.5                       # midside node 5 moves in y
    _, ax = plt.subplots()
    cfv.eldisp2(ex, ey, ed, [1, 1, 0], sfac=1.0)
    pts = _all_polygons(ax)[0]
    for xn, yn in [(1, 0.5), (0, 0), (2, 2)]:
        assert np.min(np.hypot(pts[:, 0] - xn, pts[:, 1] - yn)) < 1e-12


def test_eldisp2_unsupported():
    with pytest.raises(ValueError):
        cfv.eldisp2(np.zeros((1, 5)), np.zeros((1, 5)), np.ones((1, 10)))


def test_dispbeam2_and_secforce2_return_sfac():
    edi = np.column_stack((np.zeros(5), [0, -1, -2, -1, 0]))
    assert np.isclose(cfv.dispbeam2([0, 4], [0, 0], edi), 0.1*4/2)
    assert cfv.dispbeam2([0, 4], [0, 0], edi, sfac=0.3) == 0.3
    assert np.isclose(cfv.secforce2([0, 2], [0, 0], [2.0, 1.0, -1.0]),
                      0.2*2/2.0)


def test_eliso2():
    ex, ey, nodes, coords = element_arrays()
    ed = 2*ex + ey                       # linear field
    _, ax = plt.subplots()
    cs = cfv.eliso2(ex, ey, ed, isov=[1, 2, 3, 4])
    assert list(cs.levels) == [1, 2, 3, 4]
    for level, segs in zip(cs.levels, cs.allsegs):
        assert segs, f"no isoline at {level}"
        for seg in segs:
            assert np.allclose(2*seg[:, 0] + seg[:, 1], level)
    cfv.colorbar()


def test_eliso2_single_color_and_alias():
    ex, ey, nodes, coords = element_arrays()
    _, ax = plt.subplots()
    cs = cfv.eliso2(ex, ey, ex + ey, isov=4, plotpar=[2, 4])
    assert np.allclose(cs.get_edgecolor()[0][:3], [1, 0, 0])
    assert cfv.eliso2_mpl is cfv.eliso2


def test_draw_element_flux_invalid_input():
    coords, edof, _ = quad_mesh()
    with pytest.raises(ValueError):
        cfv.draw_element_flux(np.zeros((2, 2)), coords, edof, 1, 3)
    with pytest.raises(ValueError):
        cfv.draw_element_flux(np.zeros(edof.shape[0]), coords, edof, 1, 3)
    coords3d = np.column_stack((coords, np.zeros(coords.shape[0])))
    with pytest.raises(ValueError):
        cfv.draw_element_flux(np.zeros((edof.shape[0], 2)), coords3d, edof, 1, 3)
