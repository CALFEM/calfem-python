# -*- coding: utf-8 -*-
"""
Tests for the element flux drawing functions elflux2 and draw_element_flux.
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


def test_draw_element_flux_invalid_input():
    coords, edof, _ = quad_mesh()
    with pytest.raises(ValueError):
        cfv.draw_element_flux(np.zeros((2, 2)), coords, edof, 1, 3)
    with pytest.raises(ValueError):
        cfv.draw_element_flux(np.zeros(edof.shape[0]), coords, edof, 1, 3)
    coords3d = np.column_stack((coords, np.zeros(coords.shape[0])))
    with pytest.raises(ValueError):
        cfv.draw_element_flux(np.zeros((edof.shape[0], 2)), coords3d, edof, 1, 3)
