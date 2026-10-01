# -*- coding: utf-8 -*-
"""
Tests for the plotly visualisation module calfem.vis_plotly.
"""

import sys
from pathlib import Path

# Add src directory to system path to use local calfem package
script_dir = Path(__file__).parent.parent
src_dir = script_dir / "src"
sys.path.insert(0, str(src_dir))

import numpy as np
import pytest

go = pytest.importorskip("plotly.graph_objects")

import calfem.geometry as cfg
import calfem.vis_plotly as cfv


def quad_mesh(nx=3, ny=2, dofs_per_node=1):
    """Structured quad mesh on [0, 3] x [0, 2], dofs numbered node-wise."""
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
    return coords, edof, nodes


def hex_mesh():
    coords = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                       [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
                       [2, 0, 0], [2, 1, 0], [2, 0, 1], [2, 1, 1]], float)
    nodes = np.array([[0, 1, 2, 3, 4, 5, 6, 7],
                      [1, 8, 9, 2, 5, 10, 11, 6]])
    dofs = np.arange(1, 3*coords.shape[0] + 1).reshape(-1, 3)
    return coords, dofs[nodes].reshape(2, -1), nodes


def n_polylines(trace):
    """Number of NaN/None separated polylines in a trace."""
    x = np.array([np.nan if v is None else v for v in trace.x], float)
    return int(np.isnan(x).sum())


@pytest.fixture(autouse=True)
def fresh_figures():
    cfv.close_all()
    yield
    cfv.close_all()


# ------------------------------------------------------- figure handling

def test_figure_and_show(monkeypatch):
    shown = []
    monkeypatch.setattr(cfv, "_show_in_browser",
                        lambda figs, *a, **k: shown.extend(figs))
    f1 = cfv.figure()
    assert cfv.gcf() is f1
    f2 = cfv.figure(fig_size=(4, 3), title="T")
    assert cfv.gcf() is f2
    assert f2.layout.width == 400 and f2.layout.height == 300
    assert f2.layout.title.text == "T"
    cfv.figure(show=False)
    cfv.show()
    assert shown == [f1, f2]
    cfv.show()
    assert len(shown) == 2


def test_show_writes_html_files_in_scripts(monkeypatch, tmp_path):
    opened = []
    monkeypatch.setattr(cfv, "_open_file", lambda f: opened.append(f))
    monkeypatch.delenv("PLOTLY_RENDERER", raising=False)
    coords, edof, _ = quad_mesh()
    cfv.figure()
    cfv.draw_mesh(coords, edof, 1, 3)
    cfv.figure()
    cfv.draw_mesh(coords, edof, 1, 3)
    files = cfv.show(directory=tmp_path)
    assert [f.name for f in files] == ["figure_1.html", "figure_2.html"]
    assert all(f.stat().st_size > 0 for f in files)
    assert (tmp_path / "plotly.min.js").exists()
    assert opened == files
    assert cfv.show() == []                    # nothing left to show


def test_show_reports_when_browser_cannot_open(monkeypatch, tmp_path,
                                               capsys):
    def fail(filename):
        raise OSError("no application associated")
    monkeypatch.setattr(cfv, "_open_file", fail)
    monkeypatch.delenv("PLOTLY_RENDERER", raising=False)
    cfv.figure()
    files = cfv.show(directory=tmp_path)
    assert files[0].exists()
    out = capsys.readouterr().out
    assert str(tmp_path) in out
    assert "could not be opened automatically" in out


def test_open_file_ignores_broken_browser_variable(monkeypatch, tmp_path):
    """On Windows/macOS the OS opens the file, BROWSER is not used."""
    calls = []
    monkeypatch.setenv("BROWSER", str(tmp_path / "missing" / "browser.sh"))
    if sys.platform == "win32":
        monkeypatch.setattr("os.startfile", lambda f: calls.append(f),
                            raising=False)
    elif sys.platform == "darwin":
        monkeypatch.setattr("subprocess.run",
                            lambda args, **k: calls.append(args[1]))
    else:
        pytest.skip("Linux uses webbrowser, which follows BROWSER")
    filename = tmp_path / "figure_1.html"
    cfv._open_file(filename)
    assert calls == [str(filename)]


def test_show_uses_plotly_renderer_when_requested(monkeypatch):
    shown = []
    monkeypatch.setattr(go.Figure, "show",
                        lambda self, *a, **k: shown.append(k.get("renderer")))
    monkeypatch.setattr(cfv, "_show_in_browser",
                        lambda *a, **k: pytest.fail("wrote HTML files"))
    cfv.figure()
    cfv.show(renderer="json")
    assert shown == ["json"]

    monkeypatch.setenv("PLOTLY_RENDERER", "json")
    cfv.figure()
    cfv.show()
    assert shown == ["json", None]


def test_show_inline_in_notebooks(monkeypatch):
    shown = []
    monkeypatch.setattr(go.Figure, "show",
                        lambda self, *a, **k: shown.append(self))
    monkeypatch.setattr(cfv, "_in_notebook", lambda: True)
    monkeypatch.delenv("PLOTLY_RENDERER", raising=False)
    fig = cfv.figure()
    cfv.show()
    assert len(shown) == 1 and shown[0] is fig


def test_gcf_creates_figure_and_close():
    fig = cfv.gcf()
    assert isinstance(fig, go.Figure)
    cfv.close()
    assert cfv.gcf() is not fig


def test_color_conversion():
    assert cfv._color((0, 0, 0)) == "rgb(0,0,0)"
    assert cfv._color((1.0, 0.5, 0.0)) == "rgb(255,128,0)"
    assert cfv._color((0.2, 0.4, 0.6, 0.5)) == "rgba(51,102,153,0.5)"
    assert cfv._color("r") == "red"
    assert cfv._color("#123456") == "#123456"


def test_save_figure_html(tmp_path):
    coords, edof, _ = quad_mesh()
    cfv.draw_mesh(coords, edof, 1, 3)
    filename = tmp_path / "mesh.html"
    cfv.save_figure(filename)
    assert filename.stat().st_size > 0


def test_title_axis_text():
    cfv.title("Hello")
    cfv.axis("off")
    cfv.text("A", [1.0, 2.0])
    fig = cfv.gcf()
    assert fig.layout.title.text == "Hello"
    assert fig.layout.xaxis.visible is False
    assert fig.layout.annotations[0].text == "A"


# ------------------------------------------------------- mesh drawing

def test_draw_mesh_2d():
    coords, edof, nodes = quad_mesh(dofs_per_node=2)
    fig = cfv.draw_mesh(coords, edof, 2, 3, title="Mesh", show_nodes=True)
    lines, node_markers = fig.data
    assert n_polylines(lines) == nodes.shape[0]
    # First polyline is the closed first element
    xy = np.column_stack((lines.x[:5], lines.y[:5]))
    assert np.allclose(xy, coords[np.r_[nodes[0], nodes[0][0]]])
    assert len(node_markers.x) == coords.shape[0]
    assert fig.layout.yaxis.scaleanchor == "x"
    assert fig.layout.title.text == "Mesh"


def test_draw_mesh_3d_hexahedra():
    coords, edof, nodes = hex_mesh()
    fig = cfv.draw_mesh(coords, edof, 3, 5)
    mesh3d, edges = fig.data
    assert isinstance(mesh3d, go.Mesh3d)
    assert len(mesh3d.i) == 2*6*2          # 2 elements, 6 faces, 2 triangles
    assert isinstance(edges, go.Scatter3d)
    assert n_polylines(edges) == 20        # unique edges of two hexahedra
    assert fig.layout.scene.aspectmode == "data"


def test_draw_elements():
    coords, edof, nodes = quad_mesh()
    ex = coords[nodes, 0]
    ey = coords[nodes, 1]
    fig = cfv.draw_elements(ex, ey, line_style="dashed", filled=True)
    assert n_polylines(fig.data[0]) == nodes.shape[0]
    assert fig.data[0].line.dash == "dash"
    assert fig.data[0].fill == "toself"


# ------------------------------------------------------- results

def test_draw_element_values_2d():
    coords, edof, nodes = quad_mesh(nx=6, ny=4)
    values = np.arange(nodes.shape[0], dtype=float)
    fig = cfv.draw_element_values(values, coords, edof, 1, 3,
                                  colorbar_title="v", n_colors=8)
    *fills, edges, hover = fig.data

    # Each element is drawn in exactly one color band
    assert all(f.fill == "toself" for f in fills)
    assert sum(n_polylines(f) for f in fills) == nodes.shape[0]
    assert len(fills) == 8
    assert n_polylines(edges) == nodes.shape[0]
    # Hover/colorbar trace at the centroids
    assert np.allclose(hover.marker.color, values)
    assert np.allclose(hover.x, coords[nodes, 0].mean(axis=1))
    assert hover.marker.showscale
    assert hover.marker.colorbar.title.text == "v"


def test_draw_element_values_displaced():
    coords, edof, nodes = quad_mesh(dofs_per_node=2)
    a = np.tile([0.5, 0.0], coords.shape[0])
    fig = cfv.draw_element_values(np.ones(nodes.shape[0]), coords, edof, 2, 3,
                                  displacements=a, magnfac=2.0,
                                  draw_elements=False)
    hover = fig.data[-1]
    assert np.allclose(hover.x, coords[nodes, 0].mean(axis=1) + 1.0)


def test_draw_element_values_3d():
    coords, edof, _ = hex_mesh()
    fig = cfv.draw_element_values([1.0, 2.0], coords, edof, 3, 5, clim=(0, 4))
    mesh3d = fig.data[0]
    assert mesh3d.intensitymode == "cell"
    # Each triangle gets the value of its element, 12 triangles per element
    intensity = np.asarray(mesh3d.intensity)
    tri_x = np.asarray(mesh3d.x)[np.column_stack((mesh3d.i, mesh3d.j,
                                                  mesh3d.k))].mean(axis=1)
    assert np.allclose(intensity[tri_x < 1.0], 1.0)
    assert np.allclose(intensity[tri_x > 1.0], 2.0)
    assert np.sum(intensity == 1.0) == 12 and np.sum(intensity == 2.0) == 12
    assert (mesh3d.cmin, mesh3d.cmax) == (0, 4)


def test_draw_displacements_auto_scale():
    coords, edof, nodes = quad_mesh(dofs_per_node=2)
    a = np.zeros(2*coords.shape[0])
    a[0] = 1e-3                                   # node 0 moves in x
    fig = cfv.draw_displacements(a, coords, edof, 2, 3,
                                 draw_undisplaced_mesh=True)
    undisplaced, displaced = fig.data
    # Largest displacement = 0.1 * model size (3.0)
    assert np.isclose(displaced.x[0] - undisplaced.x[0], 0.3)


def test_draw_displacements_given_magnfac():
    coords, edof, nodes = quad_mesh(dofs_per_node=2)
    a = np.tile([0.0, 0.1], coords.shape[0])
    fig = cfv.draw_displacements(a, coords, edof, 2, 3, magnfac=10.0)
    assert np.isclose(fig.data[0].y[0], coords[nodes[0][0], 1] + 1.0)


def test_colorbar_title():
    coords, edof, nodes = quad_mesh()
    cfv.draw_element_values(np.arange(nodes.shape[0]), coords, edof, 1, 3)
    cfv.colorbar("MPa")
    assert cfv.gcf().data[-1].marker.colorbar.title.text == "MPa"


# ------------------------------------------------------- flux

def test_draw_element_flux():
    coords, edof, nodes = quad_mesh()
    flux = np.tile([2.0, -1.0], (nodes.shape[0], 1))
    scale = cfv.draw_element_flux(flux, coords, edof, 1, 3, scale=0.25)
    arrows, hover = cfv.gcf().data
    assert scale == 0.25
    # First arrow: shaft centred at the centroid, length scale*|q|
    c = coords[nodes[0]].mean(axis=0)
    assert np.allclose([arrows.x[0], arrows.y[0]], c - [0.25, -0.125])
    assert np.allclose([arrows.x[1], arrows.y[1]], c + [0.25, -0.125])
    assert np.allclose(hover.customdata[:, 3], np.sqrt(5))


def test_draw_element_flux_matches_vis_mpl_scale():
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import calfem.vis_mpl as cfv_mpl
    coords, edof, nodes = quad_mesh(nx=5, ny=3)
    flux = np.random.default_rng(0).standard_normal((nodes.shape[0], 2))
    s_plotly = cfv.draw_element_flux(flux, coords, edof, 1, 3)
    s_mpl = cfv_mpl.draw_element_flux(flux, coords, edof, 1, 3)
    matplotlib.pyplot.close("all")
    assert np.isclose(s_plotly, s_mpl)


def test_draw_element_flux_color_by_magnitude():
    coords, edof, nodes = quad_mesh()
    flux = np.column_stack((np.arange(nodes.shape[0]) + 1.0,
                            np.zeros(nodes.shape[0])))
    cfv.draw_element_flux(flux, coords, edof, 1, 3, color_by_magnitude=True,
                          draw_elements=True)
    hover = cfv.gcf().data[-1]
    assert hover.marker.showscale
    assert np.allclose(hover.marker.color, flux[:, 0])


def test_draw_element_flux_invalid_input():
    coords, edof, nodes = quad_mesh()
    with pytest.raises(ValueError):
        cfv.draw_element_flux(np.zeros((2, 2)), coords, edof, 1, 3)
    coords3d = np.column_stack((coords, np.zeros(coords.shape[0])))
    with pytest.raises(ValueError):
        cfv.draw_element_flux(np.zeros((nodes.shape[0], 2)), coords3d,
                              edof, 1, 3)


# ------------------------------------------------------- nodal values

@pytest.mark.parametrize("func, trace_type, coloring", [
    (cfv.draw_nodal_values_contourf, go.Contour, "fill"),
    (cfv.draw_nodal_values_contour, go.Contour, "lines"),
    (cfv.draw_nodal_values_shaded, go.Heatmap, None),
])
def test_draw_nodal_values(func, trace_type, coloring):
    coords, edof, nodes = quad_mesh(nx=4, ny=3)
    values = coords[:, 0] + 2*coords[:, 1]     # linear, interpolated exactly
    fig = func(values, coords, edof, resolution=50)
    trace = fig.data[0]
    assert isinstance(trace, trace_type)
    if coloring:
        assert trace.contours.coloring == coloring
    X, Y = np.meshgrid(trace.x, trace.y)
    Z = np.asarray(trace.z, float)
    assert not np.isnan(Z).any()               # whole grid inside the mesh
    assert np.allclose(Z, X + 2*Y)


def test_draw_nodal_values_levels_and_mesh():
    coords, edof, nodes = quad_mesh()
    values = coords[:, 0]
    fig = cfv.draw_nodal_values_contourf(values, coords, edof, levels=6,
                                         dofs_per_node=1, el_type=3,
                                         draw_elements=True)
    contour, mesh = fig.data
    assert contour.contours.start == 0.0
    assert contour.contours.end == 3.0
    assert np.isclose(contour.contours.size, 0.5)
    assert n_polylines(mesh) == nodes.shape[0]


def test_draw_nodal_values_outside_mesh_is_nan():
    """Grid points in a hole (here a missing element) are NaN."""
    coords, edof, nodes = quad_mesh(nx=3, ny=3)
    edof = np.delete(edof, 4, axis=0)           # remove the center element
    fig = cfv.draw_nodal_values_shaded(coords[:, 0], coords, edof,
                                       resolution=31)
    trace = fig.data[0]
    X, Y = np.meshgrid(trace.x, trace.y)
    Z = np.asarray(trace.z, float)
    center = (X > 1.1) & (X < 1.9) & (Y > 0.8) & (Y < 1.2)
    assert np.isnan(Z[center]).all()
    assert not np.isnan(Z[~center & ((X < 0.9) | (X > 2.1))]).any()


def test_draw_nodal_values_missing_mesh_arguments():
    coords, edof, _ = quad_mesh()
    with pytest.raises(ValueError):
        cfv.draw_nodal_values_contourf(coords[:, 0], coords, edof,
                                       draw_elements=True)


# ------------------------------------------------------- figure numbers

def test_numbered_figures(monkeypatch):
    shown = []
    monkeypatch.setattr(cfv, "_show_in_browser",
                        lambda figs, *a, **k: shown.extend(figs))
    f1 = cfv.figure(1)
    f2 = cfv.figure(2)
    assert cfv.figure(1) is f1
    assert cfv.gcf() is f1
    cfv.show()
    assert len(shown) == 2 and shown[0] is f1 and shown[1] is f2
    # After show() the numbers start over with new figures
    assert cfv.figure(1) is not f1


def test_empty_figures_are_all_shown(monkeypatch):
    """Figures are tracked by identity, plotly figures compare by content."""
    shown = []
    monkeypatch.setattr(cfv, "_show_in_browser",
                        lambda figs, *a, **k: shown.extend(figs))
    cfv.figure()
    cfv.figure()
    cfv.show()
    assert len(shown) == 2


def test_close_numbered_figure(monkeypatch):
    shown = []
    monkeypatch.setattr(cfv, "_show_in_browser",
                        lambda figs, *a, **k: shown.extend(figs))
    f1 = cfv.figure(1)
    cfv.figure(2)
    cfv.close(2)
    cfv.show()
    assert len(shown) == 1 and shown[0] is f1


def test_axis_limits():
    cfv.axis([-1, 5, -2, 3])
    fig = cfv.gcf()
    assert tuple(fig.layout.xaxis.range) == (-1, 5)
    assert tuple(fig.layout.yaxis.range) == (-2, 3)


# ------------------------------------------------------- classic CALFEM

def element_arrays(nx=3, ny=2):
    coords, edof, nodes = quad_mesh(nx, ny)
    return coords[nodes, 0], coords[nodes, 1], nodes, coords


def test_pltstyle2():
    assert cfv.pltstyle2([1, 2, 1]) == ("blue", "solid", "black",
                                        "circle-open")
    assert cfv.pltstyle2([2, 4, 0]) == ("red", "dash", None, None)
    with pytest.raises(ValueError):
        cfv.pltstyle2([4, 1, 1])
    with pytest.raises(ValueError):
        cfv.pltstyle2([1, 5, 1])


def test_eldraw2_with_element_numbers():
    ex, ey, nodes, _ = element_arrays()
    nel = ex.shape[0]
    fig = cfv.eldraw2(ex, ey, [2, 4, 1], elnum=np.arange(1, nel + 1))
    lines, marks, numbers = fig.data
    assert n_polylines(lines) == nel
    assert lines.line.color == "red" and lines.line.dash == "dash"
    assert len(marks.x) == ex.size
    assert list(numbers.text) == [str(i) for i in range(1, nel + 1)]
    assert np.allclose(numbers.x, ex.mean(axis=1))


def test_eldraw2_single_element_and_no_marks():
    fig = cfv.eldraw2(np.array([0.0, 2.0]), np.array([0.0, 1.0]), [1, 1, 0])
    assert len(fig.data) == 1
    assert list(fig.data[0].x[:2]) == [0.0, 2.0]


def test_eldraw2_8_node_outline_order():
    ex = np.array([[0, 2, 2, 0, 1, 2, 1, 0]], float)
    ey = np.array([[0, 0, 2, 2, 0, 1, 2, 1]], float)
    fig = cfv.eldraw2(ex, ey, [1, 1, 0])
    x = np.asarray(fig.data[0].x[:9], float)
    assert np.allclose(x, [0, 1, 2, 2, 2, 1, 0, 0, 0])


def test_scalfact2():
    ex, ey, _, _ = element_arrays()
    ed = np.array([[0.0, -0.5], [0.25, 0.1]])
    assert np.isclose(cfv.scalfact2(ex, ey, ed, 0.2), 0.2*3.0/0.5)


def test_eldisp2_quads():
    ex, ey, nodes, coords = element_arrays()
    a = np.column_stack((0.01*coords[:, 0], np.zeros(coords.shape[0])))
    ed = a[nodes].reshape(nodes.shape[0], -1)
    sfac = cfv.eldisp2(ex, ey, ed, [1, 2, 1], sfac=10.0)
    lines, marks = cfv.gcf().data
    assert sfac == 10.0
    assert np.allclose(marks.x, (ex*1.1).ravel())
    assert n_polylines(lines) == ex.shape[0]


def test_eldisp2_auto_scale():
    ex, ey, nodes, coords = element_arrays()
    ed = np.zeros((nodes.shape[0], 8))
    ed[0, 1] = -2e-3
    sfac = cfv.eldisp2(ex, ey, ed)
    assert np.isclose(sfac, 0.1*3.0/2e-3)


def test_eldisp2_bars_and_beams():
    ex = np.array([[0.0, 1.0], [1.0, 1.0]])
    ey = np.array([[0.0, 0.0], [0.0, 1.0]])
    ed_bar = np.array([[0, 0, 0.1, 0.0], [0.1, 0.0, 0.1, 0.05]])
    cfv.eldisp2(ex, ey, ed_bar, [1, 1, 0], sfac=1.0)
    x = np.asarray(cfv.gcf().data[0].x, float)
    assert np.allclose(x[:2], [0.0, 1.1])

    cfv.close_all()
    ed_beam = np.array([[0, 0, 0, 0, 0.1, 0], [0, 0.1, 0, 0, 0.1, 0]])
    cfv.eldisp2(ex, ey, ed_beam, [1, 1, 1], sfac=2.0)
    curve, marks = cfv.gcf().data
    y = np.asarray(curve.y, float)
    # Deflected beam shape ends in the displaced node
    first = y[:np.argmax(np.isnan(y))]
    assert np.isclose(first[0], 0.0) and np.isclose(first[-1], 0.2)
    assert len(first) > 2
    assert np.allclose(marks.y, [0.0, 0.2, 0.2, 1.2])


def test_eldisp2_8_node_edges_pass_through_nodes():
    ex = np.array([[0, 2, 2, 0, 1, 2, 1, 0]], float)
    ey = np.array([[0, 0, 2, 2, 0, 1, 2, 1]], float)
    ed = np.zeros((1, 16))
    ed[0, 9] = 0.5                         # midside node 5 moves in y
    cfv.eldisp2(ex, ey, ed, [1, 1, 0], sfac=1.0)
    trace = cfv.gcf().data[0]
    pts = np.column_stack((trace.x, trace.y)).astype(float)
    pts = pts[np.isfinite(pts[:, 0])]
    for xn, yn in [(1, 0.5), (0, 0), (2, 2)]:
        assert np.min(np.hypot(pts[:, 0] - xn, pts[:, 1] - yn)) < 1e-12


def test_eldisp2_unsupported():
    with pytest.raises(ValueError):
        cfv.eldisp2(np.zeros((1, 5)), np.zeros((1, 5)), np.ones((1, 10)))


def test_dispbeam2():
    edi = np.column_stack((np.zeros(5), [0, -1, -2, -1, 0]))
    sfac = cfv.dispbeam2([0, 4], [1, 1], edi, [1, 1, 1], sfac=0.1)
    curve, marks = cfv.gcf().data
    assert sfac == 0.1
    assert np.allclose(curve.x, [0, 1, 2, 3, 4])
    assert np.allclose(curve.y, 1 + 0.1*edi[:, 1])
    assert np.isclose(cfv.dispbeam2([0, 4], [0, 0], edi), 0.1*4/2)


def test_secforce2():
    es = np.array([2.0, 1.0, -1.0])
    sfac = cfv.secforce2([0, 2], [0, 0], es, [2, 1])
    diagram, hover, element = cfv.gcf().data
    assert np.isclose(sfac, 0.2*2/2.0)
    # Positive values drawn below a beam along the x-axis
    assert np.allclose(hover.y, -sfac*es)
    assert np.allclose(hover.customdata[:, 1], es)
    assert list(element.x) == [0, 2]
    assert diagram.line.color == "blue"
    with pytest.raises(ValueError):
        cfv.secforce2([0, 2], [0, 0], es, [2, 1], eci=[0, 1])


def test_scalgraph2():
    cfv.scalgraph2(2.0, [5, 1, -1], plotpar=4)
    line, label = cfv.gcf().data
    assert line.x[:2] == (1, 11.0)
    assert line.line.color == "red"
    assert label.text == ("5",)
    with pytest.raises(ValueError):
        cfv.scalgraph2(1.0, [1, 2])


def test_elflux2_matches_draw_element_flux():
    ex, ey, nodes, coords = element_arrays()
    edof = nodes + 1
    flux = np.random.default_rng(3).standard_normal((nodes.shape[0], 2))
    s1 = cfv.elflux2(ex, ey, flux, plotcolor=[2])
    arrows1, _ = cfv.gcf().data
    cfv.close_all()
    s2 = cfv.draw_element_flux(flux, coords, edof, 1, 3)
    arrows2, _ = cfv.gcf().data
    assert np.isclose(s1, s2)
    assert np.allclose(np.asarray(arrows1.x, float),
                       np.asarray(arrows2.x, float), equal_nan=True)
    assert arrows1.line.color == "blue"


def test_eliso2():
    ex, ey, nodes, coords = element_arrays()
    ed = 2*ex + ey                         # linear field, exact
    fig = cfv.eliso2(ex, ey, ed, isov=5)
    trace = fig.data[0]
    assert trace.contours.coloring == "lines"
    X, Y = np.meshgrid(trace.x, trace.y)
    assert np.allclose(np.asarray(trace.z, float), 2*X + Y)

    cfv.close_all()
    fig = cfv.eliso2(ex, ey, ed, isov=[1, 2, 3], plotpar=[2, 4])
    trace = fig.data[0]
    assert trace.showscale is False
    assert trace.line.dash == "dash"
    assert cfv.eliso2_mpl is cfv.eliso2


# ------------------------------------------------------- geometry

def test_draw_geometry():
    g = cfg.Geometry()
    for p in [[0, 0], [2, 0], [2, 1], [0, 1], [1, 0.5], [1.2, 0.5],
              [1, 0.7], [0.8, 0.5], [1, 0.3]]:
        g.point(p)
    g.spline([0, 1])
    g.spline([1, 2], marker=7)
    g.spline([2, 3])
    g.spline([3, 0])
    for arc in ([5, 4, 6], [6, 4, 7], [7, 4, 8], [8, 4, 5]):
        g.circle(arc)
    fig = cfv.draw_geometry(g, title="Geometry")
    curves, points, labels = fig.data
    assert n_polylines(curves) == 8
    assert len(points.x) == 9
    # Arc points lie on the circle
    xy = np.array([curves.x, curves.y], dtype=float).T
    arcs = xy[np.isfinite(xy[:, 0]) & (np.abs(xy[:, 1] - 0.5) < 0.21)
              & (np.abs(xy[:, 0] - 1) < 0.21)]
    assert np.allclose(np.hypot(arcs[:, 0] - 1, arcs[:, 1] - 0.5), 0.2)
    assert "1[7]" in labels.text
    assert fig.layout.title.text == "Geometry"
