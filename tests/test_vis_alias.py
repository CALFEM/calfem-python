# -*- coding: utf-8 -*-
"""
Tests for calfem.vis, an alias for calfem.vis_mpl, and the deprecated
visvis based module calfem.vis_visvis.
"""

import ast
import importlib
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

import calfem.vis as cfv
import calfem.vis_mpl as cfv_mpl


# Names of the old visvis module that only make sense with visvis
VISVIS_ONLY = {"visApp", "showAndWaitNative",
               "waitDisplayNative"}


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


def quad_mesh():
    xs, ys = np.meshgrid(np.linspace(0, 3, 4), np.linspace(0, 2, 3))
    coords = np.column_stack((xs.ravel(), ys.ravel()))
    nodes = np.array([[j*4 + i, j*4 + i + 1, j*4 + i + 5, j*4 + i + 4]
                      for j in range(2) for i in range(3)])
    return coords, nodes


def test_vis_is_vis_mpl():
    assert cfv is cfv_mpl
    import calfem
    assert calfem.vis is cfv_mpl
    from calfem import vis
    assert vis is cfv_mpl


def test_vis_shares_module_state():
    coords, nodes = quad_mesh()
    cfv.figure()
    cfv.draw_element_values(np.arange(6.0), coords, nodes + 1, 1, 3)
    # The mappable set through calfem.vis is used by calfem.vis_mpl
    assert cfv_mpl.colorbar() is not None


def test_old_vis_api_available():
    tree = ast.parse((src_dir / "calfem" / "vis_visvis.py")
                     .read_text(encoding="utf-8"))
    names = set()
    for node in tree.body:
        if isinstance(node, ast.FunctionDef):
            names.add(node.name)
        elif isinstance(node, ast.Assign):
            names.update(t.id for t in node.targets
                         if isinstance(t, ast.Name))
    public = {n for n in names if not n.startswith("_")} - VISVIS_ONLY
    missing = sorted(n for n in public if not hasattr(cfv, n))
    assert missing == []


def test_draw_nodal_values_compat():
    coords, nodes = quad_mesh()
    edof = np.column_stack([np.column_stack((2*nodes[:, k] + 1,
                                             2*nodes[:, k] + 2))
                            for k in range(4)])   # 2 dofs per node
    values = coords[:, 0]
    _, ax = plt.subplots()
    tpc = cfv.drawNodalValues(values, coords, edof, 2, 3, clim=(0, 6),
                              title="T")
    assert tpc.get_clim() == (0, 6)
    assert ax.get_title() == "T"
    assert tpc.get_array().size == coords.shape[0]
    cfv.colorbar()


def test_elval2():
    coords, nodes = quad_mesh()
    ex, ey = coords[nodes, 0], coords[nodes, 1]
    _, ax = plt.subplots()
    pc = cfv.elval2(ex, ey, np.arange(6.0), showMesh=True)
    assert len(pc.get_paths()) == 6
    assert np.allclose(pc.get_array(), np.arange(6.0))


def test_compat_aliases():
    assert cfv.drawElementValues is cfv.draw_element_values
    assert cfv.drawGeometry is cfv.draw_geometry
    assert cfv.eldraw2_mpl is cfv.eldraw2
    cfv.figure()
    cfv.addLabel("A", [0, 0])
    cfv.showGrid(True)


def test_vis_visvis_is_deprecated():
    sys.modules.pop("calfem.vis_visvis", None)
    try:
        import visvis  # noqa: F401
        have_visvis = True
    except ImportError:
        have_visvis = False

    if have_visvis:
        with pytest.warns(FutureWarning, match="deprecated"):
            importlib.import_module("calfem.vis_visvis")
    else:
        with pytest.warns(FutureWarning, match="deprecated"):
            with pytest.raises(ImportError, match=r"calfem-python\[visvis\]"):
                importlib.import_module("calfem.vis_visvis")


def test_shapes_does_not_need_visvis():
    sys.modules.pop("calfem.shapes", None)
    importlib.import_module("calfem.shapes")
    tree = ast.parse((src_dir / "calfem" / "shapes.py")
                     .read_text(encoding="utf-8"))
    imported = {a.name for n in ast.walk(tree)
                if isinstance(n, ast.Import) for a in n.names}
    assert "calfem.vis" not in imported
