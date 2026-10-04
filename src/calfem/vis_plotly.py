# -*- coding: utf-8 -*-
"""
CALFEM Visualisation module (plotly)

Interactive visualisation in Jupyter notebooks, Google Colab and the web
browser. The functions have the same names and arguments as the ones in
calfem.vis_mpl, so switching backend only requires changing the import::

    import calfem.vis_plotly as cfv

    cfv.figure()
    cfv.draw_element_values(von_mises, coords, edof, dofs_per_node, el_type)
    cfv.show()

The classic CALFEM functions working on element arrays, eldraw2, eldisp2,
dispbeam2, secforce2, scalgraph2, scalfact2, elflux2 and eliso2, are also
available with the same arguments as in MATLAB CALFEM and vis_mpl.

Like matplotlib there is a current figure that the draw functions add to.
Figures can also be numbered, figure(1), as in matplotlib. figure() starts
a new figure and show() displays all figures created since the last call
to show(): inline in notebooks, and in scripts as HTML files opened in the
web browser (the location of the files is printed). The current plotly
figure is available with gcf() for further customisation, and
save_figure() writes it to HTML or an image.

Requires plotly: pip install calfem-python[plotly]
"""

import numpy as np

try:
    import plotly.graph_objects as go
    from plotly.colors import sample_colorscale
except ImportError as e:
    raise ImportError(
        "calfem.vis_plotly requires plotly. "
        "Install it with: pip install calfem-python[plotly]"
    ) from e


_current_figure = None
_pending_figures = []
_numbered_figures = {}

DEFAULT_COLORSCALE = "Viridis"

_COLOR_CHARS = {
    "r": "red", "g": "green", "b": "blue", "y": "yellow",
    "c": "cyan", "m": "magenta", "k": "black", "w": "white",
}

_LINE_STYLES = {
    "solid": "solid", "-": "solid",
    "dashed": "dash", "--": "dash",
    "dotted": "dot", ":": "dot",
    "dashdot": "dashdot", "-.": "dashdot",
}


# ------------------------------------------------------- figure handling

def figure(figure=None, show=True, fig_size=None, title=None):
    """
    Create a new figure and make it the current figure.

    Parameters
    ----------
    figure : int or plotly.graph_objects.Figure, optional
        Figure number, as in matplotlib: an existing figure with this
        number is made current, otherwise a new one is created. A plotly
        figure is made current.
    show : bool, optional
        Include the figure when show() is called. Default True.
    fig_size : tuple, optional
        Figure size (width, height) in inches, as in vis_mpl (100 pixels
        per inch). Default is plotly's automatic size.
    title : str, optional
        Figure title.

    Returns
    -------
    fig : plotly.graph_objects.Figure
    """
    global _current_figure

    number = None
    if isinstance(figure, (int, np.integer)):
        number = int(figure)
        figure = _numbered_figures.get(number)

    if figure is None:
        figure = go.Figure()
        figure.update_layout(
            template="plotly_white",
            showlegend=False,
            margin=dict(l=40, r=40, t=60, b=40),
        )
        if number is not None:
            _numbered_figures[number] = figure
    if fig_size is not None:
        figure.update_layout(width=int(fig_size[0]*100),
                             height=int(fig_size[1]*100))
    if title is not None:
        figure.update_layout(title=title)

    _current_figure = figure
    if show and not _contains(_pending_figures, figure):
        _pending_figures.append(figure)
    return figure


def gcf():
    """Return the current figure, creating one if needed."""
    if _current_figure is None:
        figure()
    return _current_figure


def clf():
    """Remove all traces and annotations from the current figure."""
    fig = gcf()
    fig.data = []
    fig.layout.annotations = []


def _contains(figures, fig):
    """Membership by identity, plotly figures compare by content."""
    return any(f is fig for f in figures)


def close(fig=None):
    """Close a figure (default the current one) so show() skips it. fig
    can also be a figure number."""
    global _current_figure
    if isinstance(fig, (int, np.integer)):
        fig = _numbered_figures.get(int(fig))
    fig = _current_figure if fig is None else fig
    _pending_figures[:] = [f for f in _pending_figures if f is not fig]
    for number in [n for n, f in _numbered_figures.items() if f is fig]:
        del _numbered_figures[number]
    if fig is _current_figure:
        _current_figure = None


def close_all():
    """Close all figures."""
    global _current_figure
    _pending_figures.clear()
    _numbered_figures.clear()
    _current_figure = None


def _in_notebook():
    """True when running in a Jupyter kernel (Jupyter, Colab, VS Code
    notebooks, ...)."""
    try:
        from IPython import get_ipython
    except ImportError:
        return False
    shell = get_ipython()
    return shell is not None and "IPKernelApp" in shell.config


def _open_file(filename):
    """
    Open a local file with the default application (the web browser for
    HTML files).

    On Windows and macOS the operating system is asked directly. Python's
    webbrowser module is avoided there, since it follows the BROWSER
    environment variable, which e.g. the VS Code terminal sets to a helper
    script that may not work, and then fails silently.
    """
    import subprocess
    import sys
    import webbrowser

    if sys.platform == "win32":
        import os
        os.startfile(str(filename))
    elif sys.platform == "darwin":
        subprocess.run(["open", str(filename)], check=True)
    else:
        if not webbrowser.open(filename.resolve().as_uri()):
            raise RuntimeError("no web browser found")


def _show_in_browser(figures, directory=None, open_browser=True):
    """Write the figures to HTML files and open them in the web browser.
    Returns the file names."""
    import tempfile
    from pathlib import Path

    if directory is None:
        directory = tempfile.mkdtemp(prefix="calfem_plotly_")
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)

    filenames = []
    failed = None
    for i, fig in enumerate(figures, start=1):
        filename = directory / f"figure_{i}.html"
        # plotly.min.js is written once to the directory and shared by the
        # figures, which keeps the files small and works offline.
        fig.write_html(filename, include_plotlyjs="directory")
        filenames.append(filename)
        if open_browser and failed is None:
            try:
                _open_file(filename)
            except Exception as e:
                failed = e

    message = (f"calfem.vis_plotly: {len(filenames)} figure(s) written to "
               f"{directory}")
    if not open_browser:
        print(message + ".")
    elif failed is None:
        print(message + " and opened in the web browser.")
    else:
        print(message + f". They could not be opened automatically "
              f"({failed}), open them in a web browser.")
    return filenames


def show(renderer=None, directory=None):
    """
    Show all figures created since the last call to show().

    In Jupyter, Colab and VS Code notebooks the figures are shown inline.
    In scripts they are written to HTML files, which are opened in the web
    browser. The location of the files is printed, so that they can also be
    opened manually.

    Parameters
    ----------
    renderer : str, optional
        Use this plotly renderer instead, e.g. "browser", "notebook" or
        "png". The PLOTLY_RENDERER environment variable has the same effect.
    directory : str or Path, optional
        Directory for the HTML files when running a script. Default a new
        temporary directory.

    Returns
    -------
    filenames : list of Path
        The HTML files written, empty if the figures were shown by plotly.
    """
    import os
    global _current_figure

    figures = list(_pending_figures)
    _pending_figures.clear()
    _numbered_figures.clear()
    _current_figure = None

    if not figures:
        return []

    if renderer is None and "PLOTLY_RENDERER" not in os.environ \
            and not _in_notebook():
        return _show_in_browser(figures, directory)

    for fig in figures:
        fig.show(renderer=renderer)
    return []


def show_and_wait(renderer=None, directory=None):
    """Same as show(), for compatibility with vis_mpl."""
    return show(renderer, directory)


showAndWait = show_and_wait


def save_figure(filename, fig=None, **kwargs):
    """
    Save a figure (default the current one). A .html file gives an
    interactive page, other extensions (.png, .svg, .pdf, ...) a static
    image, which requires the kaleido package.
    """
    fig = gcf() if fig is None else fig
    filename = str(filename)
    if filename.lower().endswith((".html", ".htm")):
        fig.write_html(filename, **kwargs)
    else:
        fig.write_image(filename, **kwargs)


def title(text):
    """Set the title of the current figure."""
    gcf().update_layout(title=text)


def axis(option):
    """
    Change the axes of the current figure.

    Parameters
    ----------
    option : str or list
        "equal" for equal scaling, "off" to hide the axes, "on" to show
        them, or the axis limits [xmin, xmax, ymin, ymax].
    """
    fig = gcf()
    if not isinstance(option, str):
        xmin, xmax, ymin, ymax = option
        fig.update_xaxes(range=[xmin, xmax])
        fig.update_yaxes(range=[ymin, ymax])
    elif option == "equal":
        fig.update_yaxes(scaleanchor="x", scaleratio=1)
    elif option == "off":
        fig.update_xaxes(visible=False)
        fig.update_yaxes(visible=False)
        fig.update_scenes(xaxis_visible=False, yaxis_visible=False,
                          zaxis_visible=False)
    elif option == "on":
        fig.update_xaxes(visible=True)
        fig.update_yaxes(visible=True)
        fig.update_scenes(xaxis_visible=True, yaxis_visible=True,
                          zaxis_visible=True)


def colorbar(title=None, **kwargs):
    """
    Colorbars are added automatically by the draw functions. This sets the
    title (and other colorbar properties) of the most recent colorbar.
    """
    fig = gcf()
    for trace in reversed(fig.data):
        if _set_colorbar(trace, title, kwargs):
            return


def _set_colorbar(trace, title, kwargs):
    props = dict(kwargs)
    if title is not None:
        props["title"] = title
    if hasattr(trace, "colorbar") and trace.showscale is not False:
        trace.colorbar.update(props)
        return True
    marker = getattr(trace, "marker", None)
    if marker is not None and getattr(marker, "showscale", None) is True:
        marker.colorbar.update(props)
        return True
    return False


def text(text, pos, angle=0, color="black", font_size=12, **kwargs):
    """Add a text label at pos = [x, y] in the current figure."""
    gcf().add_annotation(x=pos[0], y=pos[1], text=text, textangle=-angle,
                         showarrow=False, font=dict(color=_color(color),
                                                    size=font_size),
                         **kwargs)


add_text = text


# ------------------------------------------------------- helpers

def _color(c):
    """Convert a vis_mpl color (RGB tuple 0-1, 'rgbycmkw' char or any
    plotly color) to a plotly color."""
    if c is None:
        return None
    if isinstance(c, str):
        return _COLOR_CHARS.get(c, c)
    c = tuple(c)
    if all(0 <= v <= 1 for v in c[:3]):
        rgb = [int(round(255*v)) for v in c[:3]]
    else:
        rgb = [int(v) for v in c[:3]]
    if len(c) == 4:
        return f"rgba({rgb[0]},{rgb[1]},{rgb[2]},{c[3]})"
    return f"rgb({rgb[0]},{rgb[1]},{rgb[2]})"


def _setup_2d(fig, title=None):
    fig.update_yaxes(scaleanchor="x", scaleratio=1)
    if title:
        fig.update_layout(title=title)


def _setup_3d(fig, title=None):
    fig.update_scenes(aspectmode="data")
    if fig.layout.scene.camera.eye.x is None:
        fig.update_scenes(camera_eye=dict(x=1.6, y=1.6, z=1.2))
    if title:
        fig.update_layout(title=title)


def _element_nodes(edof, dofs_per_node, el_type):
    """Zero-based node indices per element, corner nodes only."""
    edof = np.asarray(edof, dtype=int)
    nodes = (edof[:, 0::dofs_per_node] - 1) // dofs_per_node
    if el_type == 9:            # 6-node triangle
        nodes = nodes[:, :3]
    elif el_type == 16:         # 8-node quadrilateral
        nodes = nodes[:, :4]
    elif el_type not in (1, 2, 3, 4, 5):
        raise ValueError(f"Element type {el_type} not supported.")
    return nodes


# Faces of 3D elements (Gmsh node order)
_HEX_FACES = np.array([[0, 3, 2, 1], [0, 1, 5, 4], [4, 5, 6, 7],
                       [2, 6, 5, 1], [2, 3, 7, 6], [0, 4, 7, 3]])
_TET_FACES = np.array([[0, 1, 2], [0, 3, 2], [1, 3, 2], [0, 3, 1]])


def _element_faces(nodes, el_type):
    """Faces (rows of node indices) and owning element of each face."""
    if el_type == 5:
        faces = nodes[:, _HEX_FACES].reshape(-1, 4)
        owner = np.repeat(np.arange(nodes.shape[0]), 6)
    elif el_type == 4:
        faces = nodes[:, _TET_FACES].reshape(-1, 3)
        owner = np.repeat(np.arange(nodes.shape[0]), 4)
    else:
        faces = nodes
        owner = np.arange(nodes.shape[0])
    return faces, owner


def _triangulate(faces):
    """Fan triangulation of polygonal faces. Returns triangles and the
    face index of each triangle."""
    n = faces.shape[1]
    tris = np.vstack([faces[:, [0, k, k + 1]] for k in range(1, n - 1)])
    owner = np.tile(np.arange(faces.shape[0]), n - 2)
    return tris, owner


def _polylines(polys, closed=True):
    """x, y(, z) arrays for an (m, n, dim) array of polylines, separated
    by NaN so that they can be drawn with a single trace."""
    polys = np.asarray(polys, dtype=float)
    if closed:
        polys = np.concatenate((polys, polys[:, :1, :]), axis=1)
    gap = np.full((polys.shape[0], 1, polys.shape[2]), np.nan)
    pts = np.concatenate((polys, gap), axis=1).reshape(-1, polys.shape[2])
    return [pts[:, k] for k in range(pts.shape[1])]


def _unique_edges(faces):
    """Unique edges of polygonal faces as an (m, 2) array."""
    n = faces.shape[1]
    edges = np.vstack([faces[:, [k, (k + 1) % n]] for k in range(n)])
    edges = np.unique(np.sort(edges, axis=1), axis=0)
    return edges


def _value_range(values, clim):
    if clim is not None:
        return float(clim[0]), float(clim[1])
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return 0.0, 1.0
    return float(finite.min()), float(finite.max())


def _color_bins(values, cmin, cmax, n_colors):
    if cmax > cmin:
        idx = np.floor((values - cmin)/(cmax - cmin)*n_colors).astype(int)
    else:
        idx = np.zeros(values.shape, dtype=int)
    return np.clip(idx, 0, n_colors - 1)


def _add_colored_polygons(fig, polys, values, colorscale, clim, n_colors,
                          edge_color=None, edge_width=1.0,
                          colorbar_title=None, hover_name="value"):
    """Filled polygons colored by values. The polygons are grouped into
    n_colors color bands, giving one trace per band. An invisible marker
    trace at the centroids carries the colorbar and hover information."""
    values = np.asarray(values, dtype=float).ravel()
    cmin, cmax = _value_range(values, clim)
    bins = _color_bins(values, cmin, cmax, n_colors)
    colors = sample_colorscale(colorscale,
                               list((np.arange(n_colors) + 0.5)/n_colors))

    for k in np.unique(bins):
        x, y = _polylines(polys[bins == k])
        fig.add_trace(go.Scatter(
            x=x, y=y, mode="lines", fill="toself", fillcolor=colors[k],
            line=dict(color=colors[k], width=1.0), hoverinfo="skip",
        ))

    if edge_color is not None:
        x, y = _polylines(polys)
        fig.add_trace(go.Scatter(
            x=x, y=y, mode="lines", hoverinfo="skip",
            line=dict(color=_color(edge_color), width=edge_width),
        ))

    centroids = polys.mean(axis=1)
    fig.add_trace(go.Scatter(
        x=centroids[:, 0], y=centroids[:, 1], mode="markers",
        marker=dict(size=8, opacity=0, color=values, colorscale=colorscale,
                    cmin=cmin, cmax=cmax, showscale=True,
                    colorbar=dict(title=colorbar_title)),
        customdata=np.arange(1, values.size + 1),
        hovertemplate=f"element %{{customdata}}<br>{hover_name} = "
                      "%{marker.color:.4g}<extra></extra>",
    ))


def _nodal_displacements(a, n_nodes, dofs_per_node, ndim):
    """Nodal displacements as an N-by-ndim array, from either the global
    displacement vector (dofs_per_node values per node) or an N-by-ndim
    array."""
    a = np.asarray(a, dtype=float)
    if dofs_per_node is not None and a.size == n_nodes*dofs_per_node:
        return a.reshape(n_nodes, dofs_per_node)[:, :ndim]
    return a.reshape(n_nodes, ndim)


def _auto_magnfac(coords, u, rel=0.1):
    size = np.max(coords.max(axis=0) - coords.min(axis=0))
    umax = np.max(np.linalg.norm(u, axis=1))
    return rel*size/umax if umax > 0 else 1.0


# ------------------------------------------------------- mesh drawing

def draw_mesh(
    coords,
    edof,
    dofs_per_node,
    el_type,
    title=None,
    color=(0, 0, 0),
    face_color=(0.8, 0.8, 0.8),
    node_color=(0, 0, 0),
    filled=False,
    show_nodes=False,
    line_width=1.0,
):
    """
    Draws the mesh in 2D or 3D.

    Parameters
    ----------
    coords : array_like
        An N-by-2 or N-by-3 array. Row i contains the x,y,z coordinates of
        node i.
    edof : array_like
        An E-by-L array. Element topology. (E is the number of elements and
        L is the number of dofs per element)
    dofs_per_node : int
        Dofs per node.
    el_type : int
        Element type (Gmsh numbering): 2 triangles, 3 quadrangles,
        4 tetrahedra, 5 hexahedra, 9 6-node triangles, 16 8-node quadrangles.
    title : str, optional
        Title of the figure.
    color : tuple or str, optional
        Color of the element edges. Default black.
    face_color : tuple or str, optional
        Color of the faces if filled is True (always used in 3D).
    node_color : tuple or str, optional
        Color of the nodes if show_nodes is True.
    filled : bool, optional
        Fill the elements with face_color. Default False.
    show_nodes : bool, optional
        Draw the nodes. Default False.
    line_width : float, optional
        Width of the element edges. Default 1.

    Returns
    -------
    fig : plotly.graph_objects.Figure
    """
    coords = np.asarray(coords, dtype=float)
    nodes = _element_nodes(edof, dofs_per_node, el_type)
    fig = gcf()

    if coords.shape[1] == 3 and el_type in (4, 5):
        _draw_mesh_3d(fig, coords, nodes, el_type, color, face_color,
                      line_width)
        if show_nodes:
            fig.add_trace(go.Scatter3d(
                x=coords[:, 0], y=coords[:, 1], z=coords[:, 2],
                mode="markers", marker=dict(size=2, color=_color(node_color)),
            ))
        _setup_3d(fig, title)
        return fig

    polys = coords[nodes][:, :, :2]
    x, y = _polylines(polys, closed=el_type != 1)
    fig.add_trace(go.Scatter(
        x=x, y=y, mode="lines", hoverinfo="skip",
        line=dict(color=_color(color), width=line_width),
        fill="toself" if filled else None,
        fillcolor=_color(face_color) if filled else None,
    ))
    if show_nodes:
        fig.add_trace(go.Scatter(
            x=coords[:, 0], y=coords[:, 1], mode="markers",
            marker=dict(size=4, color=_color(node_color)),
            customdata=np.arange(1, coords.shape[0] + 1),
            hovertemplate="node %{customdata}<br>(%{x:.4g}, %{y:.4g})"
                          "<extra></extra>",
        ))
    _setup_2d(fig, title)
    return fig


drawMesh = draw_mesh


def _draw_mesh_3d(fig, coords, nodes, el_type, color, face_color,
                  line_width, values=None, colorscale=None, clim=None,
                  colorbar_title=None):
    faces, face_owner = _element_faces(nodes, el_type)
    tris, tri_face = _triangulate(faces)

    mesh_args = dict(
        x=coords[:, 0], y=coords[:, 1], z=coords[:, 2],
        i=tris[:, 0], j=tris[:, 1], k=tris[:, 2],
        flatshading=True, hoverinfo="skip",
    )
    if values is None:
        mesh_args["color"] = _color(face_color)
    else:
        values = np.asarray(values, dtype=float).ravel()
        cmin, cmax = _value_range(values, clim)
        mesh_args.update(
            intensity=values[face_owner[tri_face]], intensitymode="cell",
            colorscale=colorscale, cmin=cmin, cmax=cmax, showscale=True,
            colorbar=dict(title=colorbar_title),
        )
    fig.add_trace(go.Mesh3d(**mesh_args))

    if color is not None:
        edges = _unique_edges(faces)
        x, y, z = _polylines(coords[edges], closed=False)
        fig.add_trace(go.Scatter3d(
            x=x, y=y, z=z, mode="lines", hoverinfo="skip",
            line=dict(color=_color(color), width=line_width*2),
        ))


def draw_elements(
    ex,
    ey,
    title="",
    color=(0, 0, 0),
    face_color=(0.8, 0.8, 0.8),
    node_color=(0, 0, 0),
    line_style="solid",
    filled=False,
    closed=True,
    show_nodes=False,
    line_width=1.0,
):
    """
    Draws elements given by element coordinate arrays ex, ey (one row per
    element), e.g. from coordxtr.

    Parameters
    ----------
    ex, ey : array_like
        Element x- and y-coordinates.
    title : str, optional
        Title of the figure.
    color : tuple or str, optional
        Color of the element edges. Default black.
    face_color : tuple or str, optional
        Color of the faces if filled is True.
    node_color : tuple or str, optional
        Color of the nodes if show_nodes is True.
    line_style : str, optional
        "solid", "dashed", "dotted" or "dashdot". Default "solid".
    filled : bool, optional
        Fill the elements. Default False.
    closed : bool, optional
        Draw the elements as closed polygons. Default True.
    show_nodes : bool, optional
        Draw the nodes. Default False.
    line_width : float, optional
        Width of the element edges. Default 1.

    Returns
    -------
    fig : plotly.graph_objects.Figure
    """
    ex = np.atleast_2d(np.asarray(ex, dtype=float))
    ey = np.atleast_2d(np.asarray(ey, dtype=float))
    polys = np.stack((ex, ey), axis=2)
    fig = gcf()

    x, y = _polylines(polys, closed=closed)
    fig.add_trace(go.Scatter(
        x=x, y=y, mode="lines", hoverinfo="skip",
        line=dict(color=_color(color), width=line_width,
                  dash=_LINE_STYLES.get(line_style, line_style)),
        fill="toself" if filled else None,
        fillcolor=_color(face_color) if filled else None,
    ))
    if show_nodes:
        fig.add_trace(go.Scatter(
            x=ex.ravel(), y=ey.ravel(), mode="markers", hoverinfo="skip",
            marker=dict(size=4, color=_color(node_color)),
        ))
    _setup_2d(fig, title)
    return fig


def draw_node_circles(
    ex,
    ey,
    title="",
    color=(0, 0, 0),
    face_color=(0.8, 0.8, 0.8),
    filled=False,
    marker_type="o",
):
    """Draws markers at the element nodes given by ex, ey."""
    symbols = {"o": "circle", "s": "square", "^": "triangle-up",
               "x": "x", "+": "cross", "d": "diamond"}
    symbol = symbols.get(marker_type, marker_type)
    if not filled:
        symbol += "-open"
    fig = gcf()
    fig.add_trace(go.Scatter(
        x=np.ravel(ex), y=np.ravel(ey), mode="markers", hoverinfo="skip",
        marker=dict(size=7, symbol=symbol, color=_color(color)),
    ))
    _setup_2d(fig, title)
    return fig


# ------------------------------------------------------- results

def draw_element_values(
    values,
    coords,
    edof,
    dofs_per_node,
    el_type,
    displacements=None,
    draw_elements=True,
    draw_undisplaced_mesh=False,
    magnfac=1.0,
    title=None,
    color=(0, 0, 0),
    node_color=(0, 0, 0),
    colorscale=DEFAULT_COLORSCALE,
    clim=None,
    colorbar_title=None,
    n_colors=64,
):
    """
    Draws scalar element values in 2D or 3D. Hovering an element shows its
    number and value.

    Parameters
    ----------
    values : array_like
        One scalar value per element.
    coords : array_like
        An N-by-2 or N-by-3 array. Row i contains the x,y,z coordinates of
        node i.
    edof : array_like
        An E-by-L array. Element topology.
    dofs_per_node : int
        Dofs per node.
    el_type : int
        Element type (Gmsh numbering), see draw_mesh.
    displacements : array_like, optional
        Nodal displacements, the global displacement vector or an N-by-2/3
        array. The elements are drawn in the displaced configuration.
    draw_elements : bool, optional
        Draw the element edges. Default True.
    draw_undisplaced_mesh : bool, optional
        Also draw the undisplaced mesh. Default False.
    magnfac : float, optional
        Magnification factor for the displacements. Default 1.
    title : str, optional
        Title of the figure.
    color : tuple or str, optional
        Color of the element edges.
    node_color : tuple or str, optional
        Not used, kept for compatibility with vis_mpl.
    colorscale : str or list, optional
        plotly colorscale. Default "Viridis".
    clim : tuple, optional
        Value range (min, max) of the colorscale. Default the data range.
    colorbar_title : str, optional
        Title of the colorbar.
    n_colors : int, optional
        Number of color bands used for 2D meshes. Default 64.

    Returns
    -------
    fig : plotly.graph_objects.Figure
    """
    coords = np.asarray(coords, dtype=float)
    nodes = _element_nodes(edof, dofs_per_node, el_type)
    fig = gcf()

    if draw_undisplaced_mesh:
        draw_mesh(coords, edof, dofs_per_node, el_type,
                  color=(0.6, 0.6, 0.6))

    ndim = coords.shape[1]
    if displacements is not None:
        u = _nodal_displacements(displacements, coords.shape[0],
                                 dofs_per_node, ndim)
        coords = coords + magnfac*u

    if ndim == 3 and el_type in (4, 5):
        _draw_mesh_3d(fig, coords, nodes, el_type,
                      color if draw_elements else None, None, 1.0,
                      values=values, colorscale=colorscale, clim=clim,
                      colorbar_title=colorbar_title)
        _setup_3d(fig, title)
        return fig

    polys = coords[nodes][:, :, :2]
    _add_colored_polygons(
        fig, polys, values, colorscale, clim, n_colors,
        edge_color=color if draw_elements else None,
        colorbar_title=colorbar_title,
    )
    _setup_2d(fig, title)
    return fig


drawElementValues = draw_element_values


def draw_displacements(
    a,
    coords,
    edof,
    dofs_per_node,
    el_type,
    draw_undisplaced_mesh=False,
    magnfac=-1.0,
    magscale=0.1,
    title=None,
    color=(0.3, 0.3, 0.3),
    node_color=(0, 0, 0),
):
    """
    Draws the displaced mesh in 2D or 3D.

    Parameters
    ----------
    a : array_like
        Global displacement vector or an N-by-2/3 array of nodal
        displacements.
    coords : array_like
        An N-by-2 or N-by-3 array of node coordinates.
    edof : array_like
        An E-by-L array. Element topology.
    dofs_per_node : int
        Dofs per node.
    el_type : int
        Element type (Gmsh numbering), see draw_mesh.
    draw_undisplaced_mesh : bool, optional
        Also draw the undisplaced mesh in light gray. Default False.
    magnfac : float, optional
        Magnification factor for the displacements. If negative (default),
        it is chosen so that the largest displacement is magscale times
        the model size.
    magscale : float, optional
        Largest displacement relative to the model size when magnfac is
        automatic. Default 0.1.
    title : str, optional
        Title of the figure.
    color : tuple or str, optional
        Color of the displaced mesh.
    node_color : tuple or str, optional
        Not used, kept for compatibility with vis_mpl.

    Returns
    -------
    fig : plotly.graph_objects.Figure
    """
    coords = np.asarray(coords, dtype=float)
    fig = gcf()

    if draw_undisplaced_mesh:
        draw_mesh(coords, edof, dofs_per_node, el_type,
                  color=(0.75, 0.75, 0.75))

    u = _nodal_displacements(a, coords.shape[0], dofs_per_node,
                             coords.shape[1])
    if magnfac is None or magnfac < 0:
        magnfac = _auto_magnfac(coords, u, magscale)

    draw_mesh(coords + magnfac*u, edof, dofs_per_node, el_type,
              color=color, title=title)
    return fig


drawDisplacements = draw_displacements


def _flux_scale_factor(ex, ey, flux, krel=0.8):
    nel = ex.shape[0]
    dx = np.max(ex, axis=1) - np.min(ex, axis=1)
    dy = np.max(ey, axis=1) - np.min(ey, axis=1)
    lm = np.sum(np.sqrt(dx**2 + dy**2))/nel
    qm = np.sum(np.sqrt(flux[:, 0]**2 + flux[:, 1]**2))/nel
    return 0.0 if qm == 0 else lm*krel/qm


def _arrow_lines(x0, y0, u, v, head_size=0.25, head_angle=0.4):
    """Arrows centered at (x0, y0) as NaN separated polylines."""
    xs, ys = x0 - u/2, y0 - v/2
    xe, ye = x0 + u/2, y0 + v/2
    ang = np.arctan2(v, u)
    length = np.hypot(u, v)*head_size
    xl = xe - length*np.cos(ang - head_angle)
    yl = ye - length*np.sin(ang - head_angle)
    xr = xe - length*np.cos(ang + head_angle)
    yr = ye - length*np.sin(ang + head_angle)
    nan = np.full_like(x0, np.nan)
    x = np.column_stack((xs, xe, xl, nan, xe, xr, nan)).ravel()
    y = np.column_stack((ys, ye, yl, nan, ye, yr, nan)).ravel()
    return x, y


def draw_element_flux(
    flux,
    coords,
    edof,
    dofs_per_node,
    el_type,
    scale=None,
    color=(0, 0, 0),
    color_by_magnitude=False,
    cmap=DEFAULT_COLORSCALE,
    draw_elements=False,
    element_color=(0.6, 0.6, 0.6),
    arrow_width=1.5,
    title=None,
    n_colors=32,
):
    """
    Draws element flux (flow) vectors as arrows at the element centroids
    of a 2D mesh. Hovering an arrow shows the flux components.

    Parameters
    ----------
    flux : array_like
        An E-by-2 array. Row i contains the flux vector [qx, qy] of
        element i.
    coords : array_like
        An N-by-2 array of node coordinates.
    edof : array_like
        An E-by-L array. Element topology.
    dofs_per_node : int
        Dofs per node.
    el_type : int
        Element type (Gmsh numbering), see draw_mesh.
    scale : float, optional
        Scale factor = arrow length / flux magnitude. Default automatic,
        mean arrow length 0.8 times the mean element size.
    color : tuple or str, optional
        Arrow color, if color_by_magnitude is False.
    color_by_magnitude : bool, optional
        Color the arrows by flux magnitude and add a colorbar.
    cmap : str or list, optional
        plotly colorscale used when color_by_magnitude is True.
    draw_elements : bool, optional
        Draw the mesh under the arrows. Default False.
    element_color : tuple or str, optional
        Color of the mesh if draw_elements is True.
    arrow_width : float, optional
        Line width of the arrows in pixels. Default 1.5.
    title : str, optional
        Title of the figure.
    n_colors : int, optional
        Number of color bands when color_by_magnitude is True.

    Returns
    -------
    scale : float
        Scale factor used for the arrows.
    """
    flux = np.asarray(flux, dtype=float)
    coords = np.asarray(coords, dtype=float)
    if coords.shape[1] != 2:
        raise ValueError("draw_element_flux only supports 2D meshes.")
    if flux.ndim != 2 or flux.shape[1] < 2:
        raise ValueError("flux must be an E-by-2 array of [qx, qy].")

    nodes = _element_nodes(edof, dofs_per_node, el_type)
    if flux.shape[0] != nodes.shape[0]:
        raise ValueError(
            "Check size of flux! There must be one row for each element."
        )

    ex = coords[nodes, 0]
    ey = coords[nodes, 1]
    if scale is None:
        scale = _flux_scale_factor(ex, ey, flux)

    fig = gcf()
    if draw_elements:
        draw_mesh(coords, edof, dofs_per_node, el_type, color=element_color)

    x0 = ex.mean(axis=1)
    y0 = ey.mean(axis=1)
    u = scale*flux[:, 0]
    v = scale*flux[:, 1]
    magnitude = np.hypot(flux[:, 0], flux[:, 1])

    if color_by_magnitude:
        cmin, cmax = _value_range(magnitude, None)
        bins = _color_bins(magnitude, cmin, cmax, n_colors)
        colors = sample_colorscale(cmap,
                                   list((np.arange(n_colors) + 0.5)/n_colors))
        for k in np.unique(bins):
            sel = bins == k
            x, y = _arrow_lines(x0[sel], y0[sel], u[sel], v[sel])
            fig.add_trace(go.Scatter(
                x=x, y=y, mode="lines", hoverinfo="skip",
                line=dict(color=colors[k], width=arrow_width),
            ))
        marker = dict(size=10, opacity=0, color=magnitude, colorscale=cmap,
                      cmin=cmin, cmax=cmax, showscale=True)
    else:
        x, y = _arrow_lines(x0, y0, u, v)
        fig.add_trace(go.Scatter(
            x=x, y=y, mode="lines", hoverinfo="skip",
            line=dict(color=_color(color), width=arrow_width),
        ))
        marker = dict(size=10, opacity=0)

    fig.add_trace(go.Scatter(
        x=x0, y=y0, mode="markers", marker=marker,
        customdata=np.column_stack((np.arange(1, x0.size + 1),
                                    flux[:, 0], flux[:, 1], magnitude)),
        hovertemplate="element %{customdata[0]}<br>"
                      "qx = %{customdata[1]:.4g}<br>"
                      "qy = %{customdata[2]:.4g}<br>"
                      "|q| = %{customdata[3]:.4g}<extra></extra>",
    ))
    _setup_2d(fig, title)
    return scale


# ------------------------------------------------------- nodal values

def _topo_to_tri(edof):
    """Triangles (zero-based) from a node based element topology."""
    nodes = np.asarray(edof, dtype=int) - 1
    n = nodes.shape[1]
    if n == 3:
        return nodes
    if n == 4:
        return np.vstack((nodes[:, [0, 1, 2]], nodes[:, [0, 2, 3]]))
    if n == 8:
        pattern = [[0, 4, 7], [4, 1, 5], [5, 2, 6], [6, 3, 7],
                   [4, 6, 7], [4, 5, 6]]
        return np.vstack([nodes[:, p] for p in pattern])
    raise ValueError("Element topology not supported.")


def _grid_values(values, coords, edof, resolution):
    """Interpolate nodal values on a regular grid, NaN outside the mesh."""
    from matplotlib.tri import Triangulation, LinearTriInterpolator

    coords = np.asarray(coords, dtype=float)
    x, y = coords[:, 0], coords[:, 1]
    tri = Triangulation(x, y, _topo_to_tri(edof))
    interp = LinearTriInterpolator(tri, np.asarray(values, float).ravel())

    width, height = np.ptp(x), np.ptp(y)
    if width >= height:
        nx = resolution
        ny = max(2, int(round(resolution*height/width)))
    else:
        ny = resolution
        nx = max(2, int(round(resolution*width/height)))
    # Grid points exactly on the boundary can be classified as outside the
    # mesh due to round-off, so the grid is moved slightly inwards.
    eps = 1e-9*max(width, height)
    xi = np.linspace(x.min() + eps, x.max() - eps, nx)
    yi = np.linspace(y.min() + eps, y.max() - eps, ny)
    X, Y = np.meshgrid(xi, yi)
    Z = interp(X, Y).filled(np.nan)
    return xi, yi, Z


def _contour_settings(values, levels):
    vmin, vmax = _value_range(np.asarray(values, float).ravel(), None)
    if np.ndim(levels) == 0:
        n = int(levels)
        size = (vmax - vmin)/n if vmax > vmin else 1.0
        return dict(start=vmin, end=vmax, size=size)
    levels = np.sort(np.asarray(levels, float))
    size = np.diff(levels).mean() if levels.size > 1 else 1.0
    return dict(start=levels[0], end=levels[-1], size=size)


def _draw_nodal_values(values, coords, edof, coloring, levels, title,
                       dofs_per_node, el_type, draw_elements, colorscale,
                       resolution, colorbar_title):
    fig = gcf()
    xi, yi, Z = _grid_values(values, coords, edof, resolution)

    if coloring == "heatmap":
        fig.add_trace(go.Heatmap(
            x=xi, y=yi, z=Z, colorscale=colorscale, zsmooth="best",
            colorbar=dict(title=colorbar_title),
            hovertemplate="(%{x:.4g}, %{y:.4g})<br>%{z:.4g}<extra></extra>",
        ))
    else:
        fig.add_trace(go.Contour(
            x=xi, y=yi, z=Z, colorscale=colorscale, connectgaps=False,
            autocontour=False,
            contours=dict(coloring=coloring,
                          **_contour_settings(values, levels)),
            line=dict(width=1.5 if coloring == "lines" else 0.5),
            colorbar=dict(title=colorbar_title),
            hovertemplate="(%{x:.4g}, %{y:.4g})<br>%{z:.4g}<extra></extra>",
        ))

    if draw_elements:
        if dofs_per_node is not None and el_type is not None:
            draw_mesh(coords, edof, dofs_per_node, el_type,
                      color=(0.2, 0.2, 0.2), line_width=0.5)
        else:
            raise ValueError("dofs_per_node and el_type must be specified "
                             "to draw the mesh.")
    _setup_2d(fig, title)
    return fig


def draw_nodal_values_contourf(
    values,
    coords,
    edof,
    levels=12,
    title=None,
    dofs_per_node=None,
    el_type=None,
    draw_elements=False,
    colorscale=DEFAULT_COLORSCALE,
    resolution=300,
    colorbar_title=None,
):
    """
    Draws a filled contour plot of nodal values on a 2D mesh.

    The values are interpolated linearly over the elements onto a regular
    grid with resolution points along the longest side, since plotly has no
    contour plot for unstructured meshes.

    Parameters
    ----------
    values : array_like
        One value per node.
    coords : array_like
        An N-by-2 array of node coordinates.
    edof : array_like
        Element topology with one dof per node, i.e. node numbers starting
        at 1 (3, 4 or 8 nodes per element).
    levels : int or array_like, optional
        Number of contour levels or the contour levels. Default 12.
    title : str, optional
        Title of the figure.
    dofs_per_node, el_type : int, optional
        Needed to draw the mesh if draw_elements is True.
    draw_elements : bool, optional
        Draw the mesh on top of the contours. Default False.
    colorscale : str or list, optional
        plotly colorscale. Default "Viridis".
    resolution : int, optional
        Number of grid points along the longest side. Default 300.
    colorbar_title : str, optional
        Title of the colorbar.

    Returns
    -------
    fig : plotly.graph_objects.Figure
    """
    return _draw_nodal_values(values, coords, edof, "fill", levels, title,
                              dofs_per_node, el_type, draw_elements,
                              colorscale, resolution, colorbar_title)


def draw_nodal_values_contour(
    values,
    coords,
    edof,
    levels=12,
    title=None,
    dofs_per_node=None,
    el_type=None,
    draw_elements=False,
    colorscale=DEFAULT_COLORSCALE,
    resolution=300,
    colorbar_title=None,
):
    """
    Draws contour lines of nodal values on a 2D mesh. See
    draw_nodal_values_contourf for the parameters.
    """
    return _draw_nodal_values(values, coords, edof, "lines", levels, title,
                              dofs_per_node, el_type, draw_elements,
                              colorscale, resolution, colorbar_title)


def draw_nodal_values_shaded(
    values,
    coords,
    edof,
    title=None,
    dofs_per_node=None,
    el_type=None,
    draw_elements=False,
    colorscale=DEFAULT_COLORSCALE,
    resolution=300,
    colorbar_title=None,
):
    """
    Draws smoothly shaded nodal values on a 2D mesh. See
    draw_nodal_values_contourf for the parameters.
    """
    return _draw_nodal_values(values, coords, edof, "heatmap", None, title,
                              dofs_per_node, el_type, draw_elements,
                              colorscale, resolution, colorbar_title)


draw_nodal_values = draw_nodal_values_contourf


# ------------------------------------------------------- geometry

def draw_geometry(
    geometry,
    draw_points=True,
    label_points=True,
    label_curves=True,
    title=None,
    font_size=11,
    N=20,
    rel_margin=0.05,
    draw_axis=False,
    axes=None,
):
    """
    Draws the geometry (points and curves) of a calfem.geometry.Geometry.

    Parameters
    ----------
    geometry : calfem.geometry.Geometry
        The geometry to draw.
    draw_points : bool, optional
        Draw the points. Default True.
    label_points : bool, optional
        Label the points as ID[marker]. Default True.
    label_curves : bool, optional
        Label the curves as ID(elements on curve)[marker]. Default True.
    title : str, optional
        Title of the figure.
    font_size : int, optional
        Size of the labels. Default 11.
    N : int, optional
        Number of points per curve segment. Default 20.
    rel_margin : float, optional
        Not used, plotly adds margins automatically.
    draw_axis : bool, optional
        Show the axes. Default False.
    axes : optional
        Not used, kept for compatibility with vis_mpl.

    Returns
    -------
    fig : plotly.graph_objects.Figure
    """
    from calfem.vis_mpl import (_catmullspline, _bspline, _circleArc,
                                _ellipseArc)

    fig = gcf()
    is_3d = geometry.is3D
    curve_x, curve_y, curve_z = [], [], []
    labels = []

    for ID, (curve_name, point_ids, marker, el_on_curve, _, _) in \
            geometry.curves.items():
        points = np.asarray(geometry.getPointCoords(point_ids), dtype=float)
        if curve_name == "Spline":
            P = _catmullspline(points, N)
        elif curve_name == "BSpline":
            P = _bspline(points, N)
        elif curve_name == "Circle":
            P = _circleArc(*points, pointsOnCurve=N)
        elif curve_name == "Ellipse":
            P = _ellipseArc(*points, pointsOnCurve=N)
        else:
            continue
        P = np.asarray(P)
        curve_x += list(P[:, 0]) + [None]
        curve_y += list(P[:, 1]) + [None]
        if P.shape[1] > 2:
            curve_z += list(P[:, 2]) + [None]

        label = str(ID)
        label += f"({el_on_curve})" if el_on_curve is not None else ""
        label += f"[{marker}]" if marker != 0 else ""
        labels.append((P[int(P.shape[0]*7.0/12), :], label))

    point_xyz = np.asarray(geometry.getPointCoords(), dtype=float)
    point_labels = [
        str(ID) + (f"[{marker}]" if marker != 0 else "")
        for ID, (xyz, el_size, marker) in geometry.points.items()
    ]

    if is_3d:
        fig.add_trace(go.Scatter3d(
            x=curve_x, y=curve_y, z=curve_z, mode="lines",
            line=dict(color="black", width=3), hoverinfo="skip",
        ))
        if draw_points:
            fig.add_trace(go.Scatter3d(
                x=point_xyz[:, 0], y=point_xyz[:, 1], z=point_xyz[:, 2],
                mode="markers+text" if label_points else "markers",
                text=point_labels, marker=dict(size=4, color="red"),
                textfont=dict(size=font_size, color="purple"),
            ))
        _setup_3d(fig, title)
        return fig

    fig.add_trace(go.Scatter(
        x=curve_x, y=curve_y, mode="lines", hoverinfo="skip",
        line=dict(color="black", width=1.5),
    ))
    if draw_points:
        fig.add_trace(go.Scatter(
            x=point_xyz[:, 0], y=point_xyz[:, 1],
            mode="markers+text" if label_points else "markers",
            text=point_labels, textposition="top right",
            textfont=dict(size=font_size, color="purple"),
            marker=dict(size=7, color="#1f77b4"),
            hovertemplate="point %{text}<br>(%{x:.4g}, %{y:.4g})"
                          "<extra></extra>",
        ))
    if label_curves and labels:
        fig.add_trace(go.Scatter(
            x=[p[0] for p, _ in labels], y=[p[1] for p, _ in labels],
            mode="text", text=[t for _, t in labels],
            textposition="middle right", textfont=dict(size=font_size),
            hoverinfo="skip",
        ))
    if not draw_axis:
        fig.update_xaxes(visible=False)
        fig.update_yaxes(visible=False)
    _setup_2d(fig, title)
    return fig


drawGeometry = draw_geometry


# ------------------------------------------------------- classic CALFEM
#
# Functions with the same names and arguments as in MATLAB CALFEM and
# calfem.vis_mpl, working on element coordinate arrays ex, ey (one row per
# element, e.g. from coordxtr) and plotpar codes:
#
#   linetype  1 solid, 2 dashed, 3 dotted
#   linecolor 1 black, 2 blue, 3 magenta, 4 red
#   nodemark  0 none, 1 circle, 2 star, 3 point

_PLOTPAR_COLORS = {1: "black", 2: "blue", 3: "magenta", 4: "red"}
_PLOTPAR_DASH = {1: "solid", 2: "dash", 3: "dot"}
_PLOTPAR_MARKERS = {0: None, 1: ("circle-open", 7), 2: ("asterisk-open", 8),
                    3: ("circle", 4)}


def _plotpar_color(code, name="plotpar"):
    if code not in _PLOTPAR_COLORS:
        raise ValueError(f"Invalid color code {code} in {name}, "
                         "use 1 (black), 2 (blue), 3 (magenta) or 4 (red).")
    return _PLOTPAR_COLORS[code]


def pltstyle2(plotpar):
    """
    Translate CALFEM plotpar = [linetype, linecolor, nodemark] to plotly
    styles.

    Returns
    -------
    line_color : str
    line_style : str
        plotly dash style.
    node_color : str or None
    node_symbol : str or None
        plotly marker symbol, None for no node marks.
    """
    if len(plotpar) != 3:
        raise ValueError("plotpar needs to contain 3 values.")
    p1, p2, p3 = plotpar
    if p1 not in _PLOTPAR_DASH:
        raise ValueError("Invalid value for plotpar[0].")
    if p3 not in _PLOTPAR_MARKERS:
        raise ValueError("Invalid value for plotpar[2].")
    line_color = _plotpar_color(p2, "plotpar[1]")
    marker = _PLOTPAR_MARKERS[p3]
    if marker is None:
        return line_color, _PLOTPAR_DASH[p1], None, None
    return line_color, _PLOTPAR_DASH[p1], "black", marker[0]


def _element_rows(ex, ey, *arrays):
    """ex, ey (and other element arrays) as 2D arrays, one row per
    element."""
    ex = np.asarray(ex, dtype=float)
    ey = np.asarray(ey, dtype=float)
    if ex.shape != ey.shape:
        raise ValueError("Check size of ex, ey dimensions.")
    single = ex.ndim == 1
    out = [np.atleast_2d(ex), np.atleast_2d(ey)]
    for a in arrays:
        a = np.asarray(a, dtype=float)
        out.append(a.reshape(1, -1) if single else np.atleast_2d(a))
    return out


def _add_node_marks(fig, x, y, plotpar_mark):
    marker = _PLOTPAR_MARKERS[plotpar_mark]
    if marker is None:
        return
    symbol, size = marker
    fig.add_trace(go.Scatter(
        x=np.ravel(x), y=np.ravel(y), mode="markers", hoverinfo="skip",
        marker=dict(symbol=symbol, size=size, color="black",
                    line=dict(width=1)),
    ))


def eldraw2(ex, ey, plotpar=[1, 2, 1], elnum=None):
    """
    Draw the undeformed 2D mesh for a number of elements of the same type.

    Supported elements are bars and beams (2 nodes), triangles (3 nodes),
    quadrilaterals (4 nodes) and 8-node isoparametric elements.

    Parameters
    ----------
    ex, ey : array_like
        Element node coordinates, one row per element.
    plotpar : list, optional
        [linetype, linecolor, nodemark]. Default [1, 2, 1], solid blue
        lines with circles at the nodes.
    elnum : array_like, optional
        Element numbers, drawn at the element centers.

    Returns
    -------
    fig : plotly.graph_objects.Figure
    """
    ex, ey = _element_rows(ex, ey)
    line_color, line_style, _, _ = pltstyle2(plotpar)
    nen = ex.shape[1]
    # 8-node elements: corner and midside nodes in order along the boundary
    order = [0, 4, 1, 5, 2, 6, 3, 7] if nen == 8 else list(range(nen))

    fig = gcf()
    x, y = _polylines(np.stack((ex[:, order], ey[:, order]), axis=2),
                      closed=nen > 2)
    fig.add_trace(go.Scatter(
        x=x, y=y, mode="lines", hoverinfo="skip",
        line=dict(color=line_color, dash=line_style, width=1),
    ))
    _add_node_marks(fig, ex, ey, plotpar[2])

    if elnum is not None and len(elnum) > 0:
        if len(elnum) != ex.shape[0]:
            raise ValueError("elnum must contain one number per element.")
        fig.add_trace(go.Scatter(
            x=ex.mean(axis=1), y=ey.mean(axis=1), mode="text",
            text=[str(int(n)) for n in np.ravel(elnum)], hoverinfo="skip",
            textfont=dict(color=line_color, size=11),
        ))
    _setup_2d(fig)
    return fig


def scalfact2(ex, ey, ed, rat=0.2):
    """
    Determine a scale factor for drawing computational results, such as
    displacements, section forces or flux: rat * largest model extent /
    largest absolute value in ed.
    """
    ex, ey = _element_rows(ex, ey)
    dl_max = max(np.ptp(ex), np.ptp(ey))
    ed_max = float(np.max(np.abs(ed)))
    return rat*dl_max/ed_max


def _quad8_edges(x, y, n_points=7):
    """Points along the (curved) edges of 8-node elements, one row per
    element, from the element shape functions. n_points per edge should be
    odd, so that the midside nodes are included."""
    def shape(t, s):
        return np.array([
            -0.25*(1 - t)*(1 - s)*(1 + t + s),
            -0.25*(1 + t)*(1 - s)*(1 - t + s),
            -0.25*(1 + t)*(1 + s)*(1 - t - s),
            -0.25*(1 - t)*(1 + s)*(1 + t - s),
            0.5*(1 - t*t)*(1 - s),
            0.5*(1 + t)*(1 - s*s),
            0.5*(1 - t*t)*(1 + s),
            0.5*(1 - t)*(1 - s*s),
        ])

    g = np.linspace(-1, 1, n_points)
    path = ([(t, -1) for t in g[:-1]] + [(1, s) for s in g[:-1]]
            + [(t, 1) for t in g[::-1][:-1]] + [(-1, s) for s in g[::-1]])
    N = np.array([shape(t, s) for t, s in path])      # (npts, 8)
    return x @ N.T, y @ N.T


def _beam2crd(ex, ey, ed, mag):
    from calfem.core import beam2crd
    excd, eycd = beam2crd(ex, ey, ed, mag)
    return np.atleast_2d(excd), np.atleast_2d(eycd)


def eldisp2(ex, ey, ed, plotpar=[2, 1, 1], sfac=None):
    """
    Draw the deformed 2D mesh for a number of elements of the same type.

    Supported elements are bars (2 nodes, 4 dofs), beams (2 nodes, 6 dofs,
    drawn with the deflected shape), triangles (3 nodes), quadrilaterals
    (4 nodes) and 8-node isoparametric elements (curved edges).

    Parameters
    ----------
    ex, ey : array_like
        Element node coordinates, one row per element.
    ed : array_like
        Element displacements, one row per element.
    plotpar : list, optional
        [linetype, linecolor, nodemark]. Default [2, 1, 1], dashed black
        lines with circles at the nodes.
    sfac : float, optional
        Scale factor for the displacements. Default automatic, the largest
        displacement is 0.1 times the largest model extent.

    Returns
    -------
    sfac : float
        Scale factor used.
    """
    ex, ey, ed = _element_rows(ex, ey, ed)
    if ed.shape[0] != ex.shape[0]:
        raise ValueError("Check size of ed/ex dimensions.")
    line_color, line_style, _, _ = pltstyle2(plotpar)
    nen = ex.shape[1]
    ned = ed.shape[1]

    if sfac is None:
        ed_max = float(np.max(np.abs(ed)))
        dl_max = max(np.ptp(ex), np.ptp(ey))
        sfac = 0.1*dl_max/ed_max if ed_max > 0 else 1.0
    k = sfac

    if nen == 2 and ned == 4:                           # bars
        x = ex + k*ed[:, [0, 2]]
        y = ey + k*ed[:, [1, 3]]
        xc, yc, closed = x, y, False
    elif nen == 2 and ned == 6:                         # beams
        x = ex + k*ed[:, [0, 3]]
        y = ey + k*ed[:, [1, 4]]
        xc, yc = _beam2crd(ex, ey, ed, k)
        closed = False
    elif nen in (3, 4) and ned == 2*nen:                # triangles, quads
        x = ex + k*ed[:, 0::2]
        y = ey + k*ed[:, 1::2]
        xc, yc, closed = x, y, True
    elif nen == 8 and ned == 16:                        # 8-node elements
        x = ex + k*ed[:, 0::2]
        y = ey + k*ed[:, 1::2]
        xc, yc = _quad8_edges(x, y)
        closed = False
    else:
        raise ValueError("Element type is not supported.")

    fig = gcf()
    px, py = _polylines(np.stack((xc, yc), axis=2), closed=closed)
    fig.add_trace(go.Scatter(
        x=px, y=py, mode="lines", hoverinfo="skip",
        line=dict(color=line_color, dash=line_style, width=1),
    ))
    _add_node_marks(fig, x, y, plotpar[2])
    _setup_2d(fig)
    return sfac


def dispbeam2(ex, ey, edi, plotpar=[2, 1, 1], sfac=None):
    """
    Draw the displacement diagram for a two dimensional beam element.

    Parameters
    ----------
    ex, ey : array_like
        Element node coordinates [x1, x2], [y1, y2].
    edi : array_like
        Displacements [[u1, v1], [u2, v2], ...] in local coordinates in
        evenly distributed points along the beam, e.g. from beam2s.
    plotpar : list, optional
        [linetype, linecolor, nodemark]. Default [2, 1, 1].
    sfac : float, optional
        Scale factor. Default automatic, the largest displacement is 0.1
        times the element length.

    Returns
    -------
    sfac : float
        Scale factor used.
    """
    ex = np.asarray(ex, dtype=float).ravel()
    ey = np.asarray(ey, dtype=float).ravel()
    edi = np.asarray(edi, dtype=float)
    if ex.shape != ey.shape:
        raise ValueError("Check size of ex, ey dimensions.")
    if edi.ndim != 2 or edi.shape[1] != 2:
        raise ValueError("Check size of edi dimension.")
    line_color, line_style, _, _ = pltstyle2(plotpar)

    d = np.array([ex[1] - ex[0], ey[1] - ey[0]])
    L = np.hypot(*d)
    n = d/L
    if sfac is None:
        sfac = 0.1*L/np.max(np.abs(edi))

    s = np.linspace(0.0, L, edi.shape[0])
    u = sfac*edi
    xc = ex[0] + s*n[0] + u[:, 0]*n[0] - u[:, 1]*n[1]
    yc = ey[0] + s*n[1] + u[:, 0]*n[1] + u[:, 1]*n[0]

    fig = gcf()
    fig.add_trace(go.Scatter(
        x=xc, y=yc, mode="lines", hoverinfo="skip",
        line=dict(color=line_color, dash=line_style, width=1),
    ))
    _add_node_marks(fig, xc[[0, -1]], yc[[0, -1]], plotpar[2])
    _setup_2d(fig)
    return sfac


def secforce2(ex, ey, es, plotpar=[2, 1], sfac=None, eci=None):
    """
    Draw the section force diagram for a two dimensional bar or beam
    element. Hovering the diagram shows the values.

    Parameters
    ----------
    ex, ey : array_like
        Element node coordinates [x1, x2], [y1, y2].
    es : array_like
        Section force in evaluation points along the element, e.g. one
        column of the section forces from beam2s.
    plotpar : list, optional
        [linecolor, elementcolor]. Default [2, 1], blue diagram on a black
        element.
    sfac : float, optional
        Scale factor. Default automatic, the largest section force is 0.2
        times the element length.
    eci : array_like, optional
        Local x-coordinates of the evaluation points. Default evenly
        distributed.

    Returns
    -------
    sfac : float
        Scale factor used.
    """
    ex = np.asarray(ex, dtype=float).ravel()
    ey = np.asarray(ey, dtype=float).ravel()
    es = np.asarray(es, dtype=float).ravel()
    if ex.shape != ey.shape:
        raise ValueError("Check size of ex, ey dimensions.")
    line_color = _plotpar_color(plotpar[0], "plotpar[0]")
    element_color = _plotpar_color(plotpar[1], "plotpar[1]")

    d = np.array([ex[1] - ex[0], ey[1] - ey[0]])
    L = np.hypot(*d)
    n = d/L
    if sfac is None:
        es_max = np.max(np.abs(es))
        sfac = 0.2*L/es_max if es_max > 0 else 1.0
    if eci is None:
        eci = np.linspace(0.0, L, es.size)
    eci = np.asarray(eci, dtype=float).ravel()
    if eci.size != es.size:
        raise ValueError("Check size of eci dimension.")

    # Points on the element and on the diagram
    xb = ex[0] + eci*n[0]
    yb = ey[0] + eci*n[1]
    xd = xb + sfac*es*n[1]
    yd = yb - sfac*es*n[0]

    nan = np.full(es.size, np.nan)
    fig = gcf()
    fig.add_trace(go.Scatter(                          # diagram and stripes
        x=np.concatenate((xd, [np.nan],
                          np.column_stack((xb, xd, nan)).ravel())),
        y=np.concatenate((yd, [np.nan],
                          np.column_stack((yb, yd, nan)).ravel())),
        mode="lines", hoverinfo="skip",
        line=dict(color=line_color, width=1),
    ))
    fig.add_trace(go.Scatter(                          # hover values
        x=xd, y=yd, mode="markers", marker=dict(size=6, opacity=0),
        customdata=np.column_stack((eci, es)),
        hovertemplate="x = %{customdata[0]:.4g}<br>"
                      "value = %{customdata[1]:.4g}<extra></extra>",
    ))
    fig.add_trace(go.Scatter(                          # element
        x=ex, y=ey, mode="lines", hoverinfo="skip",
        line=dict(color=element_color, width=2),
    ))
    _setup_2d(fig)
    return sfac


def scalgraph2(sfac, magnitude, plotpar=2):
    """
    Draw a graphic scale.

    Parameters
    ----------
    sfac : float
        Scale factor, e.g. from eldisp2, secforce2 or scalfact2.
    magnitude : array_like
        [Ref] or [Ref, x, y]. The scale has a length equivalent to Ref and
        starts at (x, y), default (0, -0.5).
    plotpar : int, optional
        Line color, 1 black, 2 blue (default), 3 magenta, 4 red.

    Returns
    -------
    fig : plotly.graph_objects.Figure
    """
    magnitude = np.ravel(magnitude)
    if magnitude.size == 1:
        N, x, y = magnitude[0], 0.0, -0.5
    elif magnitude.size == 3:
        N, x, y = magnitude
    else:
        raise ValueError("Check size of magnitude input argument.")
    color = _plotpar_color(plotpar)

    L = N*sfac
    h = L/20
    fig = gcf()
    fig.add_trace(go.Scatter(
        x=[x, x + L, None, x, x, None, x + L, x + L],
        y=[y, y, None, y - h, y + h, None, y - h, y + h],
        mode="lines", hoverinfo="skip", line=dict(color=color, width=1),
    ))
    fig.add_trace(go.Scatter(
        x=[x + 1.1*L], y=[y], mode="text", text=[f"{N:g}"],
        textposition="middle right", hoverinfo="skip",
    ))
    _setup_2d(fig)
    return fig


def elflux2(ex, ey, es, plotcolor=None, sfac=None, ax=None):
    """
    Display element flow arrows for 2D triangular or quadrilateral
    elements. Hovering an arrow shows the flow components.

    Parameters
    ----------
    ex, ey : array_like
        Element node coordinates, one row per element (3 or 4 nodes).
    es : array_like
        Element flow vectors [qx, qy], one row per element.
    plotcolor : list, optional
        [arrowcolor], 1 black (default), 2 blue, 3 magenta, 4 red.
    sfac : float, optional
        Scale factor = arrow length / flow magnitude. Default automatic.
    ax : optional
        Not used, kept for compatibility with vis_mpl.

    Returns
    -------
    sfac : float
        Scale factor used.

    See Also
    --------
    draw_element_flux : Equivalent function using coords/edof input.
    """
    ex, ey, es = _element_rows(ex, ey, es)
    if es.shape[0] != ex.shape[0]:
        raise ValueError("Check size of flow input argument! "
                         "There must be one row for each element.")
    if ex.shape[1] not in (3, 4):
        raise ValueError("Only 3- and 4-node elements are supported.")
    if es.shape[1] < 2:
        raise ValueError("es must contain the x- and y-components of the "
                         "flow.")

    color = _PLOTPAR_COLORS.get(int((plotcolor or [1])[0]), "black")
    if sfac is None:
        sfac = _flux_scale_factor(ex, ey, es)

    x0 = ex.mean(axis=1)
    y0 = ey.mean(axis=1)
    x, y = _arrow_lines(x0, y0, sfac*es[:, 0], sfac*es[:, 1])

    fig = gcf()
    fig.add_trace(go.Scatter(
        x=x, y=y, mode="lines", hoverinfo="skip",
        line=dict(color=color, width=1.5),
    ))
    fig.add_trace(go.Scatter(
        x=x0, y=y0, mode="markers", marker=dict(size=10, opacity=0),
        customdata=es[:, :2],
        hovertemplate="qx = %{customdata[0]:.4g}<br>"
                      "qy = %{customdata[1]:.4g}<extra></extra>",
    ))
    _setup_2d(fig)
    return sfac


def eliso2(ex, ey, ed, isov=10, plotpar=None, resolution=300):
    """
    Draw isolines from element nodal values for 2D triangular,
    quadrilateral or 8-node elements.

    Parameters
    ----------
    ex, ey : array_like
        Element node coordinates, one row per element.
    ed : array_like
        Element nodal values, one row per element, e.g. element
        temperatures from extract_eldisp.
    isov : int or array_like, optional
        Number of isolines or the isoline values. Default 10.
    plotpar : list, optional
        [linetype, linecolor] for single colored isolines. Default
        isolines colored by value with a colorbar.
    resolution : int, optional
        Number of grid points along the longest side, see
        draw_nodal_values_contour.

    Returns
    -------
    fig : plotly.graph_objects.Figure
    """
    ex, ey, ed = _element_rows(ex, ey, ed)
    if ed.shape != ex.shape:
        raise ValueError("ed must contain one value per element node.")

    # Merge the element nodes into a node based mesh
    points = np.column_stack((ex.ravel(), ey.ravel()))
    tol = 1e-9*max(np.ptp(points[:, 0]), np.ptp(points[:, 1]))
    keys = np.round(points/tol).astype(np.int64)
    _, first, inverse = np.unique(keys, axis=0, return_index=True,
                                  return_inverse=True)
    coords = points[first]
    values = ed.ravel()[first]
    topo = inverse.reshape(ex.shape) + 1

    fig = _draw_nodal_values(values, coords, topo, "lines", isov, None,
                             None, None, False, DEFAULT_COLORSCALE,
                             resolution, None)
    if plotpar is not None:
        color = _plotpar_color(plotpar[1], "plotpar[1]")
        fig.data[-1].update(colorscale=[[0, color], [1, color]],
                            showscale=False,
                            line=dict(dash=_PLOTPAR_DASH.get(plotpar[0],
                                                             "solid")))
    return fig


# Name used in vis_mpl
eliso2_mpl = eliso2
