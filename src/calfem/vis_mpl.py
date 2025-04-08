# -*- coding: utf-8 -*-
"""
CALFEM Visualisation module (matplotlib)

Contains all the functions implementing visualisation routines.
"""

import os

from matplotlib.transforms import Transform
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.collections
import matplotlib.path as mpp
import matplotlib.patches as patches
import matplotlib as mpl
import matplotlib.tri as tri

from calfem.core import beam2crd
import calfem.core as cfc

try:
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import (
        NavigationToolbar2QT as NavigationToolbar,
    )
except:
    print("Could not import Matplotlib backends. Probarbly due to missing Qt.")

from numpy import sin, cos, pi
from math import atan2

import logging as cflog

try:
    from numpy.lib.function_base import place
except:
    pass


g_figures = []


def error(msg):
    """Log error message"""
    cflog.error(msg)


def info(msg):
    """Log information message"""
    cflog.info(msg)


def figure_class():
    """Return visvis Figure class."""
    return None


figureClass = figure_class

cfv_def_mappable = None

def set_mappable(mappable):
    global cfv_def_mappable
    cfv_def_mappable = mappable

cfv_block_at_show = True

def set_block_at_show(show_block):
    cfv_block_at_show = show_block

def ioff():
    plt.ioff()

def ion():
    plt.ion()


def colorbar(**kwargs):
    """Add a colorbar to current figure"""
    global cfv_def_mappable
    if cfv_def_mappable != None:
        cbar = plt.colorbar(mappable=cfv_def_mappable, ax=plt.gca(), **kwargs)
        cfv_def_mappable = None
        return cbar
    else:
        return plt.colorbar(**kwargs)


def axis(*args, **kwargs):
    """Define axis of figure (Matplotlib passthrough)"""
    plt.axis(*args, **kwargs)


def title(*args, **kwargs):
    """Define title of figure (Matplotlib passthrough)"""
    plt.title(*args, **kwargs)


def figure(figure=None, show=True, fig_size=(6, 5.33)):
    """Create a visvis figure with extras."""
    f = None

    if figure == None:
        f = plt.figure(figsize=fig_size)
    else:
        try:
            f = plt.figure(figure.number)
        except:
            f = plt.figure(figsize=fig_size)

    if f is not None:
        g_figures.append(f)

    return f


def figure_widget(fig, parent=None):
    widget = FigureCanvas(fig)
    #widget.axes = fig.add_subplot(111)
    if parent != None:
        widget.setParent(parent)
    toolbar = NavigationToolbar(widget, widget)
    return widget


def close_all():
    """Close all visvis windows."""
    plt.close("all")


closeAll = close_all

def clf():
    """Clear visvis figure"""
    plt.clf()

def close(fig=None):
    """Close visvis figure"""
    if fig == None:
        plt.close()
    else:
        plt.close(fig)


def gca():
    """Get current axis of the current visvis figure."""
    return plt.gca()


def gcf():
    return plt.gcf()


def subplot(*args):
    """Create a visvis subplot."""
    return plt.subplot(*args)


def camera3d():
    """Get visvis 3D camera."""
    return None


def show_and_wait():
    """Wait for plot to show"""
    global cfv_block_at_show

    if "CFV_NO_BLOCK" in os.environ:
        plt.show(block=False)
    else:
        plt.show(block=cfv_block_at_show)


showAndWait = show_and_wait


def show_and_wait_mpl():
    """Wait for plot to show"""
    plt.show()


def show():
    """Use in Qt applications"""
    plt.show()


showAndWaitMpl = show_and_wait_mpl


def set_figure_dpi(dpi):
    mpl.rcParams["figure.dpi"] = dpi


def text(text, pos, angle=0, **kwargs):
    return plt.text(pos[0], pos[1], text, **kwargs)


add_text = text
addText = text
label = text


def ce2vf(coords, edof, dofs_per_node, el_type):
    """Duplicate code. Extracts verts, faces and verticesPerFace from input."""

    if np.shape(coords)[1] == 2:
        is_3d = False
        # pad with zeros to make 3D
        verts = np.hstack((coords, np.zeros([np.shape(coords)[0], 1])))
    elif np.shape(coords)[1] == 3:
        is_3d = True
        verts = coords
    else:
        raise ValueError("coords must be N-by-2 or N-by-3 array")

    if el_type in [2, 4]:  # elements with triangular faces
        vertices_per_face = 3
    elif el_type in [3, 5, 16]:  # elements with rectangular faces
        vertices_per_face = 4
    else:  # [NOTE] This covers all element types available in CALFEM plus tetrahedrons. If more element types are added it is necessary to include them here and below.
        raise ValueError("element type not implemented")

    faces = (edof[:, 0::dofs_per_node] - 1) / dofs_per_node
    # 'faces' here are actually lists of nodes in elements, not in faces necessarily if the elements are in 3D. This case is handled below.

    if el_type in [4, 5]:  # if hexahedrons or tetrahedrons:
        if el_type == 5:
            G = np.array(
                [
                    [0, 3, 2, 1],
                    [0, 1, 5, 4],
                    [4, 5, 6, 7],
                    [2, 6, 5, 1],
                    [2, 3, 7, 6],
                    [0, 4, 7, 3],
                ]
            )  # G is an array that is used to decomposes hexahedrons into its component faces.
        # The numbers are from the node orders (see p94 in the Gmsh manual) and each row makes one face.
        elif el_type == 4:
            G = np.array(
                [[0, 1, 2], [0, 3, 2], [1, 3, 2], [0, 3, 1]]
            )  # This G decomposes tetrahedrons into faces
        faces = np.vstack([faces[i, G] for i in range(faces.shape[0])])
    elif el_type == 16:  # if 8-node-quads:
        # The first 4 nodes are the corners of the high order quad.
        faces = faces[:, 0:4]

    return verts, np.asarray(faces, dtype=int), vertices_per_face, is_3d


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
):
    """
    Draws wire mesh of model in 2D or 3D. Returns the Mesh object that represents
    the mesh.
    Args:
        coords:
            An N-by-2 or N-by-3 array. Row i contains the x,y,z coordinates of node i.
        edof:
            An E-by-L array. Element topology. (E is the number of elements and L is the number of dofs per element)
        dofs_per_nodes:
            Integer. Dofs per node.
        el_type:
            Integer. Element Type. See Gmsh manual for details. Usually 2 for triangles or 3 for quadrangles.
        axes:
            Matplotlib Axes. The Axes where the model will be drawn. If unspecified the current Axes will be used, or a new Axes will be created if none exist.
        axes_adjust:
            Boolean. True if the view should be changed to show the whole model. Default True.
        title:
            String. Changes title of the figure. Default "Mesh".
        color:
            3-tuple or char. Color of the wire. Defaults to black (0,0,0). Can also be given as a character in 'rgbycmkw'.
        face_color:
            3-tuple or char. Color of the faces. Defaults to white (1,1,1). Parameter filled must be True or faces will not be drawn at all.
        filled:
            Boolean. Faces will be drawn if True. Otherwise only the wire is drawn. Default False.
    """

    verts, faces, vertices_per_face, is_3d = ce2vf(coords, edof, dofs_per_node, el_type)

    y = verts[:, 0]
    z = verts[:, 1]

    values = np.zeros(faces.shape[0], float)

    def quatplot(y, z, quatrangles, values=[], ax=None, **kwargs):
        if not ax:
            ax = plt.gca()
        yz = np.c_[y, z]
        v = yz[quatrangles]
        if filled:
            pc = matplotlib.collections.PolyCollection(
                v, facecolor=face_color, **kwargs
            )
        else:
            pc = matplotlib.collections.PolyCollection(v, facecolor="none", **kwargs)

        ax.add_collection(pc)
        ax.autoscale()
        return pc

    ax = plt.gca()
    ax.set_aspect("equal")

    pc = quatplot(y, z, faces, values, ax=ax, edgecolor=color)

    if show_nodes:
        ax.plot(y, z, marker="o", ls="", color=node_color)

    if title != None:
        ax.set(title=title)


drawMesh = draw_mesh


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
):
    """
    Draws wire mesh of model in 2D or 3D. Returns the Mesh object that represents
    the mesh.
    Args:
        coords:
            An N-by-2 or N-by-3 array. Row i contains the x,y,z coordinates of node i.
        edof:
            An E-by-L array. Element topology. (E is the number of elements and L is the number of dofs per element)
        dofs_per_nodes:
            Integer. Dofs per node.
        el_type:
            Integer. Element Type. See Gmsh manual for details. Usually 2 for triangles or 3 for quadrangles.
        axes:
            Matplotlib Axes. The Axes where the model will be drawn. If unspecified the current Axes will be used, or a new Axes will be created if none exist.
        axes_adjust:
            Boolean. True if the view should be changed to show the whole model. Default True.
        title:
            String. Changes title of the figure. Default "Mesh".
        color:
            3-tuple or char. Color of the wire. Defaults to black (0,0,0). Can also be given as a character in 'rgbycmkw'.
        face_color:
            3-tuple or char. Color of the faces. Defaults to white (1,1,1). Parameter filled must be True or faces will not be drawn at all.
        filled:
            Boolean. Faces will be drawn if True. Otherwise only the wire is drawn. Default False.
    """

    # ex = [
    #        [x1_1, x2_1, xn_1],
    #        ...
    #        [x1_m, x2_m, xn_m]
    #      ]

    # ex = [
    #        [y1_1, y2_1, yn_1],
    #        ...
    #        [y1_m, y2_m, yn_m]
    #      ]

    if ex.ndim != 1:
        nnodes = ex.shape[1]
        nel = ex.shape[0]
    else:
        nnodes = ex.shape[0]
        nel = 1

    polys = []

    for elx, ely in zip(ex, ey):
        pg = np.zeros((nnodes, 2))
        pg[:, 0] = elx
        pg[:, 1] = ely

        polys.append(pg)

    ax = plt.gca()

    if filled:
        pc = matplotlib.collections.PolyCollection(
            polys,
            facecolor=face_color,
            edgecolor=color,
            linestyle=line_style,
            closed=closed,
        )
    else:
        pc = matplotlib.collections.PolyCollection(
            polys,
            facecolor="none",
            edgecolor=color,
            linestyle=line_style,
            closed=closed,
        )

    ax.add_collection(pc)
    ax.autoscale()

    ax.set_aspect("equal")

    if title != None:
        ax.set(title=title)


def draw_node_circles(
    ex,
    ey,
    title="",
    color=(0, 0, 0),
    face_color=(0.8, 0.8, 0.8),
    filled=False,
    marker_type="o",
):
    """
    Draws wire mesh of model in 2D or 3D. Returns the Mesh object that represents
    the mesh.
    Args:
        coords:
            An N-by-2 or N-by-3 array. Row i contains the x,y,z coordinates of node i.
        edof:
            An E-by-L array. Element topology. (E is the number of elements and L is the number of dofs per element)
        dofs_per_nodes:
            Integer. Dofs per node.
        el_type:
            Integer. Element Type. See Gmsh manual for details. Usually 2 for triangles or 3 for quadrangles.
        axes:
            Matplotlib Axes. The Axes where the model will be drawn. If unspecified the current Axes will be used, or a new Axes will be created if none exist.
        axes_adjust:
            Boolean. True if the view should be changed to show the whole model. Default True.
        title:
            String. Changes title of the figure. Default "Mesh".
        color:
            3-tuple or char. Color of the wire. Defaults to black (0,0,0). Can also be given as a character in 'rgbycmkw'.
        face_color:
            3-tuple or char. Color of the faces. Defaults to white (1,1,1). Parameter filled must be True or faces will not be drawn at all.
        filled:
            Boolean. Faces will be drawn if True. Otherwise only the wire is drawn. Default False.
    """

    # ex = [
    #        [x1_1, x2_1, xn_1],
    #        ...
    #        [x1_m, x2_m, xn_m]
    #      ]

    # ex = [
    #        [y1_1, y2_1, yn_1],
    #        ...
    #        [y1_m, y2_m, yn_m]
    #      ]

    nel = ex.shape[0]
    nnodes = ex.shape[1]

    nodes = []

    x = []
    y = []

    for elx, ely in zip(ex, ey):
        for xx, yy in zip(elx, ely):
            x.append(xx)
            y.append(yy)

    ax = plt.gca()

    if filled:
        ax.scatter(x, y, color=color, marker=marker_type)
    else:
        ax.scatter(x, y, edgecolor=color, color="none", marker=marker_type)

    ax.autoscale()

    ax.set_aspect("equal")

    if title != None:
        ax.set(title=title)


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
):
    """
    Draws scalar element values in 2D or 3D.

    Args:
        values:
            An N-by-1 array or a list of scalars. The Scalar values of the elements. ev[i] should be the value of element i.

        coords:
            An N-by-2 or N-by-3 array. Row i contains the x,y,z coordinates of node i.

        edof:
            An E-by-L array. Element topology. (E is the number of elements and L is the number of dofs per element)

        dofs_per_node:
            Integer. Dofs per node.

        el_type:
            Integer. Element Type. See Gmsh manual for details. Usually 2 for triangles or 3 for quadrangles.

        displacements:
            An N-by-2 or N-by-3 array. Row i contains the x,y,z displacements of node i.

        draw_elements:
            Boolean. True if mesh wire should be drawn. Default True.

        draw_undisplaced_mesh:
            Boolean. True if the wire of the undisplaced mesh should be drawn on top of the displaced mesh. Default False. Use only if displacements != None.

        magnfac:
            Float. Magnification factor. Displacements are multiplied by this value. Use this to make small displacements more visible.

        title:
            String. Changes title of the figure. Default "Element Values".
    """

    if draw_undisplaced_mesh:
        draw_mesh(coords, edof, dofs_per_node, el_type, color=(0.5, 0.5, 0.5))

    if displacements is not None:
        if displacements.shape[1] != coords.shape[1]:
            displacements = np.reshape(displacements, (-1, coords.shape[1]))
            coords = np.asarray(coords + magnfac * displacements)

    verts, faces, vertices_per_face, is_3d = ce2vf(coords, edof, dofs_per_node, el_type)

    y = verts[:, 0]
    z = verts[:, 1]

    def quatplot(y, z, quatrangles, values=[], ax=None, **kwargs):
        if not ax:
            ax = plt.gca()
        yz = np.c_[y, z]
        v = yz[quatrangles]
        pc = matplotlib.collections.PolyCollection(v, **kwargs)

        pc.set_array(np.asarray(values))
        ax.add_collection(pc)
        ax.autoscale()
        return pc

    fig = plt.gcf()
    ax = plt.gca()
    ax.set_aspect("equal")

    if draw_elements:
        pc = quatplot(y, z, faces, values, ax=ax, edgecolor=color)
    else:
        pc = quatplot(y, z, faces, values, ax=ax, edgecolor=None)

    # pc = quatplot(y,z, np.asarray(edof-1), values, ax=ax,
    #         edgecolor="crimson", cmap="rainbow")

    set_mappable(pc)

    if title != None:
        ax.set(title=title)


def draw_displacements(
    a,
    coords,
    edof,
    dofs_per_node,
    el_type,
    draw_undisplaced_mesh=False,
    magnfac=-1.0,
    magscale=0.25,
    title=None,
    color=(0, 0, 0),
    node_color=(0, 0, 0),
):
    """
    Draws scalar element values in 2D or 3D. Returns the world object
    elementsWobject that represents the mesh.

    Args:
        ev:
            An N-by-1 array or a list of scalars. The Scalar values of the elements. ev[i] should be the value of element i.
        coords:
            An N-by-2 or N-by-3 array. Row i contains the x,y,z coordinates of node i.
        edof:
            An E-by-L array. Element topology. (E is the number of elements and L is the number of dofs per element)
        dofs_per_node:
            Integer. Dofs per node.
        el_type:
            Integer. Element Type. See Gmsh manual for details. Usually 2 for triangles or 3 for quadrangles.
        displacements:
            An N-by-2 or N-by-3 array. Row i contains the x,y,z  displacements of node i.
        axes:
            Matlotlib Axes. The Axes where the model will be drawn. If unspecified the current Axes will be used, or a new Axes will be created if none exist.
        draw_undisplaced_mesh:
            Boolean. True if the wire of the undisplaced mesh should be drawn on top of the displaced mesh. Default False. Use only if displacements != None.
        magnfac:
            Float. Magnification factor. Displacements are multiplied by this value. Use this to make small displacements more visible.
        title:
            String. Changes title of the figure. Default "Element Values".
    """

    if draw_undisplaced_mesh:
        draw_mesh(coords, edof, dofs_per_node, el_type, color=(0.8, 0.8, 0.8))

    if a is not None:
        if a.shape[1] != coords.shape[1]:
            a = np.reshape(a, (-1, coords.shape[1]))

            x_max = np.max(coords[:, 0])
            x_min = np.min(coords[:, 0])

            y_max = np.max(coords[:, 1])
            y_min = np.min(coords[:, 1])

            x_size = x_max - x_min
            y_size = y_max - y_min

            if x_size > y_size:
                max_size = x_size
            else:
                max_size = y_size

            if magnfac < 0:
                magnfac = 0.25 * max_size

            coords = np.asarray(coords + magnfac * a)

    verts, faces, vertices_per_face, is_3d = ce2vf(coords, edof, dofs_per_node, el_type)

    y = verts[:, 0]
    z = verts[:, 1]

    values = []

    def quatplot(y, z, quatrangles, values=[], ax=None, **kwargs):
        if not ax:
            ax = plt.gca()
        yz = np.c_[y, z]
        v = yz[quatrangles]
        pc = matplotlib.collections.PolyCollection(v, **kwargs)

        ax.add_collection(pc)
        ax.autoscale()
        return pc

    ax = plt.gca()
    ax.set_aspect("equal")

    pc = quatplot(
        y, z, faces, values, ax=ax, edgecolor=(0.3, 0.3, 0.3), facecolor="none"
    )

    if title != None:
        ax.set(title=title)


def create_ordered_polys(geom, N=10):
    """Creates ordered polygons from the geometry definition"""

    N = 10

    o_polys = []

    for id, (surf_name, curve_ids, holes, _, _, _) in geom.surfaces.items():
        polygon = np.empty((0, 3), float)

        polys = []

        for curve_id in curve_ids:
            curve_name, curve_points, _, _, _, _ = geom.curves[curve_id]
            points = geom.get_point_coords(curve_points)

            if curve_name == "Spline":
                P = _catmullspline(points, N)
            if curve_name == "BSpline":
                P = _bspline(points, N)
            if curve_name == "Circle":
                P = _circleArc(*points, pointsOnCurve=N)
            if curve_name == "Ellipse":
                P = _ellipseArc(*points, pointsOnCurve=N)

            polys.append(P)

        ordered_polys = []

        ordered_polys.append(polys.pop())

        while len(polys) != 0:
            p0 = ordered_polys[-1]
            for p in polys:
                if np.allclose(p0[-1], p[0]):
                    ordered_polys.append(polys.pop())
                    break
                elif np.allclose(p0[-1], p[-1]):
                    ordered_polys.append(np.flipud(polys.pop()))
                    break

        for p in ordered_polys:
            polygon = np.concatenate((polygon, p))

        o_polys.append(polygon)

    return o_polys


def draw_ordered_polys(o_polys):
    for poly in o_polys:
        ax = plt.gca()
        path = mpp.Path(poly[:, 0:2])
        patch = patches.PathPatch(path, facecolor="orange", lw=1)
        ax.add_patch(patch)


def point_in_geometry(o_polys, point):
    for poly in o_polys:
        path = mpp.Path(poly[:, 0:2])
        inside = path.contains_points([point])

        if inside:
            return True

    return False


def topo_to_tri(edof):
    """Converts 2d element topology to triangle topology to be used
    with the matplotlib functions tricontour and tripcolor."""

    if edof.shape[1] == 3:
        return edof
    elif edof.shape[1] == 4:
        new_edof = np.zeros((edof.shape[0] * 2, 3), int)
        new_edof[0::2, 0] = edof[:, 0]
        new_edof[0::2, 1] = edof[:, 1]
        new_edof[0::2, 2] = edof[:, 2]
        new_edof[1::2, 0] = edof[:, 2]
        new_edof[1::2, 1] = edof[:, 3]
        new_edof[1::2, 2] = edof[:, 0]
        return new_edof
    elif edof.shape[1] == 8:
        new_edof = np.zeros((edof.shape[0] * 6, 3), int)
        new_edof[0::6, 0] = edof[:, 0]
        new_edof[0::6, 1] = edof[:, 4]
        new_edof[0::6, 2] = edof[:, 7]
        new_edof[1::6, 0] = edof[:, 4]
        new_edof[1::6, 1] = edof[:, 1]
        new_edof[1::6, 2] = edof[:, 5]
        new_edof[2::6, 0] = edof[:, 5]
        new_edof[2::6, 1] = edof[:, 2]
        new_edof[2::6, 2] = edof[:, 6]
        new_edof[3::6, 0] = edof[:, 6]
        new_edof[3::6, 1] = edof[:, 3]
        new_edof[3::6, 2] = edof[:, 7]
        new_edof[4::6, 0] = edof[:, 4]
        new_edof[4::6, 1] = edof[:, 6]
        new_edof[4::6, 2] = edof[:, 7]
        new_edof[5::6, 0] = edof[:, 4]
        new_edof[5::6, 1] = edof[:, 5]
        new_edof[5::6, 2] = edof[:, 6]
        return new_edof
    else:
        error("Element topology not supported.")


def draw_nodal_values_contourf(
    values,
    coords,
    edof,
    levels=12,
    title=None,
    dofs_per_node=None,
    el_type=None,
    draw_elements=False,
):
    """Draws element nodal values as filled contours. Element topologies
    supported are triangles, 4-node quads and 8-node quads."""

    edof_tri = topo_to_tri(edof)

    ax = plt.gca()
    ax.set_aspect("equal")

    x, y = coords.T
    v = np.asarray(values)
    plt.tricontourf(x, y, edof_tri - 1, v.ravel(), levels)

    if draw_elements:
        if dofs_per_node != None and el_type != None:
            draw_mesh(coords, edof, dofs_per_node, el_type, color=(0.2, 0.2, 0.2))
        else:
            info("dofs_per_node and el_type must be specified to draw the mesh.")

    if title != None:
        ax.set(title=title)


def draw_nodal_values_contour(
    values,
    coords,
    edof,
    levels=12,
    title=None,
    dofs_per_node=None,
    el_type=None,
    draw_elements=False,
):
    """Draws element nodal values as filled contours. Element topologies
    supported are triangles, 4-node quads and 8-node quads."""

    edof_tri = topo_to_tri(edof)

    ax = plt.gca()
    ax.set_aspect("equal")

    x, y = coords.T
    v = np.asarray(values)
    plt.tricontour(x, y, edof_tri - 1, v.ravel(), levels)

    if draw_elements:
        if dofs_per_node != None and el_type != None:
            draw_mesh(coords, edof, dofs_per_node, el_type, color=(0.2, 0.2, 0.2))
        else:
            info("dofs_per_node and el_type must be specified to draw the mesh.")

    if title != None:
        ax.set(title=title)


def draw_nodal_values_shaded(
    values,
    coords,
    edof,
    title=None,
    dofs_per_node=None,
    el_type=None,
    draw_elements=False,
):
    """Draws element nodal values as shaded triangles. Element topologies
    supported are triangles, 4-node quads and 8-node quads."""

    edof_tri = topo_to_tri(edof)

    ax = plt.gca()
    ax.set_aspect("equal")

    x, y = coords.T
    v = np.asarray(values)
    plt.tripcolor(x, y, edof_tri - 1, v.ravel(), shading="gouraud")

    if draw_elements:
        if dofs_per_node != None and el_type != None:
            draw_mesh(coords, edof, dofs_per_node, el_type, color=(0.2, 0.2, 0.2))
        else:
            info("dofs_per_node and el_type must be specified to draw the mesh.")

    if title != None:
        ax.set(title=title)


draw_nodal_values = draw_nodal_values_contourf


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
    Draws the geometry (points and curves) in geoData
    Args:
        geoData:
            GeoData object. Geodata contains geometric information of the model.
        axes:
            Matplotlib Axes. The Axes where the model will be drawn. If unspecified the current Axes will be used, or a new Axes will be created if none exist.
        axes_adjust:
            Boolean. If True the view will be changed to show the whole model. Default True.
        draw_points:
            Boolean. If True points will be drawn.
        label_points:
            Boolean. If True Points will be labeled. The format is: ID[marker]. If a point has marker==0 only the ID is written.
        label_curves:
            Boolean. If True Curves will be labeled. The format is: ID(elementsOnCurve)[marker].
        font_size:
            Integer. Size of the text in the text labels. Default 11.
        N:
            Integer. The number of discrete points per curve segment. Default 20. Increase for smoother curves. Decrease for better performance.
        rel_margin:
            Extra spacing between geometry and axis
    """

    if axes is None:
        ax = plt.gca()
    else:
        ax = axes

    ax.set_aspect("equal")
    ax.set_frame_on(draw_axis)

    if draw_points:
        P = np.array(geometry.getPointCoords())  # M-by-3 list of M points.
        # plotArgs = {'mc':'r', 'mw':5, 'lw':0, 'ms':'o', 'axesAdjust':False, 'axes':axes}
        plotArgs = {"marker": "o", "ls": ""}
        if geometry.is3D:
            plt.plot(P[:, 0], P[:, 1], P[:, 2], **plotArgs)
        else:
            plt.plot(P[:, 0], P[:, 1], **plotArgs)

        if label_points:  # Write text label at the points:
            # [[x, y, z], elSize, marker]
            for ID, (xyz, el_size, marker) in geometry.points.items():
                text = "  " + str(ID) + ("[%s]" % marker if marker != 0 else "")
                plt.text(xyz[0], xyz[1], text, fontsize=font_size, color=(0.5, 0, 0.5))

    for ID, (
        curveName,
        pointIDs,
        marker,
        elementsOnCurve,
        _,
        _,
    ) in geometry.curves.items():
        points = geometry.getPointCoords(pointIDs)
        if curveName == "Spline":
            P = _catmullspline(points, N)
        if curveName == "BSpline":
            P = _bspline(points, N)
        if curveName == "Circle":
            P = _circleArc(*points, pointsOnCurve=N)
        if curveName == "Ellipse":
            P = _ellipseArc(*points, pointsOnCurve=N)
        # plotArgs = {'lc':'k', 'ms':None, 'axesAdjust':False, 'axes':axes} #Args for plot style. Black lines with no symbols at points.

        # Args for plot style. Black lines with no symbols at points.
        plotArgs = {"color": "black"}

        if geometry.is3D:
            plt.plot(P[:, 0], P[:, 1], P[:, 2], **plotArgs)
        else:
            plt.plot(P[:, 0], P[:, 1], **plotArgs)

        if label_curves:
            # Sort of midpoint along the curve. Where the text goes.
            midP = P[int(P.shape[0] * 7.0 / 12), :].tolist()
            # Create the text for the curve. Includes ID, elementsOnCurve, and marker:
            text = " " + str(ID)
            text += "(%s)" % (elementsOnCurve) if elementsOnCurve is not None else ""
            # Something like "4(5)[8]"
            text += "[%s]" % (marker) if marker != 0 else ""
            plt.text(midP[0], midP[1], text, fontsize=font_size)

    if title != None:
        plt.title(title)

    min_x, max_x, min_y, max_y = geometry.bounding_box_2d()

    g_width = max_x - min_x
    g_height = max_y - min_y

    if g_width > g_height:
        margin = rel_margin * g_width
    else:
        margin = rel_margin * g_height

    bottom, top = ax.get_ylim()
    left, right = ax.get_xlim()
    ax.set_ylim(bottom - margin, top + margin)
    ax.set_xlim(left - margin, right + margin)

    # if axesAdjust:
    #    _adjustaxes(axes, geoData.is3D)
    # axes.daspectAuto = False
    # axes.daspect = (1,1,1)


# drawGeometry = draw_geometry


def _catmullspline(controlPoints, pointsOnEachSegment=10):
    """
    Returns points on a Catmull-Rom spline that interpolated the control points.
    Inital/end tangents are created by mirroring the second/second-to-last)
    control points in the first/last points.

    Params:
    controlPoints - Numpy array containing the control points of the spline.
                    Each row should contain the x,y,(z) values.
                    [[x1, y2],
                     [x2, y2],
                        ...
                     [xn, yn]]

    pointsOnEachSegment - The number of points on each segment of the curve.
                        If there are n control points and k samplesPerSegment,
                        then there will be (n+1)*k numeric points on the curve.
    """
    controlPoints = np.asarray(controlPoints)  # Convert to array if input is a list.
    if (controlPoints[0, :] == controlPoints[-1, :]).all():
        # If the curve is closed we extend each opposite endpoint to the other side
        CPs = np.asmatrix(
            np.vstack((controlPoints[-2, :], controlPoints, controlPoints[1, :]))
        )
    else:  # Else make mirrored endpoints:
        CPs = np.asmatrix(
            np.vstack(
                (
                    2 * controlPoints[0, :] - controlPoints[1, :],
                    controlPoints,
                    2 * controlPoints[-1, :] - controlPoints[-2, :],
                )
            )
        )
    M = 0.5 * np.matrix([[0, 2, 0, 0], [-1, 0, 1, 0], [2, -5, 4, -1], [-1, 3, -3, 1]])
    t = np.linspace(0, 1, pointsOnEachSegment)
    T = np.matrix([[1, s, pow(s, 2), pow(s, 3)] for s in t])
    return np.asarray(
        np.vstack([T * M * CPs[j - 1 : j + 3, :] for j in range(1, len(CPs) - 2)])
    )


def _bspline(controlPoints, pointsOnCurve=20):
    """
    Uniform cubic B-spline.

    Params:
    controlPoints - Control points. Numpy array. One coordinate per row.
    pointsOnCurve - number of sub points per segment

    Mirrored start- and end-points are added if the curve is not closed.
    If the curve is closed some points are duplicated to make the closed
    spline continuous.
    (See http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-curve-closed.html)

    Based on descriptions on:
    http://www.siggraph.org/education/materials/HyperGraph/modeling/splines/b_spline.htm
    http://en.wikipedia.org/wiki/B-spline#Uniform_cubic_B-splines
    """
    controlPoints = np.asarray(controlPoints)  # Convert to array if input is a list.
    if (controlPoints[0, :] == controlPoints[-1, :]).all():
        # If the curve is closed we extend each opposite endpoint to the other side
        CPs = np.asmatrix(
            np.vstack((controlPoints[-2, :], controlPoints, controlPoints[1, :]))
        )
    else:  # Else make mirrored endpoints:
        CPs = np.asmatrix(
            np.vstack(
                (
                    2 * controlPoints[0, :] - controlPoints[1, :],
                    controlPoints,
                    2 * controlPoints[-1, :] - controlPoints[-2, :],
                )
            )
        )
    M = (1.0 / 6) * np.matrix(
        [[-1, 3, -3, 1], [3, -6, 3, 0], [-3, 0, 3, 0], [1, 4, 1, 0]]
    )
    t = np.linspace(0, 1, pointsOnCurve)
    T = np.matrix([[pow(s, 3), pow(s, 2), s, 1] for s in t])

    return np.asarray(
        np.vstack([T * M * CPs[i - 1 : i + 3, :] for i in range(1, len(CPs) - 2)])
    )


def _circleArc(start, center, end, pointsOnCurve=20):
    return _ellipseArc(start, center, start, end, pointsOnCurve)


# Update function causing the error

def _ellipseArc(start, center, majAxP, end, pointsOnCurve=20):
    """Input are 3D 1-by-3 numpy arrays or vectors"""
    # First part is to find a similarity transform in 3D that transform the ellipse to
    # the XY-plane with the center at the origin and the major axis of the ellipse along the X-axis.

    # convert to arrays in case inputs are lists:
    (
        start,
        center,
        majAxP,
        end,
    ) = (
        np.asarray(start),
        np.asarray(center),
        np.asarray(majAxP),
        np.asarray(end),
    )

    zPrim = np.cross(start - center, end - center)
    zPrim = zPrim / np.linalg.norm(zPrim)
    xPrim = (majAxP - center) / np.linalg.norm(majAxP - center)
    yPrim = np.cross(zPrim, xPrim)

    # Rotation matrix from ordinary coords to system where ellipse is in the XY-plane. (Actually hstack)
    R = np.vstack((xPrim, yPrim, zPrim)).T
    # Building Transformation matrix. -center is translation vector from ellipse center to origin.
    T = np.hstack((R, np.asarray(center).reshape(-1, 1)))
    # Transformation matrix for homogenous coordinates.
    T = np.vstack((T, [0, 0, 0, 1]))

    # Convert to arrays explicitly to avoid matrix issues
    startHC = np.vstack((np.asarray(start).reshape(-1, 1), [1]))
    # start and end points as column vectors in homogenous coordinates
    endHC = np.vstack((np.asarray(end).reshape(-1, 1), [1]))

    s = np.linalg.inv(T) @ startHC
    # start and end points in the new coordinate system
    e = np.linalg.inv(T) @ endHC

    xs, ys = s[0, 0], s[1, 0]
    # Just extract x & y from the new start and endpoints
    xe, ye = e[0, 0], e[1, 0]

    a = np.sqrt((pow(ye * xs, 2) - pow(xe * ys, 2)) / (pow(ye, 2) - pow(ys, 2)))
    b = np.sqrt(
        (pow(ye * xs, 2) - pow(xe * ys, 2))
        / (
            (pow(ye, 2) - pow(ys, 2))
            * ((pow(xe, 2) - pow(xs, 2)) / (pow(ys, 2) - pow(ye, 2)))
        )
    )

    # atan2 is a function that goes from -pi to pi. It gives the signed angle from the X-axis to point (y,x)
    ts = atan2(ys / b, xs / a)
    # We can't use the (transformed) start- and endpoints directly, but we divide x and y by the
    te = atan2(ye / b, xe / a)
    # ellipse minor&major axes to get the parameter t that corresponds to the point on the ellipse.
    # See ellipse formula: x = a * cos (t), y = b * sin (t).
    # So ts and te are the parameter values of the start- and endpoints (in the transformed coordinate system).

    if ts > te:
        # swap if the start point comes before the endpoint in the parametric parameter that goes around the ellipse.
        ts, te = te, ts
    if te - ts < np.pi:
        # parameter of ellipse curve. NOT angle to point on curve (like it could be for a circle).
        times = np.linspace(ts, te, pointsOnCurve)
    # the shortest parameter distance between start- and end-point stradles the discontinuity that jumps from pi to -pi.
    else:
        # number of points on the first length.
        ps1 = round(pointsOnCurve * (pi - te) / (2 * pi - te + ts))
        # number of points on the first length.
        ps2 = round(pointsOnCurve * (ts + pi) / (2 * pi - te + ts))
        times = np.concatenate((np.linspace(te, pi, ps1), np.linspace(-pi, ts, ps2)))

    ellArc = np.array(
        [[a * cos(t), b * sin(t)] for t in times]
    ).T  # points on arc (in 2D)
    # Make 3D homogenous coords by adding rows of 0s and 1s.
    zeros_row = np.zeros((1, ellArc.shape[1]))
    ones_row = np.ones((1, ellArc.shape[1]))
    ellArc = np.vstack((ellArc, zeros_row, ones_row))
    ellArc = T @ ellArc  # Transform back to the original coordinate system
    return np.asarray(ellArc.T[:, 0:3])  # return points as an N-by-3 array.

def eldraw2(ex, ey, plotpar=[1, 2, 1], elnum=[]):
    """
    eldraw2(ex,ey,plotpar,elnum)
    eldraw2(ex,ey,plotpar)
    eldraw2(ex,ey)

     PURPOSE
       Draw the undeformed 2D mesh for a number of elements of
       the same type. Supported elements are:

       1) -> bar element              2) -> beam el.
       3) -> triangular 3 node el.    4) -> quadrilateral 4 node el.
       5) -> 8-node isopar. elemen

     INPUT
        ex,ey:.......... nen:   number of element nodes
                         nel:   number of elements
        plotpar=[ linetype, linecolor, nodemark]

                 linetype=1 -> solid    linecolor=1 -> black
                          2 -> dashed             2 -> blue
                          3 -> dotted             3 -> magenta
                                                  4 -> red

                 nodemark=1 -> circle
                          2 -> star
                          0 -> no mark

        elnum=edof(:,1) ; i.e. the first column in the topology matrix

        Rem. Default is solid white lines with circles at nodes.
    """

    if ex.shape == ey.shape:
        if ex.ndim != 1:
            nen = ex.shape[1]
        else:
            nen = ex.shape[0]

            ex = ex.reshape(1, nen)
            ey = ey.reshape(1, nen)
    else:
        raise ValueError("Check size of ex, ey dimensions.")

    line_type = plotpar[0]
    line_color = plotpar[1]
    node_mark = plotpar[2]

    # Translate CALFEM plotpar to visvis

    if line_type == 1:
        mpl_line_style = "solid"
    elif line_type == 2:
        mpl_line_style = (0, (5, 5))
    elif line_type == 3:
        mpl_line_style = "dotted"

    if line_color == 1:
        mpl_line_color = (0, 0, 0)  # 'k'
    elif line_color == 2:
        mpl_line_color = (0, 0, 1)  # 'b'
    elif line_color == 3:
        mpl_line_color = (1, 0, 1)  # 'm'
    elif line_color == 4:
        mpl_line_color = (1, 0, 0)  # 'r'

    if node_mark == 1:
        mpl_node_mark = "o"
    elif node_mark == 2:
        mpl_node_mark = "x"
    elif node_mark == 0:
        mpl_node_mark = ""

    plt.axis("equal")

    draw_element_numbers = False

    if len(elnum) == ex.shape[0]:
        draw_element_numbers = True

    draw_elements(ex, ey, color=mpl_line_color, line_style=mpl_line_style, filled=False)
    if mpl_node_mark != "":
        draw_node_circles(
            ex, ey, color=mpl_line_color, filled=False, marker_type=mpl_node_mark
        )

    return None


def scalfact2(ex, ey, ed, rat=0.2):
    """
    [sfac]=scalfact2(ex,ey,ed,rat)
    [sfac]=scalfact2(ex,ey,ed)
    -------------------------------------------------------------
    PURPOSE
    Determine scale factor for drawing computational results, such as
    displacements, section forces or flux.

    INPUT
        ex,ey:  element node coordinates

        ed:     element displacement matrix or section force matrix

        rat: relation between illustrated quantity and element size.
        If not specified, 0.2 is used.

    -------------------------------------------------------------

    LAST MODIFIED: O Dahlblom  2004-09-15
                   J Lindemann 2021-12-29 (Python)

    Copyright (c)  Division of Structural Mechanics and
                   Division of Solid Mechanics.
                   Lund University
    -------------------------------------------------------------
    """

    #  if ~((nargin==3)|(nargin==4))
    #     disp('??? Wrong number of input arguments!')
    #     return
    #  end

    if ex.shape == ey.shape:
        if ex.ndim != 1:
            nen = ex.shape[1]
        else:
            nen = ex.shape[0]
    else:
        raise ValueError("Check size of ex, ey dimensions.")

    dx_max = float(np.max(ex)) - float(np.min(ex))
    dy_max = float(np.max(ey)) - float(np.min(ey))
    dl_max = max(dx_max, dy_max)
    ed_max = float(np.max(np.max(np.abs(ed))))

    # dxmax=max(max(ex')-min(ex')); dymax=max(max(ey')-min(ey'));
    # dlmax=max(dxmax,dymax);
    # edmax=max(max(abs(ed)));

    k = rat

    return k * dl_max / ed_max


def eliso2_mpl(ex, ey, ed):
    plt.axis("equal")

    gx = []
    gy = []
    gz = []

    for elx, ely, scl in zip(ex, ey, ed):
        for x in elx:
            gx.append(x)
        for y in ely:
            gy.append(y)
        for z in ely:
            gz.append(y)

    plt.tricontour(gx, gy, gz, 5)


def pltstyle(plotpar):
    """
    -------------------------------------------------------------
     PURPOSE
       Define define linetype,linecolor and markertype character codes.

     INPUT
        plotpar=[ linetype, linecolor, nodemark ]

                 linetype=1 -> solid    linecolor=1 -> black
                          2 -> dashed             2 -> blue
                          3 -> dotted             3 -> magenta
                                                  4 -> red

                 nodemark=1 -> circle
                          2 -> star
                          0 -> no mark
     OUTPUT
         s1: linetype and color for mesh lines
         s2: type and color for node markers
    -------------------------------------------------------------

     LAST MODIFIED: Ola Dahlblom 2004-09-15
     Copyright (c)  Division of Structural Mechanics and
                    Division of Solid Mechanics.
                    Lund University
    -------------------------------------------------------------
    """
    if type(plotpar) != list:
        raise TypeError("plotpar should be a list.")
    if len(plotpar) != 3:
        raise ValueError("plotpar needs to be a list of 3 values.")

    p1, p2, p3 = plotpar

    s1 = ""
    s2 = ""

    if p1 == 1:
        s1 += "-"
    elif p1 == 2:
        s1 += "--"
    elif p1 == 3:
        s1 += ":"
    else:
        raise ValueError("Invalid value for plotpar[0].")

    if p2 == 1:
        s1 += "k"
    elif p2 == 2:
        s1 += "b"
    elif p2 == 3:
        s1 += "m"
    elif p2 == 4:
        s1 += "r"
    else:
        raise ValueError("Invalid value for plotpar[1].")

    if p3 == 1:
        s2 = "ko"
    elif p3 == 2:
        s2 = "k*"
    elif p3 == 3:
        s2 = "k."
    else:
        raise ValueError("Invalid value for plotpar[2].")

    return s1, s2


def pltstyle2(plotpar):
    """
    -------------------------------------------------------------
     PURPOSE
       Define define linetype,linecolor and markertype character codes.

     INPUT
        plotpar=[ linetype, linecolor, nodemark ]

                 linetype=1 -> solid    linecolor=1 -> black
                          2 -> dashed             2 -> blue
                          3 -> dotted             3 -> magenta
                                                  4 -> red

                 nodemark=1 -> circle
                          2 -> star
                          0 -> no mark
     OUTPUT
         s1: linetype and color for mesh lines
         s2: type and color for node markers
    -------------------------------------------------------------

     LAST MODIFIED: Ola Dahlblom 2004-09-15
     Copyright (c)  Division of Structural Mechanics and
                    Division of Solid Mechanics.
                    Lund University
    -------------------------------------------------------------
    """

    cfc.check_list_array(plotpar, "plotpar needs to be a list or an array of 3 values.")
    cfc.check_length(plotpar, 3, "plotpar needs to contain 3 values.")

    p1, p2, p3 = plotpar

    s1 = ""
    s2 = ""

    line_style = ""
    line_color = ""
    node_color = ""
    node_type = ""

    if p1 == 1:
        line_style = "solid"
    elif p1 == 2:
        line_style = (0, (5, 5))
    elif p1 == 3:
        line_style = "dotted"
    else:
        raise ValueError("Invalid value for plotpar[0].")

    if p2 == 1:
        line_color = (0, 0, 0)
    elif p2 == 2:
        line_color = (0, 0, 1)
    elif p2 == 3:
        line_color = (1, 0, 1)
    elif p2 == 4:
        line_color = (1, 0, 0)
    else:
        raise ValueError("Invalid value for plotpar[1].")

    if p3 == 1:
        node_color = (0, 0, 0)
        node_type = "o"
    elif p3 == 2:
        node_color = (0, 0, 0)
        node_type = "*"
    elif p3 == 3:
        node_color = (0, 0, 0)
        node_type = "."
    else:
        raise ValueError("Invalid value for plotpar[2].")

    return line_color, line_style, node_color, node_type


def eldisp2(ex, ey, ed, plotpar=[2, 1, 1], sfac=None):
    """
    eldisp2(ex,ey,ed,plotpar,sfac)
    [sfac]=eldisp2(ex,ey,ed,plotpar)
    [sfac]=eldisp2(ex,ey,ed)
    -------------------------------------------------------------
     PURPOSE
       Draw the deformed 2D mesh for a number of elements of
       the same type. Supported elements are:

               1) -> bar element              2) -> beam el.
               3) -> triangular 3 node el.    4) -> quadrilateral 4 node el.
               5) -> 8-node isopar. element
      INPUT
        ex,ey:.......... nen:   number of element nodes
                         nel:   number of elements
        ed:     element displacement matrix

        plotpar=[  linetype, linecolor, nodemark]

                 linetype=1 -> solid    linecolor=1 -> black
                          2 -> dashed             2 -> blue
                          3 -> dotted             3 -> magenta
                                                  4 -> red
                 nodemark=1 -> circle
                          2 -> star
                          0 -> no mark

        sfac:  scale factor for displacements

        Rem. Default if sfac and plotpar is left out is auto magnification
             and dashed black lines with circles at nodes -> plotpar=[2 1 1]
    -------------------------------------------------------------

     LAST MODIFIED: O Dahlblom 2004-10-01
                    J Lindemann 2021-12-30 (Python)

     Copyright (c)  Division of Structural Mechanics and
                    Division of Solid Mechanics.
                    Lund University
    -------------------------------------------------------------
    """

    if ex.shape == ey.shape:
        if ex.ndim != 1:
            nen = ex.shape[1]

            if ed.shape[0] != ex.shape[0]:
                raise ValueError("Check size of ed/ex dimensions.")

            ned = ed.shape[1]
        else:
            nen = ex.shape[0]
            ned = ed.shape[0]

            ex = ex.reshape(1, nen)
            ey = ey.reshape(1, nen)
            ed = ed.reshape(1, ned)
    else:
        raise ValueError("Check size of ex, ey dimensions.")

    dx_max = float(np.max(ex)) - float(np.min(ex))
    dy_max = float(np.max(ey)) - float(np.min(ey))
    dl_max = max(dx_max, dy_max)
    ed_max = float(np.max(np.max(np.abs(ed))))
    krel = 0.1

    if sfac is None:
        sfac = krel * dl_max / ed_max

    k = sfac

    line_color, line_style, node_color, node_style = pltstyle2(plotpar)

    if nen == 2:
        if ned == 4:
            x = np.transpose(ex + k * ed[:, [0, 2]])
            y = np.transpose(ey + k * ed[:, [1, 3]])
            xc = np.transpose(x)
            yc = np.transpose(y)
        elif ned == 6:
            x = np.transpose(ex + k * ed[:, [0, 3]])
            y = np.transpose(ey + k * ed[:, [1, 4]])
            exc, eyc = beam2crd(ex, ey, ed, k)
            xc = exc
            yc = eyc
    elif nen == 3:
        pass
    elif nen == 4:
        pass
    elif nen == 8:
        pass
    else:
        print("Error: Element type is not supported.")
        return

    draw_elements(
        xc, yc, color=line_color, line_style=line_style, filled=False, closed=False
    )

    if node_style != "":
        draw_node_circles(x, y, color=node_color, filled=False, marker_type=node_style)


# % ********** Bar or Beam elements *************
#     if nen==2
#        if ned==4  % -----------  Bar elements -------------
#           x=(ex+k*ed(:,[1 3]))';
#           y=(ey+k*ed(:,[2 4]))';
#           xc=x;
#           yc=y;
#        elseif ned==6  % -------- Beam elements ------------
#           x=(ex+k*ed(:,[1 4]))';
#           y=(ey+k*ed(:,[2 5]))';
#           [exc,eyc]=beam2crd(ex,ey,ed,k);
#           xc=exc';
#           yc=eyc';
#        end
# % ********* 2D triangular elements ************
#     elseif nen==3
#        x=(ex+k*ed(:,[1 3 5]))';
#        y=(ey+k*ed(:,[2 4 6]))';
#        xc=[x; x(1,:)];
#        yc=[y; y(1,:)];

# % ********* 2D quadrilateral elements *********
#     elseif nen==4
#        x=(ex+k*ed(:,[1 3 5 7]))';
#        y=(ey+k*ed(:,[2 4 6 8]))';
#        xc=[x; x(1,:)];
#        yc=[y; y(1,:)];
#  % ********* 2D 8-node quadratic elements *********
#     elseif nen==8
#        x=(ex+k*ed(:,[1 3 5 7 9 11 13 15]));
#        y=(ey+k*ed(:,[2 4 6 8 10 12 14 16]));
# %       xc=[x(1); x(5); x(2); x(6); x(3); x(7); x(4); x(8);x(1)];
# %       yc=[y(1); y(5); y(2); y(6); y(3); y(7); y(4); y(8);y(1)];
# %
# % isoparametric elements
# %
#     t=-1;
#     n=0;
#     for s=-1:0.4:1
#       n=n+1;
#       N1=-1/4*(1-t)*(1-s)*(1+t+s);
#       N2=-1/4*(1+t)*(1-s)*(1-t+s);
#       N3=-1/4*(1+t)*(1+s)*(1-t-s);
#       N4=-1/4*(1-t)*(1+s)*(1+t-s);
#       N5=1/2*(1-t*t)*(1-s);
#       N6=1/2*(1+t)*(1-s*s);
#       N7=1/2*(1-t*t)*(1+s);
#       N8=1/2*(1-t)*(1-s*s);
#       N=[ N1, N2, N3 ,N4, N5, N6, N7, N8 ];

#       x1(n,:)=N*x';
#       y1(n,:)=N*y';
#     end;
#     xc=[xc x1];
#     yc=[yc y1];
#     clear x1
#     clear y1
# %
#     s=1;
#     n=0;
#     for t=-1:0.4:1
#       n=n+1;
#       N1=-1/4*(1-t)*(1-s)*(1+t+s);
#       N2=-1/4*(1+t)*(1-s)*(1-t+s);
#       N3=-1/4*(1+t)*(1+s)*(1-t-s);
#       N4=-1/4*(1-t)*(1+s)*(1+t-s);
#       N5=1/2*(1-t*t)*(1-s);
#       N6=1/2*(1+t)*(1-s*s);
#       N7=1/2*(1-t*t)*(1+s);
#       N8=1/2*(1-t)*(1-s*s);
#       N=[ N1, N2, N3 ,N4, N5, N6, N7, N8 ];

#       x1(n,:)=N*x';
#       y1(n,:)=N*y';
#     end;
#     xc=[xc x1];
#     yc=[yc y1];
#     clear x1
#     clear y1
# %
#     t=1;
#     n=0;
#     for s=1:-0.4:-1
#       n=n+1;
#       N1=-1/4*(1-t)*(1-s)*(1+t+s);
#       N2=-1/4*(1+t)*(1-s)*(1-t+s);
#       N3=-1/4*(1+t)*(1+s)*(1-t-s);
#       N4=-1/4*(1-t)*(1+s)*(1+t-s);
#       N5=1/2*(1-t*t)*(1-s);
#       N6=1/2*(1+t)*(1-s*s);
#       N7=1/2*(1-t*t)*(1+s);
#       N8=1/2*(1-t)*(1-s*s);
#       N=[ N1, N2, N3 ,N4, N5, N6, N7, N8 ];

#       x1(n,:)=N*x';
#       y1(n,:)=N*y';
#     end;
#     xc=[xc x1];
#     yc=[yc y1];
#     clear x1
#     clear y1
# %
#     s=-1;
#     n=0;
#     for t=1:-0.4:-1
#       n=n+1;
#       N1=-1/4*(1-t)*(1-s)*(1+t+s);
#       N2=-1/4*(1+t)*(1-s)*(1-t+s);
#       N3=-1/4*(1+t)*(1+s)*(1-t-s);
#       N4=-1/4*(1-t)*(1+s)*(1+t-s);
#       N5=1/2*(1-t*t)*(1-s);
#       N6=1/2*(1+t)*(1-s*s);
#       N7=1/2*(1-t*t)*(1+s);
#       N8=1/2*(1-t)*(1-s*s);
#       N=[ N1, N2, N3 ,N4, N5, N6, N7, N8 ];

#       x1(n,:)=N*x';
#       y1(n,:)=N*y';
#     end;
#     xc=[xc x1];
#     yc=[yc y1];
#     clear x1
#     clear y1
# %
# %**********************************************************
#     else
#        disp('Sorry, this element is currently not supported!')
#        return
#     end
# % ************* plot commands *******************
#     hold on
#     axis equal
#     plot(xc,yc,s1)
#     plot(x,y,s2)
#     hold off
# %--------------------------end--------------------------------


def dispbeam2(ex, ey, edi, plotpar=[2, 1, 1], sfac=None):
    """
        dispbeam2(ex,ey,edi,plotpar,sfac)
        [sfac]=dispbeam2(ex,ey,edi)
        [sfac]=dispbeam2(ex,ey,edi,plotpar)
    ------------------------------------------------------------------------
        PURPOSE
        Draw the displacement diagram for a two dimensional beam element.

        INPUT:   ex = [ x1 x2 ]
                ey = [ y1 y2 ]	element node coordinates.

                edi = [ u1 v1;
                       u2 v2;
                                 .....] 	matrix containing the displacements
                                              in Nbr evaluation points along the beam.

                plotpar=[linetype, linecolour, nodemark]

                         linetype=1 -> solid   linecolour=1 -> black
                                  2 -> dashed             2 -> blue
                                  3 -> dotted             3 -> magenta
                                                         4 -> red
                         nodemark=0 -> no mark
                                  1 -> circle
                                  2 -> star
                                  3 -> point

                         sfac = [scalar] scale factor for displacements.

                Rem. Default if sfac and plotpar is left out is auto magnification
               and dashed black lines with circles at nodes -> plotpar=[1 1 1]
    ------------------------------------------------------------------------

        LAST MODIFIED: O Dahlblom  2015-11-18
                       O Dahlblom  2023-01-31 (Python)

        Copyright (c)  Division of Structural Mechanics and
                       Division of Solid Mechanics.
                       Lund University
    ------------------------------------------------------------------------
    """
    if ex.shape != ey.shape:
        raise ValueError("Check size of ex, ey dimensions.")

    rows, cols = edi.shape
    if cols != 2:
        raise ValueError("Check size of edi dimension.")
    Nbr = rows

    x1, x2 = ex
    y1, y2 = ey
    dx = x2 - x1
    dy = y2 - y1
    L = np.sqrt(dx * dx + dy * dy)
    nxX = dx / L
    nyX = dy / L
    n = np.array([nxX, nyX])

    line_color, line_style, node_color, node_style = pltstyle2(plotpar)

    if sfac is None:
        sfac = (0.1 * L) / (np.max(abs(edi)))

    eci = np.arange(0.0, L + L / (Nbr - 1), L / (Nbr - 1)).reshape(Nbr, 1)

    edi1 = edi * sfac
    # From local x-coordinates to global coordinates of the beam element.
    A = np.zeros(2 * Nbr).reshape(Nbr, 2)
    A[0, 0] = ex[0]
    A[0, 1] = ey[0]
    for i in range(1, Nbr):
        A[i, 0] = A[0, 0] + eci[i] * n[0]
        A[i, 1] = A[0, 1] + eci[i] * n[1]

    for i in range(0, Nbr):
        A[i, 0] = A[i, 0] + edi1[i, 0] * n[0] - edi1[i, 1] * n[1]
        A[i, 1] = A[i, 1] + edi1[i, 0] * n[1] + edi1[i, 1] * n[0]
    xc = np.array(A[:, 0])
    yc = np.array(A[:, 1])

    plt.plot(xc, yc, color=line_color, linewidth=1)

    A1 = np.array([A[0, 0], A[Nbr - 1, 0]]).reshape(1, 2)
    A2 = np.array([A[0, 1], A[Nbr - 1, 1]]).reshape(1, 2)
    draw_node_circles(A1, A2, color=node_color, filled=False, marker_type=node_style)


def secforce2(ex, ey, es, plotpar=[2, 1], sfac=None, eci=None):
    """
    secforce2(ex,ey,es,plotpar,sfac)
    secforce2(ex,ey,es,plotpar,sfac,eci)
    [sfac]=secforce2(ex,ey,es)
    [sfac]=secforce2(ex,ey,es,plotpar)
    --------------------------------------------------------------------------
    PURPOSE:
    Draw section force diagram for a two dimensional bar or beam element.

    INPUT:  ex = [ x1 x2 ]
                ey = [ y1 y2 ]	element node coordinates.

                es = [ S1;
                   S2;
                        ... ] 	vector containing the section force
                                        in Nbr evaluation points along the element.

            plotpar=[linecolour, elementcolour]

                linecolour=1 -> black      elementcolour=1 -> black
                           2 -> blue                     2 -> blue
                           3 -> magenta                  3 -> magenta
                           4 -> red                       4 -> red

                sfac = [scalar]	scale factor for section force diagrams.

            eci = [  x1;
                     x2;
                   ... ]  local x-coordinates of the evaluation points (Nbr).
                          If not given, the evaluation points are assumed to be uniformly
                          distributed
    --------------------------------------------------------------------------

    LAST MODIFIED: O Dahlblom  2019-12-16
                   O Dahlblom  2023-01-31 (Python)

    Copyright (c)  Division of Structural Mechanics and
                   Division of Solid Mechanics.
                   Lund University
    --------------------------------------------------------------------------
    """
    if ex.shape != ey.shape:
        raise ValueError("Check size of ex, ey dimensions.")

    c = len(es)
    Nbr = c

    x1, x2 = ex
    y1, y2 = ey
    dx = x2 - x1
    dy = y2 - y1
    L = np.sqrt(dx * dx + dy * dy)
    nxX = dx / L
    nyX = dy / L
    n = np.array([nxX, nyX])

    if sfac is None:
        sfac = (0.2 * L) / max(abs(es))

    if eci is None:
        eci = np.arange(0.0, L + L / (Nbr - 1), L / (Nbr - 1)).reshape(Nbr, 1)

    p1 = plotpar[0]
    if p1 == 1:
        line_color = (0, 0, 0)
    elif p1 == 2:
        line_color = (0, 0, 1)
    elif p1 == 3:
        line_color = (1, 0, 1)
    elif p1 == 4:
        line_color = (1, 0, 0)
    else:
        raise ValueError("Invalid value for plotpar[1].")
    line_style = "solid"

    p2 = plotpar[1]
    if p2 == 1:
        line_color1 = (0, 0, 0)
    elif p2 == 2:
        line_color1 = (0, 0, 1)
    elif p2 == 3:
        line_color1 = (1, 0, 1)
    elif p2 == 4:
        line_color1 = (1, 0, 0)
    else:
        raise ValueError("Invalid value for plotpar[1].")

    a = len(eci)
    if a != c:
        raise ValueError("Check size of eci dimension.")

    es = es * sfac

    # From local x-coordinates to global coordinates of the element
    A = np.zeros(2 * Nbr).reshape(Nbr, 2)
    A[0, 0] = ex[0]
    A[0, 1] = ey[0]
    for i in range(Nbr):
        A[i, 0] = A[0, 0] + eci[i] * n[0]
        A[i, 1] = A[0, 1] + eci[i] * n[1]

    B = np.array(A)

    # Plot diagram
    for i in range(0, Nbr):
        A[i, 0] = A[i, 0] + es[i] * n[1]
        A[i, 1] = A[i, 1] - es[i] * n[0]

    xc = np.array(A[:, 0])
    yc = np.array(A[:, 1])

    plt.plot(xc, yc, color=line_color, linewidth=1)

    # Plot stripes in diagram
    xs = np.zeros(2)
    ys = np.zeros(2)
    for i in range(Nbr):
        xs[0] = B[i, 0]
        xs[1] = A[i, 0]
        ys[0] = B[i, 1]
        ys[1] = A[i, 1]
        plt.plot(xs, ys, color=line_color, linewidth=1)

    # Plot element
    plt.plot(ex, ey, color=line_color1, linewidth=2)


def scalgraph2(sfac, magnitude, plotpar=2):
    """
    scalgraph2(sfac, magnitude, plotpar)
    scalgraph2(sfac, magnitude)
    -------------------------------------------------------------
    PURPOSE
    Draw a graphic scale

    INPUT:  sfac = [scalar]	scale factor.

            magnitude = [Ref x y]	The graphic scale has a length equivalent
                to Ref and starts at coordinates (x,y).
                If no coordinates are given the starting
            point will be (0,-0.5).

            plotpar=[linecolor]
                linecolor=1 -> black
                2 -> blue
                3 -> magenta
                4 -> red
    -------------------------------------------------------------

     LAST MODIFIED: O Dahlblom  2015-12-02
                    O Dahlblom  2023-01-23 (Python)

     Copyright (c)  Division of Structural Mechanics and
                    Division of Solid Mechanics.
                    Lund University
    -------------------------------------------------------------
    """
    cols = len(magnitude)
    if cols != 1 and cols != 3:
        raise ValueError("Check size of magnitude input argument.")
    if cols == 1:
        N = magnitude
        x = 0
        y = -0.5
    if cols == 3:
        N, x, y = magnitude

    L = N * sfac

    if plotpar == 1:
        line_color = (0, 0, 0)
    elif plotpar == 2:
        line_color = (0, 0, 1)
    elif plotpar == 3:
        line_color = (1, 0, 1)
    elif plotpar == 4:
        line_color = (1, 0, 0)
    else:
        raise ValueError("Invalid value for plotpar[1].")

    plt.plot([x, (x + L)], [y, y], color=line_color, linewidth=1)
    plt.plot([x, x], [(y - L / 20), (y + L / 20)], color=line_color, linewidth=1)
    plt.plot(
        [(x + L), (x + L)], [(y - L / 20), (y + L / 20)], color=line_color, linewidth=1
    )
    plt.text(x + L * 1.1, (y - L / 20), str(N))
