# -*- coding: utf-8 -*-

import visvis as vv

from visvis.wibjects.colorWibjects import Colorbar
from visvis import Colormapable
from visvis.wobjects.textures import minmax

import numpy as np
import matplotlib.pyplot as plt

import OpenGL.GL as gl #@UnresolvedImport

from numpy import sin, cos, pi
from math import atan2

import logging as cflog

globalWindows = [] # For supporting ElementView:s for eldraw ...

global globalVisVisApp
global visApp
visApp = vv.use('qt5') # use qt4


def error(msg):
    """Log error message"""
    cflog.error(msg)

def info(msg):
    """Log information message"""
    cflog.info(msg)

def figure_class():
    """Return visvis Figure class."""
    global visApp
    return visApp.GetFigureClass()

figureClass = figure_class

def figure(figure=None, show=True):
    """Create a visvis figure with extras."""
    f = None
    
    if figure == None:
        f = vv.figure()
    else:
        try:
            f = vv.figure(figure)
        except:
            f = vv.figure() 

    if show:
        f._widget.show()
        f._widget.raise_()
        
    return f

def close_all():
    """Close all visvis windows."""
    vv.closeAll()
    
closeAll = close_all
    
def clf():
    """Clear visvis figure"""
    vv.clf()
    
def gca():
    """Get current axis of the current visvis figure."""
    return vv.gca()
    
def subplot(*args):
    """Create a visvis subplot."""
    return vv.subplot(*args)
    
def camera3d():
    """Get visvis 3D camera."""
    return vv.cameras.ThreeDCamera()

def show_grid(flag = True):
    """Show grid."""
    vv.gca().axis.showGrid = flag    
    
showGrid = show_grid
        
def show_and_wait():
    """Show visvis windows and enter application loop."""
    global globalVisVisApp
    globalVisVisApp = vv.use()
    globalVisVisApp.Create()
    globalVisVisApp.Run()
    
showAndWait = show_and_wait
    
def show_and_wait_mpl():
    plt.show()
    
showAndWaitMpl = show_and_wait_mpl

def get_color_bar(axes=None):
    '''
    Returns the Colorbar.
    If axes is None the colorbar in the current axes will be found.
    If several colorbars exists in the axes the first found will be returned
    If no colorbar is found None is returned.
    '''
    #An ugly solution, but visvis seems to have no other way of getting the colorbar,
    #or most other entities that exist in the axes. 
    if axes is None:
        axes = vv.gca()
    for obj in axes.children: 
        if type(obj) == Colorbar:
            return obj  
    return None

getColorbar = get_color_bar
    
def color_bar(axes=None):
    """Short form of getColorbar"""
    return getColorbar(axes)

colorBar = color_bar
    

def _makeColorBar(text, axes=None):
    '''
    Convenience function that finds the current colorbar in the axes
    or creates a new one if one does not exist.
    The reason is that colorbars can not be deleted without clearing 
    the whole figure, and several colorbars can exist simultaneously. 
    This should be avoided.
    '''
    if axes is None:
        axes = vv.gca()
    colBar = getColorbar(axes)
    if colBar is None:
        vv.colorbar(axes).SetLabel(text) #Creates a colorbar and sets the label.
    else:
        colBar.SetLabel(text)# A colorbar already exists, Change label.

def add_label(text, pos, angle=0, font_name=None, font_size=9, color='k', bgcolor=None, axes=None):
    '''
    Adds a label inside the axes. Returns the Label object.
    Parameters:
    text    - String. The text of the label
    pos     - Tuple with two numbers. The (x,y) position of the label with origin
              at the upper left corner.
    angle   - Float or int. The rotation of the label in degrees.
    font_name- String. Either 'mono', 'sans' or 'serif'.
    font_size- Int. Size of the text. Default 9.
    color   - A 3-tuple or a character in 'rgbycmkw', etc that defines text color.
              Default 'k' (black).
    bgcolor - Background color. See color. Default None.
    axes    - Axes wherein the label is placed. If None then the current axes is
              chosen.
    '''
    if axes is None:
        axes = vv.gca()
    label = vv.Label(axes, text, font_name, font_size, color)
    label.position = pos
    label.bgcolor = bgcolor
    label.textAngle = angle
    return label

addLabel = add_label

def label(text, pos, angle=0, font_name=None, font_size=9, color='k', bgcolor=None, axes=None):
    return addLabel(text, pos, angle, font_name, font_size, color, bgcolor, axes)    

def add_text(text, pos, angle=0, font_name=None, font_size=9, color='k', bgcolor=None, axes=None):
    '''
    Adds a text in the world space. Returns the Text object.
    Parameters:
    text    - String. The text of the label
    pos     - Tuple with two or three numbers. The (x,y,z) position of the text in
              world space.
    angle   - Float or int. The rotation of the label in degrees.
    font_name- String. Either 'mono', 'sans' or 'serif'.
    font_size- Int. Size of the text. Default 9.
    color   - A 3-tuple or a character in 'rgbycmkw', etc that defines text color.
              Default 'k' (black).
    bgcolor - Background color. See color. Default None.
    axes    - Axes wherein the label is placed. If None then the current axes is 
              chosen.
    '''
    if axes is None:
        axes = vv.gca()
    text = vv.Text(axes, text, *pos, fontName=font_name, fontSize=font_size, color=color)
    text.bgcolor = bgcolor
    text.textAngle = angle
    return text

addText = add_text

def text(txt, pos, angle=0, font_name=None, font_size=9, color='k', bgcolor=None, axes=None):
    return addText(txt, pos, angle, font_name, font_size, color, bgcolor, axes)   

def draw_mesh(coords, edof, dofs_per_node, el_type, axes=None, axes_adjust=True, 
             title=None, color=(0,0,0), face_color=(1,1,1), filled=False):
    '''
    Draws wire mesh of model in 2D or 3D. Returns the Mesh object that represents
    the mesh.
    Parameters:
    coords      - An N-by-2 or N-by-3 array. Row i contains the x,y,z coordinates
                  of node i.
    edof        - An E-by-L array. Element topology. (E is the number of elements
                  and L is the number of dofs per element)
    dofs_per_node - Integer. Dofs per node.
    el_type      - Integer. Element Type. See Gmsh manual for details. Usually 2
                  for triangles or 3 for quadrangles.
    axes        - Visvis Axes. The Axes where the model will be drawn. 
                  If unspecified the current Axes will be used, or a new Axes will
                  be created if none exist.
    axes_adjust  - Boolean. True if the view should be changed to show the whole
                  model. Default True.
    title       - String. Changes title of the figure. Default "Mesh".
    color       - 3-tuple or char. Color of the wire. Defaults to black (0,0,0).
                  Can also be given as a character in 'rgbycmkw'.
    face_color   - 3-tuple or char. Color of the faces. Defaults to white (1,1,1).
                  Parameter filled must be True or faces will not be drawn at all.
    filled      - Boolean. Faces will be drawn if True. Otherwise only the wire is
                  drawn. Default False.
    '''
    #Prep:
    axes, verts, faces, verticesPerFace, is3D = _preMeshDrawPrep(axes, coords, edof, dofs_per_node, el_type)
    #Create mesh:
    m = vv.Mesh(parent=axes, vertices=verts, faces=faces, values=color, verticesPerFace=verticesPerFace)
    #Settings:
    fShade = 'plain' if filled else None
    m.faceShading, m.edgeShading = (fShade, 'plain')
    m.edgeColor = color
    m.face_color = face_color
    m.specular = 0 
    #Adjust axes:
    if axes_adjust:
        _adjustaxes(axes, is3D)
    #Set title and return:
    vv.title(title, axes)
    return m
    
drawMesh = draw_mesh
    
def draw_nodal_values(node_vals, coords, edof, dofs_per_node, el_type, clim=None, axes=None, axes_adjust=True, draw_elements=True, title=None):
    '''
    Draws scalar nodal values in 2D or 3D. Returns the Mesh object that represents
    the mesh.
    Parameters:
    node_vals    - An N-by-1 array or a list of scalars. The Scalar values at the
                  nodes. node_vals[i] should be the value of node i 
    coords      - An N-by-2 or N-by-3 array. Row i contains the x,y,z coordinates
                  of node i.
    edof        - An E-by-L array. Element topology. (E is the number of elements 
                  and L is the number of dofs per element)
    dofs_per_node - Integer. Dofs per node.
    el_type      - Integer. Element Type. See Gmsh manual for details. Usually 2 
                  for triangles or 3 for quadrangles.
    clim        - 2-tuple. Colorbar limits (min, max). Defines the value range of
                  the colorbar. Defaults to None, in which case min/max are set to
                  min/max of node_vals.
    axes        - Visvis Axes. The Axes where the model will be drawn. 
                  If unspecified the current Axes will be used, or a new Axes will
                  be created if none exist.
    axes_adjust  - Boolean. True if the view should be changed to show the whole 
                  model. Default True.
    draw_elements  - Boolean. True if mesh wire should be drawn. Default True.
    title       - String. Changes title of the figure. Default "Node Values".
    '''    
    axes, verts, faces, verticesPerFace, is3D = _preMeshDrawPrep(axes, coords, edof, dofs_per_node, el_type)
    m = vv.Mesh(parent=axes, vertices=verts, faces=faces, values=node_vals, verticesPerFace=verticesPerFace)
    
    if clim != None: #Set colorbar limits.
        m.clim = clim
        setClim = False
    else:
        setClim = True
    
    edgeSh = 'plain' if draw_elements else None
    m.faceShading, m.edgeShading = ('smooth', edgeSh)#NOTE: It seems colormap coloring breaks when faceshading='plain'. 'smooth' must be used.
    m.ambient = 1
    m.diffuse = 0
    m.specular = 0 #Disable specular. 
    m.SetValues(node_vals, setClim) #Set the values again, because it doesn't work in the constructor for unknown reasons
    
    axes.light0.ambient = 1.0
    axes.light0.diffuse = 0.0  #Only ambient light to avoid shadows
    
    m.colormap = vv.colormaps['jet']
    _makeColorBar("Node values", axes)
        
    # Adjust axes:
    if axes_adjust:
        _adjustaxes(axes, is3D)
    
    vv.title(title, axes)
    return m

drawNodalValues = draw_nodal_values
    
def draw_element_values(ev, coords, edof, dofs_per_node, el_type, displacements=None, clim=None, axes=None, 
                      axes_adjust=True, draw_elements=True, draw_undisplaced_mesh=False, magnfac=1.0, title=None):
    '''
    Draws scalar element values in 2D or 3D. Returns the world object 
    elementsWobject that represents the mesh.
    Parameters:
    ev          - An N-by-1 array or a list of scalars. The Scalar values of the
                  elements. ev[i] should be the value of element i. 
    coords      - An N-by-2 or N-by-3 array. Row i contains the x,y,z coordinates
                  of node i.
    edof        - An E-by-L array. Element topology. (E is the number of elements
                  and L is the number of dofs per element)
    dofs_per_node - Integer. Dofs per node.
    el_type      - Integer. Element Type. See Gmsh manual for details. Usually 2 
                  for triangles or 3 for quadrangles.
    displacements - An N-by-2 or N-by-3 array. Row i contains the x,y,z 
                    displacements of node i.
    clim        - 2-tuple. Colorbar limits (min, max). Defines the value range of
                  the colorbar. Defaults to None, in which case min/max are set to
                  min/max of node_vals.
    axes        - Visvis Axes. The Axes where the model will be drawn. 
                  If unspecified the current Axes will be used, or a new Axes will
                  be created if none exist.
    axes_adjust  - Boolean. True if the view should be changed to show the whole 
                  model. Default True.
    draw_elements  - Boolean. True if mesh wire should be drawn. Default True.
    draw_undisplaced_mesh - Boolean. True if the wire of the undisplaced mesh 
                  should be drawn on top of the displaced mesh. Default False. 
                  Use only if displacements != None.
    magnfac     - Float. Magnification factor. Displacements are multiplied by
                  this value. Use this to make small displacements more visible.
    title       - String. Changes title of the figure. Default "Element Values".
    '''
    #Since vis.Mesh does not allow setting different colours for different faces, we need
    # a custom world object (WObject) for this function. 
    # http://code.google.com/p/visvis/wiki/example_customWobject
    # http://code.google.com/p/visvis/wiki/creatingWibjectsAndWobjects
    

    if draw_undisplaced_mesh:
        drawMesh(coords, edof, dofs_per_node, el_type, axes, axes_adjust, color=(0.5, 0.5, 0.5))
    
    if displacements is not None:
        if displacements.shape[1] != coords.shape[1]:
            displacements = np.reshape(displacements, (-1, coords.shape[1]))
            coords = np.asarray(coords + magnfac * displacements)
    
    axes, verts, faces, verticesPerFace, is3D = _preMeshDrawPrep(axes, coords, edof, dofs_per_node, el_type)
    
    
    #This is done because 3D elements are made up of several faces.
    #TODO: Discard inner faces that are not visible.
    fPerElms = { 1:0,   2:1,   3:1,   4:4,   5:6} #TODO: Extend with more element types
    facesPerElement = fPerElms[el_type]
    #Repeat the element values so that we get the value of each face:
    faceVals = np.repeat(ev, facesPerElement, axis=0) 

    c = _elementsWobject(axes, faceVals, verts, faces, verticesPerFace, draw_elements, clim) #Creates the world object that gets drawn on screen.
    
    _makeColorBar("Element values", axes) #Finds or creates colorbar and sets the label.
    
    # Adjust axes
    if axes_adjust:
        _adjustaxes(axes, is3D)
    
    vv.title(title, axes)
    return c

drawElementValues = draw_element_values
    
def draw_displacements(displacements, coords, edof, dofs_per_node, el_type, node_vals=None, clim=None, axes=None, 
                      axes_adjust=True, draw_undisplaced_mesh=True, magnfac=1.0,  title=None):
    '''
    Draws mesh with displacements in 2D or 3D. Scalar nodal values can also be 
    drawn on the mesh. Returns the displaced Mesh object.
    Parameters:
    displacements-An N-by-1 array (or matrix). Row i contains the displacement of
                  dof i. 
                  N-by-2 or N-by-3 arrays are also accepted, in which case row i 
                  contains the x,y,z displacements of node i.
    coords      - An N-by-2 or N-by-3 array. Row i contains the x,y,z coordinates
                  of node i.
    edof        - An E-by-L array. Element topology. (E is the number of elements
                  and L is the number of dofs per element)
    dofs_per_node - Integer. Dofs per node.
    el_type      - Integer. Element Type. See Gmsh manual for details. Usually 2 
                  for triangles or 3 for quadrangles.
    node_vals    - An N-by-1 array or a list of scalars. The Scalar values at the
                  nodes. node_vals[i] should be the value of node i.
    clim        - 2-tuple. Colorbar limits (min, max). Defines the value range of
                  the colorbar. Defaults to None, in which case min/max are set 
                  to min/max of node_vals.
    axes        - Visvis Axes. The Axes where the model will be drawn. 
                  If unspecified the current Axes will be used, or a new Axes will
                  be created if none exist.
    axes_adjust  - Boolean. True if the view should be changed to show the whole 
                  model. Default True.
    draw_elements  - Boolean. True if mesh wire should be drawn. Default True.
    magnfac     - Float. Magnification factor. Displacements are multiplied by
                  this value. Use this to make small displacements more visible.
    title       - String. Changes title of the figure. Default None (in which case
                  title depends on other parameters).
    '''
    
    if displacements.shape[1] != coords.shape[1]:
        displacements = np.reshape(displacements, (-1, coords.shape[1]))
    displaced = np.asarray(coords + magnfac * displacements)
       
    if draw_undisplaced_mesh:
        drawMesh(coords, edof, dofs_per_node, el_type, axes, axes_adjust, title=title, color=(0.5, 0.5, 0.5), filled=False)
    
    if node_vals != None:
        m = drawNodalValues(node_vals, displaced, edof, dofs_per_node, el_type, clim=clim, axes=axes, axes_adjust=axes_adjust, draw_elements=True, title=title)
    else:
        m = drawMesh(displaced, edof, dofs_per_node, el_type, axes, axes_adjust, title=title) 
    
    if title != None:
        vv.title(title, axes)
    return m

drawDisplacements = draw_displacements
    
def draw_geometry(geoData, axes=None, axes_adjust=True, draw_points=True, label_points=True, label_curves=True, title=None, font_size=11, N=20):
    '''
    Draws the geometry (points and curves) in geoData
    Parameters:
    geoData    - GeoData object. Geodata contains geometric information of the 
                 model.
    axes       - Visvis Axes. The Axes where the model will be drawn. 
                 If unspecified the current Axes will be used, or a new Axes will
                 be created if none exist.
    axes_adjust - Boolean. If True the view will be changed to show the whole 
                 model. Default True.
    draw_points - Boolean. If True points will be drawn.
    label_points- Boolean. If True Points will be labeled. The format is:
                 ID[marker]. If a point has marker==0 only the ID is written. 
    label_curves- Boolean. If True Curves will be labeled. The format is: 
                 ID(elementsOnCurve)[marker].
    font_size   - Integer. Size of the text in the text labels. Default 11. 
    N          - Integer. The number of discrete points per curve segment. 
                 Default 20. Increase for smoother curves. Decrease for better 
                 performance.
    '''
    
    if axes is None:
        axes = vv.gca()
    axes.bgcolor = (0.7, 0.7, 0.7)
    
    if draw_points:
        P = np.array(geoData.getPointCoords()) #M-by-3 list of M points.
        plotArgs = {'mc':'r', 'mw':5, 'lw':0, 'ms':'o', 'axes_adjust':False, 'axes':axes}
        if geoData.is3D: 
            vv.plot(P[:,0], P[:,1], P[:,2], **plotArgs)
        else:
            vv.plot(P[:,0], P[:,1], **plotArgs)           
        
        if label_points: #Write text label at the points:
            for (ID, (xyz, el_size, marker)) in geoData.points.items(): #[[x, y, z], el_size, marker]
                text = "  " + str(ID) + ("[%s]"%marker if marker is not 0 else '')
                addText(text, xyz, font_size=font_size-1, color=(0.5,0,0.5), axes=axes)  
    
    for(ID, (curveName, pointIDs, marker, elementsOnCurve, _, _)) in geoData.curves.items():
        points = geoData.getPointCoords(pointIDs)
        if curveName == "Spline":
            P = _catmullspline(points, N)
        if curveName == "BSpline":
            P = _bspline(points, N)
        if curveName == "Circle":
            P = _circleArc(*points, pointsOnCurve=N)
        if curveName == "Ellipse":
            P = _ellipseArc(*points, pointsOnCurve=N)
        plotArgs = {'lc':'k', 'ms':None, 'axes_adjust':False, 'axes':axes} #Args for plot style. Black lines with no symbols at points.
        if geoData.is3D:
            vv.plot(P[:,0], P[:,1], P[:,2], **plotArgs)
        else:
            vv.plot(P[:,0], P[:,1], **plotArgs)
        
        if label_curves:
            midP = P[int(P.shape[0]*7.0/12), :].tolist() # Sort of midpoint along the curve. Where the text goes.
            #Create the text for the curve. Includes ID, elementsOnCurve, and marker:
            text = " "+str(ID)
            text += "(%s)"%(elementsOnCurve) if elementsOnCurve is not None else ''
            text += "[%s]"%(marker) if marker is not 0 else '' #Something like "4(5)[8]"
            addText(text, midP, font_size=font_size, axes=axes)
        
    if title != None:
        vv.title(title, axes)
    
    if axes_adjust:
        _adjustaxes(axes, geoData.is3D)
    axes.daspectAuto = False
    axes.daspect = (1,1,1)
    
drawGeometry = draw_geometry

def _preMeshDrawPrep(axes, coords, edof, dofs_per_node, el_type):
    '''Duplicate code. Extracts verts, faces and verticesPerFace from input.'''
    if axes is None:
        axes = vv.gca() #Gets current Axis or creates a new one if none exists.
    
    if np.shape(coords)[1] == 2:
        is3D = False
        verts = np.hstack((coords, np.zeros([np.shape(coords)[0],1]))) #pad with zeros to make 3D
    elif np.shape(coords)[1] == 3:
        is3D = True
        verts = coords
    else:
        raise ValueError('coords must be N-by-2 or N-by-3 array')
    
    if el_type in [2, 4]: #elements with triangular faces    
        verticesPerFace = 3
    elif el_type in [3,5,16]: #elements with rectangular faces
        verticesPerFace = 4
    else:   #[NOTE] This covers all element types available in CALFEM plus tetrahedrons. If more element types are added it is necessary to include them here and below.
        raise ValueError('element type not implemented')
    
    faces = (edof[:,0::dofs_per_node]-1)/dofs_per_node  
	#'faces' here are actually lists of nodes in elements, not in faces necessarily if the elements are in 3D. This case is handled below.   
    
    if el_type in [4,5]: #if hexahedrons or tetrahedrons:
        if  el_type == 5:
            G = np.array([[0,3,2,1],
                       [0,1,5,4],
                       [4,5,6,7],
                       [2,6,5,1],
                       [2,3,7,6],
                       [0,4,7,3]]) #G is an array that is used to decomposes hexahedrons into its component faces.
					   #The numbers are from the node orders (see p94 in the Gmsh manual) and each row makes one face.
        elif el_type == 4:
            G = np.array([[0,1,2],
                       [0,3,2],
                       [1,3,2],
                       [0,3,1]]) #This G decomposes tetrahedrons into faces
        faces = np.vstack([ faces[i, G] for i in range(faces.shape[0]) ])
    elif el_type == 16: #if 8-node-quads:
        faces = faces[:, 0:4] #The first 4 nodes are the corners of the high order quad.
        
    axes.bgcolor = (0.7, 0.7, 0.7) #background colour.
    return axes, verts, faces, verticesPerFace, is3D


def _adjustaxes(axes, is3D):
    if axes.daspectAuto is None:
            axes.daspectAuto = False
    axes.cameraType = '3d' if is3D else '2d'
    axes.SetLimits(margin=0.1)


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
    controlPoints = np.asarray(controlPoints) #Convert to array if input is a list.
    if (controlPoints[0,:] == controlPoints[-1,:]).all():
        #If the curve is closed we extend each opposite endpoint to the other side  
        CPs = np.asmatrix(np.vstack((controlPoints[-2,:],
                                controlPoints,
                                controlPoints[1,:])))
    else: #Else make mirrored endpoints:
        CPs = np.asmatrix(np.vstack((2*controlPoints[0,:] - controlPoints[1,:],
                                controlPoints,
                                2*controlPoints[-1,:] - controlPoints[-2,:])))
    M = 0.5 * np.matrix([[ 0,  2,  0,  0],[-1,  0,  1,  0],[ 2, -5,  4, -1],[-1,  3, -3,  1]])
    t = np.linspace(0, 1, pointsOnEachSegment)
    T = np.matrix([[1, s, pow(s,2), pow(s,3)] for s in t])
    return np.asarray( np.vstack( T * M * CPs[j-1:j+3,:] for j in range( 1, len(CPs)-2 ) ) )


def _bspline(controlPoints, pointsOnCurve=20):
    '''
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
    '''
    controlPoints = np.asarray(controlPoints) #Convert to array if input is a list.
    if (controlPoints[0,:] == controlPoints[-1,:]).all():
        #If the curve is closed we extend each opposite endpoint to the other side  
        CPs = np.asmatrix(np.vstack((controlPoints[-2,:],
                                controlPoints,
                                controlPoints[1,:])))
    else:#Else make mirrored endpoints:
        CPs = np.asmatrix(np.vstack((2*controlPoints[0,:] - controlPoints[1,:],
                                controlPoints,
                                2*controlPoints[-1,:] - controlPoints[-2,:])))
    M = (1.0/6) * np.matrix([[-1,  3, -3, 1],
                          [ 3, -6,  3, 0],
                          [-3,  0,  3, 0],
                          [ 1,  4,  1, 0]])
    t = np.linspace(0, 1, pointsOnCurve)
    T = np.matrix([[pow(s,3), pow(s,2), s, 1] for s in t])
    return np.asarray( np.vstack( T * M * CPs[i-1 : i+3, :] for i in range( 1, len(CPs)-2 ) ) )

def _circleArc(start, center, end, pointsOnCurve=20):
    return _ellipseArc(start, center, start, end, pointsOnCurve)

def _ellipseArc(start, center, majAxP, end, pointsOnCurve=20):
    '''Input are 3D 1-by-3 numpy arrays or vectors'''
    #First part is to find a similarity transform in 3D that transform the ellipse to
    #the XY-plane with the center at the origin and the major axis of the ellipse along the X-axis.
    
    #convert to arrays in case inputs are lists:
    start, center, majAxP, end, = np.asarray(start), np.asarray(center), np.asarray(majAxP), np.asarray(end)
    
    zPrim = np.cross(start-center, end-center)
    zPrim = zPrim / np.linalg.norm(zPrim)
    xPrim = (majAxP-center) / np.linalg.norm(majAxP-center)
    yPrim = np.cross(zPrim, xPrim)
    
    R = np.vstack((xPrim, yPrim, zPrim)).T #Rotation matrix from ordinary coords to system where ellipse is in the XY-plane. (Actually hstack)
    T = np.hstack((R, np.asmatrix(center).T))   #Building Transformation matrix. -center is translation vector from ellipse center to origin.
    T = np.mat( np.vstack((T, [0,0,0,1])) ) #Transformation matrix for homogenous coordinates.
    
    startHC = np.vstack((np.matrix(start).T, [1])) #
    endHC = np.vstack((np.matrix(end).T, [1]))     # start and end points as column vectors in homogenous coordinates
    
    s = np.linalg.inv(T) * startHC #
    e = np.linalg.inv(T) * endHC   # start and end points in the new coordinate system
    
    xs, ys = s[0,0], s[1,0] #
    xe, ye = e[0,0], e[1,0] # Just extract x & y from the new start and endpoints  
    
    a = np.sqrt( (pow(ye*xs,2) - pow(xe*ys,2)) / (pow(ye,2) - pow(ys,2)) )
    b = np.sqrt( (pow(ye*xs,2) - pow(xe*ys,2)) / ( (pow(ye,2) - pow(ys,2)) * ((pow(xe,2) - pow(xs,2)) / (pow(ys,2) - pow(ye,2)) ) ) )
    
    ts = atan2(ys/b, xs/a) #atan2 is a function that goes from -pi to pi. It gives the signed angle from the X-axis to point (y,x)
    te = atan2(ye/b, xe/a) #We can't use the (transformed) start- and endpoints directly, but we divide x and y by the
    # ellipse minor&major axes to get the parameter t that corresponds to the point on the ellipse. 
    # See ellipse formula: x = a * cos (t), y = b * sin(t). 
    # So ts and te are the parameter values of the start- and endpoints (in the transformed coordinate system).
    
    if ts > te:
        ts, te = te, ts #swap if the start point comes before the endpoint in the parametric parameter that goes around the ellipse.
    if te - ts < np.pi:
        times = np.linspace(ts, te, pointsOnCurve) #parameter of ellipse curve. NOT angle to point on curve (like it could be for a circle).
    else: #the shortest parameter distance between start- and end-point stradles the discontinuity that jumps from pi to -pi.  
        ps1 = round(pointsOnCurve * (pi-te)/(2*pi-te+ts)) #number of points on the first length.
        ps2 = round(pointsOnCurve * (ts+pi)/(2*pi-te+ts)) #number of points on the first length.
        times = np.concatenate((np.linspace(te, pi, ps1), np.linspace(-pi, ts, ps2)))
    
    ellArc = np.array([[a*cos(t), b*sin(t)] for t in times]).T #points on arc (in 2D)
    ellArc = np.vstack((ellArc, np.repeat(np.matrix([[0],[1]]), ellArc.shape[1], 1))) #Make 3D homogenous coords by adding rows of 0s and 1s.
    ellArc = T * ellArc #Transform back to the original coordinate system   
    return np.asarray(ellArc.T[:,0:3]) #return points as an N-by-3 array. 

class _elementsWobject(vv.Wobject, Colormapable):
    '''
    Custom wobject for drawing element values.  
    Based on example http://code.google.com/p/visvis/wiki/example_customWobject
    '''
    #TODO: Find a way to make rendering faster. (Dump internal faces, cull backfaces, 
    #TODO: Don't draw lines separately, pass array to GL instead of looping, etc?)
    def __init__(self, parent, fVals, verts, faces, verticesPerFace, draw_elements, clim=None):
        vv.Wobject.__init__(self, parent)
        self._fVals = fVals # N-by-1 array?            values of N faces (not elements if 3D)
        self._verts = verts # M-by-3 array             coordinates of M vertices  
        self._faces = faces # N-by-3 or N-by-4 array   verts that make up N faces
        self.verticesPerFace = verticesPerFace #3 or 4. Either we have triangle faces or quads.
        self.draw_elements = draw_elements
        self._valMin = np.amin(fVals)
        self._valMax = np.amax(fVals)
        Colormapable.__init__(self)
        self._texture = None
        self.clim = minmax(self._fVals) if clim is None else clim
        self.colormap = vv.colormaps['jet']
    
    def _drawFaces(self, how):
        for (value, faceVerts) in zip(self._fVals, self._faces):
            valueIndex = int( 255 * (value-self.clim.min)/(self.clim.max-self.clim.min) ) #Turn the value into an index between 0 and 255
            if valueIndex > 255:
                valueIndex = 255
            elif valueIndex < 0:
                valueIndex = 0
            color = self.mapData[valueIndex, :] #get colour from value
            coords = self._verts[faceVerts.astype(int), :] #get coordinates of the vertices of the face.
            gl.glColor(*color)
            gl.glBegin(how)
            for (x,y,z) in coords:
                gl.glVertex(x,y,z)
            gl.glEnd()
            
    def _drawLines(self, how=gl.GL_LINE_LOOP, color=(0,0,0)):
        for faceVerts in self._faces:
            coords = self._verts[faceVerts.astype(int), :] #get coordinates of the vertices of the face.
            gl.glColor(*color)
            gl.glBegin(how)
            for (x,y,z) in coords:
                gl.glVertex(x,y,z)
            gl.glEnd()
            
    def _GetLimits(self):
        """ Tell the axes how big this object is.
        """ 
        # Get limits
        x1, x2 = minmax(self._verts[:,0])
        y1, y2 = minmax(self._verts[:,1])
        z1, z2 = minmax(self._verts[:,2])
        return vv.Wobject._GetLimits(self, x1, x2, y1, y2, z1, z2)
    
    def OnDraw(self):
        """ To draw the object.
        """ 
        if self.draw_elements:
            gl.glDisable(gl.GL_LINE_SMOOTH)
            #gl.glEnable(gl.GL_BLEND)
            gl.glLineWidth(1)
            self._drawLines(gl.GL_LINE_LOOP, (0,0,0))
        if self.verticesPerFace == 3:
            self._drawFaces(gl.GL_TRIANGLES)
        elif self.verticesPerFace == 4:
            self._drawFaces(gl.GL_QUADS)
    
    def OnDrawShape(self, clr):
        """ To draw the shape of the object.
        Only necessary if you want to be able to "pick" this object
        """
        self._drawFaces(gl.GL_TRIANGLES, clr)
    
    def OnDrawScreen(self):
        """ If the object also needs to draw in screen coordinates.
        Text needs this for instance.
        """
        pass
    
    def OnDestroyGl(self):
        """ To clean up any OpenGl resources such as textures or shaders.
        """
        pass
    
    def OnDestroy(self):
        """ To clean up any other resources.
        """
        pass
    
    #Overload _SetColormap get-setter so that we can change self.mapData, which we use 
    #to map values to colours without calling self._colormap.GetData() every draw call.
    def _SetColormap(self, value):
        self._colormap.SetMap(value)
        self.mapData = self._colormap.GetData()
        
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
    
    line_type = plotpar[0]
    line_color = plotpar[1]
    node_mark = plotpar[2]
    
    # Translate CALFEM plotpar to visvis
    
    vv_line_type = '-'
    vv_line_color = 'b'
    vv_node_mark = 'o'
    
    if line_type == 1:
        vv_line_type = '-'
    elif line_type == 2:
        vv_line_type = '--'
    elif line_type == 3:
        vv_line_type = ':'
        
    if line_color == 1:
        vv_line_color = 'k'
    elif line_color == 2:
        vv_line_color = 'b'
    elif line_color == 3:
        vv_line_color = 'm'
    elif line_color == 4:
        vv_line_color = 'r'
        
    if node_mark == 1:
        vv_node_mark = 'o'
    elif node_mark == 2:
        vv_node_mark = 'x'
    elif node_mark == 0:
        vv_node_mark = ''
        
    vv_marker_color = vv_line_color

    vv.axis('equal')
    
    draw_element_numbers = False
    
    if len(elnum) == ex.shape[0]:
        draw_element_numbers = True

    i = 0

    for elx, ely in zip(ex, ey):
        x = elx.tolist()
        x.append(elx[0])
        y = ely.tolist()
        y.append(ely[0])
        
        xm = sum(x)/len(x)
        ym = sum(y)/len(y)
        
        vv.plot(x, y, ls=vv_line_type, lc = vv_line_color, ms = vv_node_mark, mc = vv_marker_color)
        
        if draw_element_numbers:
            text(str(elnum[i]), [xm, ym])
            i += 1
            
def eldraw2_mpl(ex, ey, plotpar=[1, 2, 1], elnum=[]):
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
    
    line_type = plotpar[0]
    line_color = plotpar[1]
    node_mark = plotpar[2]
    
    # Translate CALFEM plotpar to visvis
    
    vv_line_type = '-'
    vv_line_color = 'b'
    vv_node_mark = 'o'
    
    if line_type == 1:
        vv_line_type = '-'
    elif line_type == 2:
        vv_line_type = '--'
    elif line_type == 3:
        vv_line_type = ':'
        
    if line_color == 1:
        vv_line_color = 'k'
    elif line_color == 2:
        vv_line_color = 'b'
    elif line_color == 3:
        vv_line_color = 'm'
    elif line_color == 4:
        vv_line_color = 'r'
        
    if node_mark == 1:
        vv_node_mark = 'o'
    elif node_mark == 2:
        vv_node_mark = 'x'
    elif node_mark == 0:
        vv_node_mark = ''
        
    vv_marker_color = vv_line_color

    plt.axis('equal')
    
    draw_element_numbers = False
    
    if len(elnum) == ex.shape[0]:
        draw_element_numbers = True

    i = 0

    for elx, ely in zip(ex, ey):
        x = elx.tolist()
        x.append(elx[0])
        y = ely.tolist()
        y.append(ely[0])
        
        xm = sum(x)/len(x)
        ym = sum(y)/len(y)
        
        plt.plot(x, y, vv_line_color + vv_node_mark + vv_line_type)
        #vv.plot(x, y, ls=vv_line_type, lc = vv_line_color, ms = vv_node_mark, mc = vv_marker_color)
        
        #if draw_element_numbers:
        #    text(str(elnum[i]), [xm, ym])
        #    i += 1          
            
def eliso2_mpl(ex, ey, ed):

    plt.axis('equal')
    
    print(np.shape(ex))    
    print(np.shape(ey))
    print(np.shape(ed))
    
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
        
def eldraw2_old(ex, ey):
    """
    Draw elements in 2d.
    
    Parameters:
    
        ex, ey          Element coordinates
        plotpar         (not implemented yet)
    
    """
    #if not haveWx:
    #    print("wxPython not installed.")
    #    return
    
    #class ElDispApp(wx.App):
    #    def OnInit(self):
    #        wx.InitAllImageHandlers()
    #        mainWindow = ElementView(None, -1, "")
    #        mainWindow.ex = ex
    #        mainWindow.ey = ey
    #        mainWindow.showNodalValues = False
    #        self.SetTopWindow(mainWindow)
    #        mainWindow.Show()
    #        return 1
    #
    #app = ElDispApp(0)
    #app.MainLoop()
    mainWindow = ElementView(None, -1, "")
    mainWindow.ex = ex
    mainWindow.ey = ey
    mainWindow.showNodalValues = False
    mainWindow.Show()   
    globalWindows.append(mainWindow)
    
def eliso2_old(ex, ey, ed, showMesh=False):
    """
    Draw nodal values in 2d.
    
    Parameters:
    
        ex, ey          Element coordinates
        ed              Element nodal values
        plotpar         (not implemented yet)
    
    """
    #if not haveWx:
    #    print("wxPython not installed.")
    #    return
    
    #class ElDispApp(wx.App):
    #    def OnInit(self):
    #        wx.InitAllImageHandlers()
    #        mainWindow = ElementView(None, -1, "")
    #        mainWindow.ex = ex
    #        mainWindow.ey = ey
    #        mainWindow.ed = ed
    #        mainWindow.showMesh = showMesh
    #        mainWindow.showNodalValues = True
    #        self.SetTopWindow(mainWindow)
    #        mainWindow.Show()
    #        return 1   
    #
    #app = ElDispApp(0)
    #app.MainLoop()
    mainWindow = ElementView(None, -1, "")
    mainWindow.ex = ex
    mainWindow.ey = ey
    mainWindow.ed = ed
    mainWindow.showMesh = showMesh
    mainWindow.showNodalValues = True
    mainWindow.Show()
    globalWindows.append(mainWindow)
    
def elval2(ex, ey, ev, showMesh=False):
    """
    Draw elements values in 2d.
    
    Parameters:
    
        ex, ey          Element coordinates
        ev              Element values (scalar)
        plotpar         (not implemented yet)
    
    """
    #if not haveWx:
    #    print("wxPython not installed.")
    #    return
    
    mainWindow = ElementView(None, -1, "")
    mainWindow.ex = ex
    mainWindow.ey = ey
    mainWindow.ev = ev
    mainWindow.showMesh = showMesh
    mainWindow.showElementValues = True
    mainWindow.showNodalValues = False
    mainWindow.Show()
    globalWindows.append(mainWindow)
    
def eldisp2(ex, ey, ed, magnfac=0.1, showMesh=True):
    #if not haveWx:
    #    print("wxPython not installed.")
    #    return
        
    mainWindow = ElementView(None, -1, "")
    mainWindow.dofs_per_node = 2
    mainWindow.ex = ex
    mainWindow.ey = ey
    mainWindow.ed = ed
    mainWindow.showMesh = showMesh
    mainWindow.showNodalValues = False
    mainWindow.showDisplacements = True
    mainWindow.magnfac = magnfac
    mainWindow.Show()    
    globalWindows.append(mainWindow)
    
def waitDisplayNative():
    if haveQt:
        globalQtApp.exec_()        
    else:
        globalWxApp.MainLoop()
        

def showAndWaitNative():
    if haveQt:
        globalQtApp.exec_()
    else:
        globalWxApp.MainLoop()

if __name__ == "__main__":
    print("Testing")
