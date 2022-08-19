# -*- coding: utf-8 -*-

"""
CALFEM Vedo

Module for 3D visualization in CALFEM using Vedo (https://vedo.embl.es/)

@author: Andreas Åmand
"""

import numpy as np
import vedo as v
import pyvtk
import vtk
import sys
import time
from scipy.io import loadmat
import calfem.core as cfc
import vis_vedo_utils as vdu
import vtk

# Examples using this module:
    # exv1: Spring
    # exv2: 3D Truss (beams & bars)
    # exv3: 3D Flow
    # exv4: 3D Solid
    # exv5: 2D Plate

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
# Initialization


def init_app():
    """
    Routine for calling the VedoPlotWindow class via show_and_wait.
    """
    global vedo_app
    vedo_app = VedoPlotWindow.instance()
    if vedo_app is None: vedo_app = VedoMainWindow()
    return vedo_app

class VedoPlotWindow:
    """
    Class for handling the VedoMainWindow class and making sure there is only one instance of it.
    """
    __instance = None
    @staticmethod
    def instance():
        """ Static access method. """
        if VedoPlotWindow.__instance == None: VedoPlotWindow()
        return VedoPlotWindow.__instance

    def __init__(self):
        """ Virtually private constructor. """
        if VedoPlotWindow.__instance is None:
            VedoPlotWindow.__instance = self
            self.plot_window = VedoMainWindow()

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
# Main window


class VedoMainWindow():
    """
    Class for handling the rendering and interaction of figures.
    """
    
    def __init__(self):
        # Variables to keep track of plotters
        self.fig = 0
        self.bg = []
        self.mode_2D = []
        self.mode_hover = []
        self.mode_anim = []
        self.rendered = 0
        self.plt = {}

        # Variables for handeling main objects to render
        self.geometries = [[]]
        self.meshes = [[]]
        self.nodes = [[]]

        # Variables for handeling optional objects to render
        self.msg = [[]]
        self.proj = [[]]
        self.rulers = [[]]
        self.vectors = [[]]
        self.dia_lines = [[]]
        self.dia_points = [[]]
        self.dia_axes = [[]]

        # Variables for handeling animations
        self.keyframes = [[]]
        self.keyframe_dict = [[]]
        self.loop = None
        self.dt = 500 #milliseconds

        # Mouse click callback
        self.silcont = [None]
        self.click_msg = v.Text2D("", pos="bottom-center", bg='auto', alpha=0.1, font='Calco',c='black')

        # Global settings
        v.settings.immediateRendering = False
        
        #v.settings.allowInteraction = True
        v.settings.allowInteraction = False
        #v.settings.renderLinesAsTubes = True
        v.settings.allowInteraction = True
        v.settings.useFXAA = True
        v.settings.useSSAO         = True
        v.settings.visibleGridEdges = True

        # Fonts that could be interesting to use
        #v.settings.defaultFont = 'Normografo'
        #v.settings.defaultFont = 'LogoType'
        #v.settings.defaultFont = 'Courier'
        #v.settings.defaultFont = 'Comae'
        #v.settings.defaultFont = 'Calco'

    def click(self,evt):
        """
        Routine for handling mouse clicks and highlighting elements/nodes in draw_mesh.
        NOTE: not to be run by a user
        """
        if evt.actor:
            # Silouette
            sil = evt.actor.silhouette().lineWidth(5).c('red5')
            self.plt[evt.title].remove(self.silcont.pop()).add(sil)
            self.silcont.append(sil)
            self.click_msg.text(evt.actor.name)
            
        else:
            self.click_msg.text('')
            if self.silcont != [None]:
                sil = None
                self.plt[evt.title].remove(self.silcont.pop()).add(sil)
                self.silcont.append(sil)
            return

    def render(self):
        """
        Routine for adding rendering objects in a figure.
        NOTE: not to be run by a user
        """

        while len(self.meshes) < self.fig + 1:
            self.geometries.append([])
            self.meshes.append([])
            self.nodes.append([])

            self.msg.append([])
            self.proj.append([])
            self.rulers.append([])
            self.vectors.append([])

            self.dia_lines.append([])
            self.dia_points.append([])
            self.dia_axes.append([])

            self.keyframes.append([])
            self.keyframe_dict.append([])

        while len(self.mode_2D) < self.fig + 1:
            self.bg.append('w')
            self.mode_2D.append(False)
            self.mode_hover.append(False)
            self.mode_anim.append(False)

        for plot in range(self.rendered, self.fig+1):
            if self.mode_2D[plot] is True:
                pp = True
                mode = 5
            else:
                pp = False
                mode = 0
            
            if self.mode_anim[plot] == True:
                v.settings.immediateRendering = True
                opts = dict(axes=4, interactive=0, title=f'Figure {plot+1} - CALFEM vedo visualization tool')
                keyframes = self.keyframes[self.fig]
                keyframe_dict = self.keyframe_dict[self.fig]
                nmesh = len(keyframes)
                plt = v.Plotter(**opts)

                plt += v.Text2D('Press ESC to exit', pos='bottom-middle')


                for j in range(len(self.msg[plot])):
                    if self.bg[self.fig] == 'black':
                        self.msg[plot][j].c('w')
                    plt.add(self.msg[plot][j])

                for j in range(len(self.proj[plot])):
                    if self.bg[self.fig] == 'black':
                        self.proj[plot][j].c('w')
                    plt.add(self.proj[plot][j])

                for j in range(len(self.rulers[plot])):
                    if self.bg[self.fig] == 'black':
                        self.rulers[plot][j].c('w')
                    plt.add(self.rulers[plot][j])

                for j in range(len(self.vectors[plot])):
                    plt.add(self.vectors[plot][j])

                for j in range(len(self.dia_lines[plot])):
                    plt.add(self.dia_lines[plot][j])

                for j in range(len(self.dia_points[plot])):
                    if self.bg[self.fig] == 'black':
                        self.dia_points[plot][j].c('w')
                    plt.add(self.dia_points[plot][j])

                for j in range(len(self.dia_axes[plot])):
                    if self.bg[self.fig] == 'black':
                        lst = self.dia_axes[plot][j].unpack()
                        for k in lst:
                            k.c('w')
                    plt.add(self.dia_axes[plot][j])


                plt += v.Text2D('Press ESC to exit', pos='bottom-middle')
                plt += keyframes[0][keyframe_dict[0][1]]
                plt.show(resetcam=True)
                
                while True:
                    
                    for key, val in keyframe_dict[0].items():

                        time.sleep(self.dt/1000)
                        for i in range(nmesh):
                            plt.pop()
                        

                        for i in range(nmesh):
                            plt += keyframes[i][val]
                        plt.show(resetcam=False)

                        
                    if plt.escaped: break  # if ESC is hit during the loop

                plt.close()
                v.settings.immediateRendering = False

            
            else:
                if self.bg[self.fig] != 'white':
                    opts = dict(bg=self.bg[self.fig], mode=mode, axes=4, interactive=False, new=True, title=f'Figure {plot+1} - CALFEM vedo visualization tool')
                    if self.bg[self.fig] == 'black':
                        self.click_msg.c('w')
                        for j in range(len(self.geometries[plot])):
                            for k in range(len(self.geometries[plot][j])):
                                self.geometries[plot][j][k].c('w')
                else:
                    opts = dict(mode=mode, axes=4, interactive=False, new=True, title=f'Figure {plot+1} - CALFEM vedo visualization tool')
                plt = v.show(self.geometries[plot], self.meshes[plot], self.nodes[plot], self.click_msg, **opts)
                plt.parallelProjection(value=pp)

            if self.mode_hover[plot] == True or self.mode_2D[plot] == True:
                plt.addCallback('MouseMove', self.click)
            else:
                plt.addCallback('mouse click', self.click)
            
            self.plt[f'Figure {plot+1} - CALFEM vedo visualization tool'] = plt

            for j in range(len(self.msg[plot])):
                if self.bg[self.fig] == 'black':
                    self.msg[plot][j].c('w')
                plt.add(self.msg[plot][j])

            for j in range(len(self.proj[plot])):
                if self.bg[self.fig] == 'black':
                    self.proj[plot][j].c('w')
                plt.add(self.proj[plot][j])

            for j in range(len(self.rulers[plot])):
                if self.bg[self.fig] == 'black':
                    self.rulers[plot][j].c('w')
                plt.add(self.rulers[plot][j])

            for j in range(len(self.vectors[plot])):
                plt.add(self.vectors[plot][j])

            for j in range(len(self.dia_lines[plot])):
                plt.add(self.dia_lines[plot][j])

            for j in range(len(self.dia_points[plot])):
                if self.bg[self.fig] == 'black':
                    self.dia_points[plot][j].c('w')
                plt.add(self.dia_points[plot][j])

            for j in range(len(self.dia_axes[plot])):
                if self.bg[self.fig] == 'black':
                    lst = self.dia_axes[plot][j].unpack()
                    for k in lst:
                        k.c('w')
                plt.add(self.dia_axes[plot][j])

            self.rendered += 1

        v.interactive()





### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
# Main functions for plotting geometries, meshes, animations & results from CALFEM calculations


def draw_geometry(points=None,lines=None,surfaces=None,scale=0.05,points_alpha=1,lines_alpha=1,surfaces_alpha=1):
    """
    Routine for geometry of spring, bar, flow, solid, beam or plate elements.

    :param dict points: Points from CALFEM goemetry module
    :param dict lines: Lines from CALFEM goemetry module
    :param dict surfaces: Surfaces from CALFEM goemetry module
    :param float scale: Text scale
    :param float points_alpha: Transparency for points
    :param float lines_alpha: Transparency for lines
    :param float surfaces_alpha: Transparency for surfaces
    """
    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window

    if surfaces == None and lines == None and points == None:
        print("draw_geometry: Please input either (points), (points, lines) or (points, lines, surfaces) from geometry module")
        sys.exit()
    else:
        if surfaces is not None:

            pts = []
            p_text = []
            for key,val in points.items():
                pts.append(v.Point(pos=val[0], r=12, c='black', alpha=points_alpha))
                p_text.append(v.Text3D(key, pos=val[0], s=scale, font='Normografo', hspacing=1.15, vspacing=2.15, depth=0, italic=False, justify='bottom-left', c='black', alpha=points_alpha, literal=False))

            plot_window.geometries[plot_window.fig].append(pts)
            plot_window.geometries[plot_window.fig].append(p_text)

            l = []
            l_text = []
            for key,val in lines.items():
                p1 = points[val[1][0]][0]
                p2 = points[val[1][1]][0]
                l.append(v.Lines([p1], [p2], c='k4', alpha=lines_alpha, lw=1, res=2))
                l_text.append(v.Text3D(key, pos=[(p1[0]+0.5*(p2[0]-p1[0])), (p1[1]+0.5*(p2[1]-p1[1])), (p1[2]+0.5*(p2[2]-p1[2]))], s=scale, font='Normografo', hspacing=1.15, vspacing=2.15, depth=0, italic=False, justify='bottom-left', c='black', alpha=lines_alpha, literal=False))

            plot_window.geometries[plot_window.fig].append(l)
            plot_window.geometries[plot_window.fig].append(l_text)

            npts = len(val[1])

            

            surf = []
            for key,val in surfaces.items():
                surface_lines = val[1]
                

                point_key_dict = {}
                npts = len(val[1])
                indx = 0
                for line in surface_lines:
                    point_key_dict[indx] = lines[line][1]
                    indx += 1

                point_list = []
                for i in range(npts):
                    l0 = point_key_dict[i]
                    if i == npts - 1:
                        l1 = point_key_dict[0]
                    else:
                        l1 = point_key_dict[i+1]

                    if l0[1] == l1[0] or l0[1] == l1[1]:
                        point_list.append(l0[0])
                    elif l0[0] == l1[0] or l0[0] == l1[1]:
                        point_list.append(l0[1])
                    else:
                        print('Error when rendering surface geometry')
                        sys.exit()

                coords = []
                for i in point_list:
                    coords.append(points[i][0])

                pln = v.Mesh([coords, [[0,1,2,3]]], alpha=surfaces_alpha, c='grey6')
                
                surf.append(pln)
                pln.name = f"Surface {key}"

            plot_window.geometries[plot_window.fig].append(surf)
            
        elif lines is not None and points is not None:
            pts = []
            p_text = []
            for key,val in points.items():
                pts.append(v.Point(pos=val[0], r=12, c='black', alpha=points_alpha))
                p_text.append(v.Text3D(key, pos=val[0], s=scale, font='Normografo', hspacing=1.15, vspacing=2.15, depth=0, italic=False, justify='bottom-left', c='black', alpha=points_alpha, literal=False))

            plot_window.geometries[plot_window.fig].append(pts)
            plot_window.geometries[plot_window.fig].append(p_text)

            l = []
            l_text = []
            for key,val in lines.items():
                p1 = points[val[1][0]][0]
                p2 = points[val[1][1]][0]
                l.append(v.Lines([p1], [p2], c='k4', alpha=lines_alpha, lw=1, res=2))
                l_text.append(v.Text3D(key, pos=[(p1[0]+0.5*(p2[0]-p1[0])), (p1[1]+0.5*(p2[1]-p1[1])), (p1[2]+0.5*(p2[2]-p1[2]))], s=scale, font='Normografo', hspacing=1.15, vspacing=2.15, depth=0, italic=False, justify='bottom-left', c='black', alpha=lines_alpha, literal=False))

            plot_window.geometries[plot_window.fig].append(l)
            plot_window.geometries[plot_window.fig].append(l_text)
        elif points is not None:
            pts = []
            text = []

            for key,val in points.items():
                pts.append(v.Point(pos=val[0], r=12, c='black', alpha=points_alpha))
                text.append(v.Text3D(key, pos=val[0], s=scale, font='Normografo', hspacing=1.15, vspacing=2.15, depth=0, italic=False, justify='bottom-left', c='black', alpha=points_alpha, literal=False))

            plot_window.geometries[plot_window.fig].append(pts)
            plot_window.geometries[plot_window.fig].append(text)
        elif lines is not None:
            print("draw_geometry: Please provide point coordinates along with lines")
            sys.exit()

    
'''
Element types:
    1: Spring
    2: Bar
    3: Flow
    4: Solid
    5: Beam
    6: Plate

Required input for draw_mesh(), draw_displaced_mesh() & animation()
'''


def draw_mesh(
    # Main input
    edof,
    coord,
    dof,
    element_type,

    # Other parameters
    scale = 0.02,
    alpha = 1,
    render_nodes = True,
    color = 'yellow',
    offset = [0, 0, 0],

    # BC- & F-marking
    bcPrescr = None,
    bc = None,
    bc_color = 'red',
    fPrescr = None,
    f = None,
    f_color = 'blue6',
    eq_els = None,
    eq = None,

    # Element-specific input
    spring = True,
    nseg = 2
    ):
    """
    Routine for undisplaced mesh for spring, bar, flow, solid, beam or plate elements.

    :param array edof: Element topology by degrees of freedom [nel x (n_dofs_per_node)|(n_dofs_per_node+1)*n_nodes ]
    :param array coord: Nodal coordinates [number of nodes x 3]
    :param array dof: Global degrees of freedom [number of nodes x degrees of freedom per node]
    :param int element_type: Element type [1-6]
    :param float scale: Element scale, nodes are scaled 50% larger than this value
    :param float alpha: Element and node transparency [0-1]
    :param bool render_nodes: If True, nodes are rendered
    :param str color: Element color
    :param list offset: Offset actors in 3D space [x, y, z]
    :param array bcPrescr: Degrees of freedom with boundary conditions [number of degrees of freedom with BCs x 1]
    :param array bc: Boundary conditions [number of degrees of freedom with BCs x 1]
    :param str bc_color: Color for nodes with boundary conditions applied
    :param array fPrescr: Degrees of freedom with forces [number of degrees of freedom with forces x 1]
    :param array f: Forces at degrees of freedom [number of degrees of freedom with forces x 1]
    :param str f_color: Color for nodes/elements with forces applied
    :param array eq_els: Element numbers where forces are applied [number of elements with forces x 1]
    :param array eq: Element force vector [number of elements with forces x 1 | number of elements with forces x 3]
    :param bool spring: If True, renders spring elements as coil springs
    :param int nseg: Number of points along beam elements for segmenting [number of segments + 1]

    :return array mesh: List of mesh actors
    """

    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window

    if np.size(coord, axis = 1) == 1:
        coord = np.append(coord, np.zeros((np.size(coord, axis = 0),1)), axis=1)
        coord = np.append(coord, np.zeros((np.size(coord, axis = 0),1)), axis=1)
    elif np.size(coord, axis = 1) == 2:
        coord = np.append(coord, np.zeros((np.size(coord, axis = 0),1)), axis=1)

    if 1 <= element_type <= 6:
        nel, ndof_per_el, nnode, ndim, ndof, ndof_per_n = vdu.check_input(edof,coord,dof,element_type,nseg=nseg)
    else:
        print("draw_mesh: Invalid element type, please declare 'element_type'. The element types are:\n    1 - Spring\n    2 - Bar\n    3 - Flow\n    4 - Solid\n    5 - Beam\n    6 - Plate")
        sys.exit()

    # OUTPUT FROM check_input
    # Number of elements:                       nel
    # Number of degrees of freedom per element: ndof_per_el
    # Number of nodes:                          nnode
    # Number of dimensions:                     ndim
    # Number of degrees of freedom:             ndof
    # Number of degrees of freedom per node:    ndof_per_n
    # Number of displacements:                  ndisp
    # Element/nodal values:                     val

    # Elements w/ a length (spring, bar & beam)
    if element_type == 1 or element_type == 2 or element_type == 5:           

        ncoord = np.size(coord, axis = 0)
        nel = np.size(edof, axis = 0)
        elements = []

        coord[:] += offset

        if element_type == 5:
            for i in range(nel):
                eq_dict = {}
                indx = 0
                if isinstance(eq_els, np.ndarray):
                    for j in eq_els:
                        eq_dict[j[0]] = eq[indx][0]
                        indx += 1

        for i in range(nel):
            coord1,coord2 = vdu.get_coord_from_edof(edof[i,:],dof,element_type)
            if element_type == 1 and spring == True:
                element = v.Spring([coord[coord1,0],coord[coord1,1],coord[coord1,2]],[coord[coord2,0],coord[coord2,1],coord[coord2,2]],r=1.5*scale,c=color).alpha(alpha)
                element.name = f"Spring element {i+1}"
                elements.append(element)
            elif element_type == 1 and spring == False:
                element = v.Cylinder([[coord[coord1,0],coord[coord1,1],coord[coord1,2]],[coord[coord2,0],coord[coord2,1],coord[coord2,2]]],r=scale,res=4,c=color).alpha(alpha)
                element.name = f"Spring element {i+1}"
                elements.append(element)
            elif element_type == 2:
                bar = v.Cylinder([[coord[coord1,0],coord[coord1,1],coord[coord1,2]],[coord[coord2,0],coord[coord2,1],coord[coord2,2]]],r=scale,res=4,c=color).alpha(alpha)
                bar.name = f"Bar element {i+1}"
                elements.append(bar)
            # Segmented beam
            elif element_type == 5 and nseg > 2:
                steps = np.float32(1/(nseg-1))

                dx = (coord[coord2,0]-coord[coord1,0])*steps
                dy = (coord[coord2,1]-coord[coord1,1])*steps
                dz = (coord[coord2,2]-coord[coord1,2])*steps

                for j in range(nseg-1):
                    x1 = coord[coord1,0]+dx*j
                    y1 = coord[coord1,1]+dy*j
                    z1 = coord[coord1,2]+dz*j

                    x2 = coord[coord1,0]+dx*(j+1)
                    y2 = coord[coord1,1]+dy*(j+1)
                    z2 = coord[coord1,2]+dz*(j+1)

                    if np.any(np.isin(eq_els, i, assume_unique=True)) == True:
                        beam = v.Cylinder([[x1,y1,z1],[x2,y2,z2]],r=scale,res=4,c=f_color).alpha(alpha)
                    else:
                        beam = v.Cylinder([[x1,y1,z1],[x2,y2,z2]],r=scale,res=4,c=color).alpha(alpha)

                    if i in eq_dict:
                        beam.name = f"Beam element {i+1}, seg. {j+1}, Forces: [{eq_dict[i][0]}, {eq_dict[i][1]}, {eq_dict[i][2]}, {eq_dict[i][3]}]"
                    else:
                        beam.name = f"Beam element {i+1}, seg. {j+1}"
                    elements.append(beam)
            elif element_type == 5:
                if np.any(np.isin(eq_els, i, assume_unique=True)) == True:
                    beam = v.Cylinder([[coord[coord1,0],coord[coord1,1],coord[coord1,2]],[coord[coord2,0],coord[coord2,1],coord[coord2,2]]],r=scale,res=4,c=f_color).alpha(alpha)
                else:
                    beam = v.Cylinder([[coord[coord1,0],coord[coord1,1],coord[coord1,2]],[coord[coord2,0],coord[coord2,1],coord[coord2,2]]],r=scale,res=4,c=color).alpha(alpha)
                
                if i in eq_dict:
                    beam.name = f"Beam element {i+1}, Forces: [{eq_dict[i][0]}, {eq_dict[i][1]}, {eq_dict[i][2]}, {eq_dict[i][3]}]"
                else:
                    beam.name = f"Beam element {i+1}"

                elements.append(beam)

        if render_nodes == True:
            if element_type == 1 and spring == False:
                nodes = vdu.get_node_elements(coord,scale,alpha,dof,bcPrescr,bc,bc_color,fPrescr,f,f_color,dofs_per_node=1)
            elif element_type == 1:
                nodes = vdu.get_node_elements(coord,scale*0.5,alpha,dof,bcPrescr,bc,bc_color,fPrescr,f,f_color,dofs_per_node=1)
            elif element_type == 2:
                nodes = vdu.get_node_elements(coord,scale,alpha,dof,bcPrescr,bc,bc_color,fPrescr,f,f_color,dofs_per_node=3)
            elif element_type == 5:
                nodes = vdu.get_node_elements(coord,scale,alpha,dof,bcPrescr,bc,bc_color,fPrescr,f,f_color,dofs_per_node=6)
            plot_window.meshes[plot_window.fig].extend(elements)
            plot_window.nodes[plot_window.fig].extend(nodes)
        else:
            plot_window.meshes[plot_window.fig].extend(elements)

        return elements


    # Elements w/ a volume/surface (flow, solid & plate)
    elif element_type == 3 or element_type == 4 or element_type == 6:

        meshes = []
        nel = np.size(edof, axis = 0)

        coord[:] += offset

        for i in range(nel):
            eq_dict = {}
            indx = 0
            if isinstance(eq_els, np.ndarray):
                for j in eq_els:
                    eq_dict[j[0]] = eq[indx][0]
                    indx += 1

            if element_type == 3:
                coords = vdu.get_coord_from_edof(edof[i,:],dof,3)
            elif element_type == 4:
                coords = vdu.get_coord_from_edof(edof[i,:],dof,4)
            elif element_type == 6:
                coords = vdu.get_coord_from_edof(edof[i,:],dof,6)
                new_coord = np.zeros([8,3])

                new_coord[0,0] = coord[coords[0],0]
                new_coord[1,0] = coord[coords[1],0]
                new_coord[2,0] = coord[coords[2],0]
                new_coord[3,0] = coord[coords[3],0]
                new_coord[4,0] = coord[coords[0],0]
                new_coord[5,0] = coord[coords[1],0]
                new_coord[6,0] = coord[coords[2],0]
                new_coord[7,0] = coord[coords[3],0]

                new_coord[0,1] = coord[coords[0],1]
                new_coord[1,1] = coord[coords[1],1]
                new_coord[2,1] = coord[coords[2],1]
                new_coord[3,1] = coord[coords[3],1]
                new_coord[4,1] = coord[coords[0],1]
                new_coord[5,1] = coord[coords[1],1]
                new_coord[6,1] = coord[coords[2],1]
                new_coord[7,1] = coord[coords[3],1]

            if element_type == 3 or element_type == 4:
                if np.any(np.isin(eq_els, i, assume_unique=True)) == True:
                    mesh = v.Mesh([coord[coords,:],[[0,1,2,3],[4,5,6,7],[0,3,7,4],[1,2,6,5],[0,1,5,4],[2,3,7,6]]],alpha=alpha,c=f_color).lw(1)
                else:
                    mesh = v.Mesh([coord[coords,:],[[0,1,2,3],[4,5,6,7],[0,3,7,4],[1,2,6,5],[0,1,5,4],[2,3,7,6]]],alpha=alpha).lw(1)
            elif element_type == 6:
                if np.any(np.isin(eq_els, i, assume_unique=True)) == True:
                    mesh = v.Mesh([new_coord,[[0,1,2,3]]],alpha=alpha,c=f_color).lw(1)
                else:
                    mesh = v.Mesh([new_coord,[[0,1,2,3]]],alpha=alpha).lw(1)

            if element_type == 3:
                if i in eq_dict:
                    mesh.name = f"Flow element {i+1}, Force: {eq_dict[i][0]}"
                else:
                    mesh.name = f"Flow element {i+1}"
            elif element_type == 4:
                if i in eq_dict:
                    mesh.name = f"Solid element {i+1}, Force: {eq_dict[i][0]}"
                else:
                    mesh.name = f"Solid element {i+1}"
            elif element_type == 6:
                if i in eq_dict:
                    mesh.name = f"Plate element {i+1}, Force: {eq_dict[i][0]}"
                else:
                    mesh.name = f"Plate element {i+1}"
                
            meshes.append(mesh)

        if render_nodes == True:
            if element_type == 3:
                nodes = vdu.get_node_elements(coord,scale,alpha,dof,bcPrescr,bc,bc_color,fPrescr,f,f_color,dofs_per_node=1)
            elif element_type == 4 or element_type == 6:
                nodes = vdu.get_node_elements(coord,scale,alpha,dof,bcPrescr,bc,bc_color,fPrescr,f,f_color,dofs_per_node=3)

            plot_window.meshes[plot_window.fig].extend(meshes)
            plot_window.nodes[plot_window.fig].extend(nodes)
            #print("Adding mesh to figure ",plot_window.fig+1)
        else:
            plot_window.meshes[plot_window.fig].extend(meshes)
            #print("Adding mesh to figure ",plot_window.fig+1)
        return meshes



def draw_displaced_mesh(
    # Main input
    edof,
    coord,
    dof,
    element_type,
    a = None,
    scalars = None,
    
    # Other parameters
    scale = 0.02,
    alpha = 1,
    def_scale = 1,
    render_nodes = False,
    color = 'white',
    offset = [0, 0, 0],
    only_ret = False,
    lines = False,
    wireframe = False,

    # Input for colormapping
    colormap = 'jet',
    colors = 256,
    vmax = None,
    vmin = None,
    scalar_title = 'scalar',

    # Element-specific input
    spring = True,
    nseg = 2
    ):
    """
    Routine for displaced mesh for spring, bar, flow, solid, beam or plate elements.

    :param array edof: Element topology by degrees of freedom [nel x (n_dofs_per_node)|(n_dofs_per_node+1)*n_nodes ]
    :param array coord: Nodal coordinates [number of nodes x 3]
    :param array dof: Global degrees of freedom [number of nodes x degrees of freedom per node]
    :param int element_type: Element type [1-6]
    :param array a: Global displacement vector [degrees of freedom x 1]
    :param array scalars: Scalars for colormapping [number of elements x 1 | global number of nodes x 1 | global number of elements x number of nodes per element]
    :param float scale: Element scale, nodes are scaled 50% larger than this value
    :param float alpha: Element and node transparency [0-1]
    :param float def_scale: Deformation scalefactor
    :param bool render_nodes: If True, nodes are rendered
    :param str color: Element color
    :param list offset: Offset actors in 3D space [x, y, z]
    :param bool only_ret: Only return mesh, doesn't add it to a figure
    :param bool lines: Render lines at element borders
    :param bool wireframe: Only render lines at element borders
    :param str colormap: Name of colormap
    :param int colors: Number of colors in colormap [1-256]
    :param float vmax: Manually set colormap maximum
    :param float vmin: Manually set colormap minimum
    :param str scalar_title: Set scalar title (only for use when exporting)
    :param bool spring: If True, renders spring elements as coil springs
    :param int nseg: Number of points along beam elements for segmenting [number of segments + 1]

    :return array mesh: List of mesh actors (el. type 1/2/5) OR single mesh actor (el. type 3/4/6)
    """
    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window

    values = scalars

    if np.size(coord, axis = 1) == 1:
        coord = np.append(coord, np.zeros((np.size(coord, axis = 0),1)), axis=1)
        coord = np.append(coord, np.zeros((np.size(coord, axis = 0),1)), axis=1)
    elif np.size(coord, axis = 1) == 2:
        coord = np.append(coord, np.zeros((np.size(coord, axis = 0),1)), axis=1)

    if 1 <= element_type <= 6:
        if a is None and values is None:
            nel, ndof_per_el, nnode, ndim, ndof, ndof_per_n = vdu.check_input(edof,coord,dof,element_type,nseg=nseg)
            val = None
        elif a is None:
            nel, ndof_per_el, nnode, ndim, ndof, ndof_per_n, val = vdu.check_input(edof,coord,dof,element_type,values=values,nseg=nseg)
        elif values is None:
            nel, ndof_per_el, nnode, ndim, ndof, ndof_per_n, ndisp = vdu.check_input(edof,coord,dof,element_type,a,nseg=nseg)
            val = None
        else:
            nel, ndof_per_el, nnode, ndim, ndof, ndof_per_n, ndisp, val = vdu.check_input(edof,coord,dof,element_type,a,values,nseg=nseg)
    else:
        print("draw_displaced_mesh: Invalid element type, please declare 'element_type'. The element types are:\n    1 - Spring\n    2 - Bar\n    3 - Flow\n    4 - Solid\n    5 - Beam\n    6 - Plate")
        sys.exit()

    # OUTPUT FROM check_input
    # Number of elements:                       nel
    # Number of degrees of freedom per element: ndof_per_el
    # Number of nodes:                          nnode
    # Number of dimensions:                     ndim
    # Number of degrees of freedom:             ndof
    # Number of degrees of freedom per node:    ndof_per_n
    # Number of displacements:                  ndisp
    # Element/nodal values:                     val

    if a is None:
        a = np.zeros((nnode*ndof_per_n,1))
    elif element_type == 3:
        print('draw_displaced_mesh: Element type is flow, but deformation matrix given. Deformation was set to 0 for all DOFs')
        a = np.zeros((nnode*ndof_per_n,1))

    # Elements w/ a length (spring, bar & beam)
    if element_type == 1 or element_type == 2 or element_type == 5:
        ndof = np.size(dof, axis = 0)*np.size(dof, axis = 1)
        ncoord = np.size(coord, axis = 0)
        coord[:] += offset
        def_coord = np.zeros([ncoord,3])

        if vmin == None and vmax == None:
            vmin, vmax = np.min(values), np.max(values)
        elif vmin == None:
            vmin = np.min(values)
        elif vmax == None:
            vmax = np.max(values)

        for i in range(0, ncoord):
            if element_type == 1:
                a_dx = a[i]*def_scale

                x = coord[i,0]+a_dx
                y = coord[i,1]
                z = coord[i,2]

                def_coord[i] = [x,y,z]
            elif element_type == 2 or element_type == 5:
                a_dx, a_dy, a_dz = vdu.get_a_from_coord(i,6,a,def_scale)

                x = coord[i,0]+a_dx
                y = coord[i,1]+a_dy
                z = coord[i,2]+a_dz

                def_coord[i] = [x,y,z]

        nel = np.size(edof, axis = 0)
        elements = []
        res = 4

        el_values_array = np.zeros((1,4*4))[0,:]

        for i in range(nel):
            coord1,coord2 = vdu.get_coord_from_edof(edof[i,:],dof,element_type)

            if element_type == 1 and spring == True:
                element = v.Spring([def_coord[coord1,0],def_coord[coord1,1],def_coord[coord1,2]],[def_coord[coord2,0],def_coord[coord2,1],def_coord[coord2,2]],r=scale*1.5,c=color).alpha(alpha)
                element.name = f"Spring element {i+1}"
                elements.append(element)

                if values is not None:
                    el_values_array = []
                    for j in range(14):
                        el_values_array.append(values[i])
                    element.celldata[scalar_title] = el_values_array
                    element.cmap(colormap, scalar_title, on="cells", vmin=vmin, vmax=vmax)
                
            elif element_type == 2 or element_type == 1 and spring == False:
                element = v.Cylinder([[def_coord[coord1,0],def_coord[coord1,1],def_coord[coord1,2]],[def_coord[coord2,0],def_coord[coord2,1],def_coord[coord2,2]]],r=scale,res=4,c=color).alpha(alpha)
                if element_type == 2:
                    element.name = f"Bar element {i+1}"
                elif element_type == 1 and spring == False:
                    element.name = f"Spring element {i+1}"

                elements.append(element)

                if values is not None:
                    el_values_array = []
                    for j in range(6):
                        el_values_array.append(values[i])
                    element.celldata[scalar_title] = el_values_array
                    element.cmap(colormap, scalar_title, on="cells", n=colors, vmin=vmin, vmax=vmax)

            elif element_type == 5:
                if nseg > 2:
                    steps = np.float32(1/(nseg-1))

                    dx = (def_coord[coord2,0]-def_coord[coord1,0])*steps
                    dy = (def_coord[coord2,1]-def_coord[coord1,1])*steps
                    dz = (def_coord[coord2,2]-def_coord[coord1,2])*steps

                    for j in range(nseg-1):

                        x1 = def_coord[coord1,0]+dx*j
                        y1 = def_coord[coord1,1]+dy*j
                        z1 = def_coord[coord1,2]+dz*j

                        x2 = def_coord[coord1,0]+dx*(j+1)
                        y2 = def_coord[coord1,1]+dy*(j+1)
                        z2 = def_coord[coord1,2]+dz*(j+1)

                        element = v.Cylinder([[x1,y1,z1],[x2,y2,z2]],r=scale,res=res,c=color).alpha(alpha)
                        element.name = f"Beam element {i+1}, seg. {j+1}"
                        elements.append(element)

                        if values is not None:
                            el_value1 = values[nseg*i+j]
                            el_value2 = values[nseg*i+j+1]

                            el_values_array[1] = el_value1
                            el_values_array[3] = el_value1
                            el_values_array[5] = el_value1
                            el_values_array[7] = el_value1
                            el_values_array[12] = el_value1
                            el_values_array[13] = el_value1
                            el_values_array[14] = el_value1
                            el_values_array[15] = el_value1

                            el_values_array[0] = el_value2
                            el_values_array[2] = el_value2
                            el_values_array[4] = el_value2
                            el_values_array[6] = el_value2
                            el_values_array[8] = el_value2
                            el_values_array[9] = el_value2
                            el_values_array[10] = el_value2
                            el_values_array[11] = el_value2

                            element.pointdata[scalar_title] = el_values_array
                            
                            element.cmap(colormap, scalar_title, on="points", n=colors, vmin=vmin, vmax=vmax)

                else:
                    element = v.Cylinder([[def_coord[coord1,0],def_coord[coord1,1],def_coord[coord1,2]],[def_coord[coord2,0],def_coord[coord2,1],def_coord[coord2,2]]],r=scale,res=res,c=color).alpha(alpha)
                    element.name = f"Beam element {i+1}"
                    elements.append(element)

                    if values is not None:
                        el_value1 = values[2*i]
                        el_value2 = values[2*i+1]

                        el_values_array[1] = el_value1
                        el_values_array[3] = el_value1
                        el_values_array[5] = el_value1
                        el_values_array[7] = el_value1
                        el_values_array[12] = el_value1
                        el_values_array[13] = el_value1
                        el_values_array[14] = el_value1
                        el_values_array[15] = el_value1

                        el_values_array[0] = el_value2
                        el_values_array[2] = el_value2
                        el_values_array[4] = el_value2
                        el_values_array[6] = el_value2
                        el_values_array[8] = el_value2
                        el_values_array[9] = el_value2
                        el_values_array[10] = el_value2
                        el_values_array[11] = el_value2

                        element.pointdata[scalar_title] = el_values_array

                        element.cmap(colormap, scalar_title, on="points", n=colors, vmin=vmin, vmax=vmax)        

        if only_ret == False:

            if render_nodes == True:
                if element_type == 1 and spring == False:
                    nodes = vdu.get_node_elements(def_coord,scale,alpha,dof,dofs_per_node=1)
                elif element_type == 1:
                    nodes = vdu.get_node_elements(def_coord,scale*0.5,alpha,dof,dofs_per_node=1)
                elif element_type == 2:
                    nodes = vdu.get_node_elements(def_coord,scale,alpha,dof,dofs_per_node=3)
                elif element_type == 5:
                    nodes = vdu.get_node_elements(def_coord,scale,alpha,dof,dofs_per_node=6)
                plot_window.meshes[plot_window.fig].extend(elements)
                plot_window.nodes[plot_window.fig].extend(nodes)
            else:
                plot_window.meshes[plot_window.fig].extend(elements)

        return elements



    # Elements w/ a volume/surface (flow, solid & plate)
    elif element_type == 3 or element_type == 4 or element_type == 6:
        ex,ey,ez = cfc.coordxtr(edof,coord,dof)
        coord[:] += offset

        ed = cfc.extractEldisp(edof,a)
        if element_type == 3:
            if val != 'nodal_values_by_el':
                coord2, topo, node_dofs, a_node, node_scalars = vdu.convert_to_node_topo(edof,ex,ey,ez,ed,ignore_first=False,dofs_per_node=1)
            else:
                coord2, topo, node_dofs, a_node, node_scalars = vdu.convert_to_node_topo(edof,ex,ey,ez,ed,values,ignore_first=False,dofs_per_node=1)
            def_coord = coord2 + a_node*def_scale
        elif element_type == 4:
            if val != 'nodal_values_by_el':
                coord2, topo, node_dofs, a_node, node_scalars = vdu.convert_to_node_topo(edof,ex,ey,ez,ed,ignore_first=False)
            else:
                coord2, topo, node_dofs, a_node, node_scalars = vdu.convert_to_node_topo(edof,ex,ey,ez,ed,values,ignore_first=False)
            def_coord = coord2 + a_node*def_scale
        elif element_type == 6:
            coord2, topo, node_dofs, a_node, test = vdu.convert_to_node_topo(edof,ex,ey,ez,ed,ignore_first=False)
            def_coord = coord2
            def_coord[:,2] = a_node[:,0]*def_scale

        if element_type == 3 or element_type == 4:
            ct = vtk.VTK_HEXAHEDRON
            celltypes = [ct] * nel

            ug=v.UGrid([def_coord, topo, celltypes])
            ug.points(def_coord)
            
            mesh = ug.tomesh().alpha(alpha)

        elif element_type == 6:
            mesh = v.Mesh([def_coord, topo]).alpha(alpha)

        if lines == True:
            mesh.lw(1)

        if wireframe == True:
            mesh.wireframe()

        if vmin == None and vmax == None:
            vmin, vmax = np.min(values), np.max(values)
        elif vmin == None:
            vmin = np.min(values)
        elif vmax == None:
            vmax = np.max(values)


        if element_type == 3 or element_type == 4:
            if val and val == 'el_values':
                el_values = vdu.convert_el_values(edof,values)
                mesh.celldata[scalar_title] = el_values

                mesh.cmap(colormap, scalar_title, on="cells", n=colors, vmax=vmax, vmin=vmin)
            
            elif val and val == 'nodal_values_by_el':
                nodal_values = vdu.convert_nodal_values(edof,topo,dof,values)
                mesh.pointdata[scalar_title] = nodal_values
                mesh.cmap(colormap, scalar_title, on="points", n=colors, vmax=vmax, vmin=vmin)

            elif val and val == 'nodal_values':
                mesh.pointdata[scalar_title] = values
                mesh.cmap(colormap, scalar_title, on="points", n=colors, vmax=vmax, vmin=vmin)

        elif element_type == 6:
            if val and val == 'el_values':
                mesh.celldata[scalar_title] = values
                mesh.cmap(colormap, scalar_title, on="cells", n=colors, vmax=vmax, vmin=vmin)
        
        if only_ret == False:

            if render_nodes == True:
                if element_type == 3:
                    nodes = vdu.get_node_elements(def_coord,scale,alpha,dof,dofs_per_node=1)
                elif element_type == 4 or element_type == 6:
                    nodes = vdu.get_node_elements(def_coord,scale,alpha,dof,dofs_per_node=3)
                plot_window.meshes[plot_window.fig].extend(elements)
                plot_window.nodes[plot_window.fig].extend(nodes)
            else:
                plot_window.meshes[plot_window.fig].append(mesh)

        return mesh



def animation(
    # Main input
    edof,
    coord,
    dof,
    element_type,
    a=None,
    scalars=None,

    # Animation parameters
    steps=10,
    loop=False,
    negative=False,
    dt=50,
    animate_colormap=True,

    # Other parameters
    scale=0.02,
    alpha=1,
    def_scale=1,
    only_export = False,

    # Export
    export=False,
    file='anim/CALFEM_anim',

    # Input for colormapping
    colormap='jet',
    colors = 256,
    vmax=None,
    vmin=None,
    scalar_title = 'scalar',

    # Element-specific input
    spring = True,
    nseg=2
    ):
    """
    Routine for animating spring, bar, flow, solid, beam or plate elements.

    :param array edof: Element topology by degrees of freedom [nel x (n_dofs_per_node)|(n_dofs_per_node+1)*n_nodes ]
    :param array coord: Nodal coordinates [number of nodes x 3]
    :param array dof: Global degrees of freedom [number of nodes x degrees of freedom per node]
    :param int element_type: Element type [1-6]
    :param array a: Global displacement vector [degrees of freedom x 1]
    :param array scalars: Scalars for colormapping [number of elements x 1 | global number of nodes x 1 | global number of elements x number of nodes per element]
    :param int steps: Number of steps between undeformed & deformed states - 1
    :param bool loop: If True, runs a loop of animation (see documentation)
    :param bool negative: If True, runs a loop of animation & includes negative deformation (see documentation)
    :param int dt: Time difference between steps (milliseconds)
    :param bool animate_colormap: If True, animates the colormap on actors (NOTE: only supported for el. types 3/4/6)
    :param float scale: Element scale, nodes are scaled 50% larger than this value
    :param float alpha: Element and node transparency [0-1]
    :param float def_scale: Deformation scalefactor
    :param bool only_export: Only export mesh, doesn't add it to a figure
    :param str colormap: Name of colormap
    :param int colors: Number of colors in colormap [1-256]
    :param float vmax: Manually set colormap maximum
    :param float vmin: Manually set colormap minimum
    :param str scalar_title: Set scalar title (only for use when exporting)
    :param bool spring: If True, renders spring elements as coil springs
    :param int nseg: Number of points along beam elements for segmenting [number of segments + 1]
    """
    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window

    plot_window.loop = loop
    plot_window.dt = dt

    keyframe_dict = {}

    keyframes = {}

    plot_window.mode_anim[plot_window.fig] = True

    nnode = np.size(coord,1)

    values = scalars

    if element_type == 1 or element_type == 3:
        ndof_per_n = 1
    elif element_type == 2 or element_type == 4 or element_type == 6:
        ndof_per_n = 3
    elif element_type == 5:
        ndof_per_n = 6

    if a is None:
        a = np.zeros((nnode*ndof_per_n,1))
    elif element_type == 3:
        print('animation: NOTE Element type is flow, but deformation matrix given. Deformation was set to 0 for all DOFs')
        a = np.zeros((nnode*ndof_per_n,1))

    if negative == True:
        timesteps = np.arange(0, 1+1/(steps), 1/(steps))
        unique_timesteps = timesteps
        timesteps = np.append(timesteps, np.flip(np.arange(0, 1, 1/(steps))))
        timesteps = np.append(timesteps, np.flip(np.arange(-1, 0, 1/(steps))))
        unique_timesteps = np.append(unique_timesteps, np.flip(np.arange(-1, 0, 1/(steps))))
        timesteps = np.append(timesteps, np.arange(-1+1/(steps), 0, 1/(steps)))
    elif loop == True:
        timesteps = np.arange(0, 1+1/(steps), 1/(steps))
        unique_timesteps = timesteps
        timesteps = np.append(timesteps, np.flip(np.arange(0+1/(steps), 1, 1/(steps))))
    else:
        timesteps = np.arange(0, 1+1/(steps), 1/(steps))
        unique_timesteps = timesteps

    indx = 1
    for t in timesteps:
        keyframe_dict[indx] = np.round(t,3)
        indx += 1

    if values is not None:
        if vmax is None and vmin is None:
            vmin, vmax = np.min(values), np.max(values)
        elif vmax is None:
            vmax = np.max(values)
        elif vmin is None:
            vmin = np.min(values)

    it = 0
    for t in unique_timesteps:
        #print('Creating keyframe', it+1)
        mesh = None
        
        if values is not None and animate_colormap == True:
            
            if t >= 0:
                mesh = draw_displaced_mesh(edof,coord,dof,element_type,a*t,values*t,scale=scale,alpha=alpha,def_scale=def_scale,colormap=colormap,colors=colors,vmax=vmax,vmin=vmin,scalar_title=scalar_title,nseg=nseg,spring=spring,only_ret=True)
            
            else:
                mesh = draw_displaced_mesh(edof,coord,dof,element_type,a*t,-values*t,scale=scale,alpha=alpha,def_scale=def_scale,colormap=colormap,colors=colors,vmax=vmax,vmin=vmin,scalar_title=scalar_title,nseg=nseg,spring=spring,only_ret=True)
        
        elif values is not None and animate_colormap == False:
            mesh = draw_displaced_mesh(edof,coord,dof,element_type,a*t,values,scale=scale,alpha=alpha,def_scale=def_scale,colormap=colormap,colors=colors,vmax=vmax,vmin=vmin,scalar_title=scalar_title,nseg=nseg,spring=spring,only_ret=True)
        
        else:
            mesh = draw_displaced_mesh(edof,coord,dof,element_type,a*t,scale=scale,alpha=alpha,def_scale=def_scale,nseg=nseg,spring=spring,only_ret=True)
        
        if element_type == 1 or element_type == 2 or element_type == 5:
            
            mesh = v.merge(mesh)

        keyframes[np.round(t,3)] = mesh
        
        it += 1
    if only_export == False:
        plot_window.keyframes[plot_window.fig].append(keyframes)
        plot_window.keyframe_dict[plot_window.fig].append(keyframe_dict)

    if export == True:
        for key, val in keyframe_dict.items():
            mesh = keyframes[val]
            output = file+f'_{int(key)}'
            export_vtk(output,mesh)





### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
# Plotting functions, for adding usefull things to a rendering window


def add_scalar_bar(
    label,
    pos=[0.8,0.05],
    text_pos="bottom-right",
    font_size=24,
    color='black',
    on = 'mesh'
    ):
    """
    Routine for adding a scalar bar for a colormapped displaced mesh (requires draw_displaced_mesh to be run first).

    :param string label: Labeling of scalar bar (put scalar name and unit here)
    :param list pos: Positioning of sclarbar in the rendering window
    :param string text_pos: Positioning of label in the rendering window
    :param int font_size: Font size of the label
    :param string color: text color of label
    :param string on: scalarbar for mesh or vectors ['mesh'/'vectors']
    """
    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window

    fig = plot_window.fig

    if plot_window.mode_anim[fig] == True:

        for key in plot_window.keyframes[fig][0].keys():
            plot_window.keyframes[fig][0][key].addScalarBar(pos=pos, titleFontSize=font_size)

    else:

        if on == 'mesh':
            plot_window.meshes[fig][0].addScalarBar(pos=pos, titleFontSize=font_size)
        elif on == 'vectors':
            plot_window.vectors[fig][0].addScalarBar(pos=pos, titleFontSize=font_size)

    msg = v.Text2D(label, pos=text_pos, alpha=1, c=color)
    plot_window.msg[plot_window.fig] += [msg]



def add_text(
    text,
    color='black',
    pos='top-middle',
    size=1
    ):
    """
    Routine for adding 2D text to the render window.

    :param string text: Text string
    :param string pos: Positioning of text in the rendering window
    :param int size: Font size
    :param string color: text color
    """
    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window

    msg = v.Text2D(text, pos=pos, alpha=1, c=color)
    plot_window.msg[plot_window.fig] += [msg]



def add_text_3D(text,pos=(0,0,0),color='black',size=1,alpha=1):
    """
    Routine for adding text in 3D space.

    :param string text: Text string
    :param list pos: Positioning of text in 3D space
    :param int size: Font size
    :param string color: text color
    :param float alpha: transparency
    """
    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window

    msg = v.Text3D(text, pos=pos, s=size, font='', hspacing=1.15, vspacing=2.15, depth=0, italic=False, justify='bottom-left', c=color, alpha=alpha, literal=False)
    plot_window.msg[plot_window.fig] += [msg]



def add_mesh_axes(mesh):
    """
    Routine for automatically adding axes for a mesh to a figure.

    :param list mesh: mesh for adding axes
    """
    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window

    axes = mesh.buildAxes()
    plot_window.dia_axes[plot_window.fig].append(axes)



def add_axes(
    # Ranges
    xrange=[0,1],
    yrange=[0,1],
    zrange=[0,1],

    # Titles
    xtitle = 'x',
    ytitle = 'y',
    ztitle = 'z',
    htitle = '',

    # Grids
    xyGrid=True, yzGrid=True, zxGrid=True, xyGrid2=False, yzGrid2=False, zxGrid2=False,

    xyGridTransparent=True, yzGridTransparent=True, zxGridTransparent=True, xyGrid2Transparent=True, yzGrid2Transparent=True, zxGrid2Transparent=True,

    # Other
    numberOfDivisions=10
    ):
    """
    Routine for adding custom x-, y- & z-axes to a figure.

    :param list xrange: Range of x-axis [min,max]
    :param list yrange: Range of y-axis [min,max]
    :param list zrange: Range of z-axis [min,max]
    :param string xtitle: Label of x-axis
    :param string ytitle: Label of y-axis
    :param string ztitle: Label of z-axis
    :param boolean xyGrid: Add primary grid on xy-plane
    :param boolean yzGrid: Add primary grid on yz-plane
    :param boolean zxGrid: Add primary grid on xz-plane
    :param boolean xyGrid2: Add secondary grid on xy-plane
    :param boolean yzGrid2: Add secondary grid on yz-plane
    :param boolean zxGrid2: Add secondary grid on xz-plane
    :param boolean xyGridTransparent: Transparency for xyGrid
    :param boolean yzGridTransparent: Transparency for yzGrid
    :param boolean xzGridTransparent: Transparency for zxGrid
    :param boolean xyGrid2Transparent: Transparency for xyGrid2
    :param boolean yzGrid2Transparent: Transparency for yzGrid2
    :param boolean zxGrid2Transparent: Transparency for zxGrid2
    :param int numberOfDivisions: Number of divisions for all axes
    """
    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window

    axes = v.Axes(
        numberOfDivisions=numberOfDivisions,
        xtitle=xtitle,
        ytitle=ytitle,
        ztitle=ztitle,
        htitle=htitle,
        xyGrid=xyGrid,
        yzGrid=yzGrid,
        zxGrid=zxGrid,
        xyGrid2=xyGrid,
        yzGrid2=yzGrid,
        zxGrid2=zxGrid,
        xyGridTransparent=xyGridTransparent, yzGridTransparent=yzGridTransparent, zxGridTransparent=zxGridTransparent, xyGrid2Transparent=xyGrid2Transparent, yzGrid2Transparent=yzGrid2Transparent, zxGrid2Transparent=zxGrid2Transparent,
        xrange=xrange,
        yrange=yrange,
        zrange=zrange
    )

    plot_window.dia_axes[plot_window.fig].append(axes)



def add_projection(color='black',plane='xy',offset=-1,rulers=False):
    """
    Routine for adding a 2D projection of a model on a plane.

    :param string color: projection color
    :param string plane: which plane to project to
    :param float offset: offset projection from mesh (normal direction of plane)
    :param boolean rulers: add rulers to projection
    """
    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window

    if plane == 'xy':
        assem = plot_window.meshes[plot_window.fig]
        plot_window.meshes[plot_window.fig]
        proj = assem.projectOnPlane('z').z(offset).silhouette('2d').c(color)
        plot_window.proj[plot_window.fig] += [proj]
        if rulers == True:
            ruler = v.addons.RulerAxes(proj, xtitle='', ytitle='', ztitle='', xlabel='', ylabel='', zlabel='', xpad=1, ypad=1, zpad=1, font='Normografo', s=None, italic=0, units='m', c=color, alpha=1, lw=1, precision=3, labelRotation=0, axisRotation=0, xycross=True)
            plot_window.rulers[plot_window.fig] += [ruler]
    elif plane == 'xz':
        meshes = []
        for i in range(len(plot_window.meshes[plot_window.fig])):
            meshes.append(plot_window.meshes[plot_window.fig][i].clone())
        assem = v.merge(meshes)
        proj = assem.projectOnPlane('y').y(offset).silhouette('2d').c(color)
        plot_window.proj[plot_window.fig] += [proj]
        if rulers == True:
            ruler = v.addons.RulerAxes(proj, xtitle='', ytitle='', ztitle='', xlabel='', ylabel='', zlabel='', xpad=1, ypad=1, zpad=1, font='Normografo', s=None, italic=0, units='m', c=color, alpha=1, lw=1, precision=3, labelRotation=0, axisRotation=0, xycross=True)
            plot_window.rulers[plot_window.fig] += [ruler]
    elif plane == 'yz':
        assem = v.Assembly(plot_window.meshes[plot_window.fig])
        proj = assem.projectOnPlane('x').x(offset).silhouette('2d').c(color)
        plot_window.proj[plot_window.fig] += [proj]
        if rulers == True:
            ruler = v.addons.RulerAxes(proj, xtitle='', ytitle='', ztitle='', xlabel='', ylabel='', zlabel='', xpad=0, ypad=0.1, zpad=0.1, font='Normografo', s=None, italic=0, units='m', c=color, alpha=1, lw=1, precision=3, labelRotation=0, axisRotation=0, xycross=True)
            plot_window.rulers[plot_window.fig] += [ruler]
    else:
        print("Please choose a plane to project to. Set plane to 'xy', 'xz' or 'yz'")
        sys.exit()



def add_rulers(xtitle='', ytitle='', ztitle='', xlabel='', ylabel='', zlabel='', alpha=1):
    """
    Routine for automatically adding rulers (measurments) for a model.

    :param string xtitle: title for ruler along x-axis
    :param string ytitle: title for ruler along y-axis
    :param string ztitle: title for ruler along z-axis
    :param string xlabel: label for ruler along x-axis
    :param string ylabel: label for ruler along y-axis
    :param string zlabel: label for ruler along z-axis
    :param float alpha: transparency of rulers
    """
    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window

    assem = v.merge(plot_window.meshes[plot_window.fig])

    ruler = v.addons.RulerAxes(assem, xtitle=xtitle, ytitle=ytitle, ztitle=ztitle, xlabel=xlabel, ylabel=ylabel, zlabel=zlabel, xpad=0.1, ypad=0.1, zpad=0.1, font='Normografo', s=None, italic=0, units='m', c=(0,0,0), alpha=alpha, lw=1, precision=3, labelRotation=0, axisRotation=0, xycross=False)

    plot_window.rulers[plot_window.fig] += [ruler]



def eldia(ex,ey,ez,es,eci,dir='y',scale=1,thickness=5,alpha=1,label='y',invert=True):
    """
    Routine for beam section force diagrams.

    :param array ex: node x-coordinates
    :param array ey: node y-coordinates
    :param array ez: node z-coordinates
    :param array es: section force
    :param array eci: points along beam
    :param string dir: direction of y-axis for plot ['x'/'y'/'z']
    :param float scale: Scale y-axis of plot
    :param int thickness: Thickness of plot line (point diameter is x1.5)
    :param float alpha: transparency
    :param string label: label of plot y-axis
    :param boolean invert: invert y-axis of plot
    """
    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window

    nel = np.size(ex,0)
    nseg = np.size(eci,1)
    upd_scale = scale*(1/np.max(np.absolute(es)))

    for i in range(nel):
        x = np.linspace(ex[i,0], ex[i,1], nseg)
        if dir=='x':
            if invert == True:
                x = x - es[i]*upd_scale
            else:
                x = x + es[i]*upd_scale
        y = np.linspace(ey[i,0], ey[i,1], nseg)
        if dir=='y':
            if invert == True:
                y = y - es[i]*upd_scale
            else:
                y = y + es[i]*upd_scale
        z = np.linspace(ez[i,0], ez[i,1], nseg)
        if dir=='z':
            if invert == True:
                z = z - es[i]*upd_scale
            else:
                z = z + es[i]*upd_scale

        pts = []
        pts.append(v.Point(pos=[ex[i,0],ey[i,0],ez[i,0]], r=thickness*1.5, c='black', alpha=alpha))
        
        for j in range(nseg):
            pts.append(v.Point(pos=[x[j],y[j],z[j]], r=thickness*1.5, c='black', alpha=alpha))
        pts.append(v.Point(pos=[ex[i,1],ey[i,1],ez[i,1]], r=thickness*1.5, c='black', alpha=alpha))
        
        lines = []
        for j in range(len(pts)-1):
            lines.append(v.Lines(pts[j].points(), pts[j+1].points(), c='k4', alpha=alpha, res=2).lw(0.5*thickness))

        plot_window.dia_lines[plot_window.fig].append(lines)

        graph = v.merge(pts)
        for j in range(len(pts)):
            plot_window.dia_points[plot_window.fig].append(pts[j])
        
        if invert == True:
            ticks = -np.linspace(np.max(es[i])*upd_scale, np.min(es[i])*upd_scale, 10)
        else:
            ticks = np.linspace(np.min(es[i])*upd_scale, np.max(es[i])*upd_scale, 10)
        
        labels = np.round(np.linspace(np.min(es[i]), np.max(es[i]), 10),3)

        if dir=='x':
            if invert == True:
                axes = graph.buildAxes(xInverted=True, c='black', xTitleOffset=[(ticks[1]-ticks[0]),0,0], xtitle=label, xTitleRotation=270, xrange=[-np.min(es[i])*upd_scale+(ticks[1]-ticks[0]), -np.max(es[i])*upd_scale-(ticks[1]-ticks[0])], xValuesAndLabels=[(ticks[0], labels[0]), (ticks[1], labels[1]), (ticks[2], labels[2]), (ticks[3], labels[3]), (ticks[4], labels[4]), (ticks[5], labels[5]), (ticks[6], labels[6]), (ticks[7], labels[7]), (ticks[8], labels[8]), (ticks[9], labels[9])])
            else:
                axes = graph.buildAxes(c='black', xTitleOffset=[(ticks[1]-ticks[0]),0,0], xtitle=label, xTitleRotation=270, xrange=[np.min(es[i])*upd_scale, np.max(es[i])*upd_scale+(ticks[1]-ticks[0])], xValuesAndLabels=[(ticks[0], labels[0]), (ticks[1], labels[1]), (ticks[2], labels[2]), (ticks[3], labels[3]), (ticks[4], labels[4]), (ticks[5], labels[5]), (ticks[6], labels[6]), (ticks[7], labels[7]), (ticks[8], labels[8]), (ticks[9], labels[9])])
        elif dir=='y':
            if invert == True:
                axes = graph.buildAxes(yInverted=True, c='black', yTitleOffset=[0,(ticks[1]-ticks[0]),0], ytitle=label, yTitleRotation=270, yrange=[-np.max(es[i])*upd_scale-(ticks[1]-ticks[0]), -np.min(es[i])*upd_scale+(ticks[1]-ticks[0])], yValuesAndLabels=[(ticks[0], labels[0]), (ticks[1], labels[1]), (ticks[2], labels[2]), (ticks[3], labels[3]), (ticks[4], labels[4]), (ticks[5], labels[5]), (ticks[6], labels[6]), (ticks[7], labels[7]), (ticks[8], labels[8]), (ticks[9], labels[9])])
            else:
                axes = graph.buildAxes(c='black', yTitleOffset=[0,(ticks[1]-ticks[0]),0], ytitle=label, yTitleRotation=270, yrange=[np.min(es[i])*upd_scale, np.max(es[i])*upd_scale+(ticks[1]-ticks[0])], yValuesAndLabels=[(ticks[0], labels[0]), (ticks[1], labels[1]), (ticks[2], labels[2]), (ticks[3], labels[3]), (ticks[4], labels[4]), (ticks[5], labels[5]), (ticks[6], labels[6]), (ticks[7], labels[7]), (ticks[8], labels[8]), (ticks[9], labels[9])])
        elif dir=='z':
            if invert == True:
                axes = graph.buildAxes(zInverted=True, c='black', zTitleOffset=[0,0,(ticks[1]-ticks[0])], ztitle=label, zTitleRotation=270, zrange=[-np.min(es[i])*upd_scale+(ticks[1]-ticks[0]), -np.max(es[i])*upd_scale-(ticks[1]-ticks[0])], zValuesAndLabels=[(ticks[0], labels[0]), (ticks[1], labels[1]), (ticks[2], labels[2]), (ticks[3], labels[3]), (ticks[4], labels[4]), (ticks[5], labels[5]), (ticks[6], labels[6]), (ticks[7], labels[7]), (ticks[8], labels[8]), (ticks[9], labels[9])])
            else:
                axes = graph.buildAxes(c='black', zTitleOffset=[0,0,(ticks[1]-ticks[0])], ztitle=label, zTitleRotation=270, zrange=[np.min(es[i])*upd_scale, np.max(es[i])*upd_scale+(ticks[1]-ticks[0])], zValuesAndLabels=[(ticks[0], labels[0]), (ticks[1], labels[1]), (ticks[2], labels[2]), (ticks[3], labels[3]), (ticks[4], labels[4]), (ticks[5], labels[5]), (ticks[6], labels[6]), (ticks[7], labels[7]), (ticks[8], labels[8]), (ticks[9], labels[9])])

        plot_window.dia_axes[plot_window.fig].append(axes)



def elprinc(ex,ey,ez,val,vec,ed=None,scale=.1, colormap = 'jet', unit='Pa'):
    """
    Routine for principal stresses.

    :param array ex: node x-coordinates
    :param array ey: node y-coordinates
    :param array ez: node z-coordinates
    :param array val: eigenvalues
    :param array vec: eigenvectors
    :param array ed: element displacements
    :param float scale: vector scale
    :param string colormap: Colormap for vector scalar
    :param string unit: vector unit
    """
    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window

    nel = np.size(ex,0)
    n = int((np.size(val,0)*np.size(val,1))/np.size(ex,0))

    upd_scale = (1/np.max(val))*scale

    if ed is None:
        ed = np.zeros((nel,3))

    pts = []
    points = []
    vectors = []
    text = []
    vmin = np.min(val)
    vmax = np.max(val)
    values = []
    for i in range(nel):
        x = np.average(ex[i])+ed[i,0]
        y = np.average(ey[i])+ed[i,1]
        z = np.average(ez[i])+ed[i,2]
        for j in range(n):
            points.append([x,y,z])
            text.append(f'Principal stress {j+1} at El. {i+1}: {np.round(val[i,j],3)} {unit}')
            values.append([val[i,j],val[i,j],val[i,j],val[i,j],val[i,j],val[i,j]])
            vector = vec[i, :, j]*val[i,j]*upd_scale
            vectors.append(vector)

    plot = vdu.vectors(points, vectors, c='k', alpha=1, shaftLength=0.8, shaftWidth=0.05, headLength=0.25, headWidth=0.2, fill=True, text=text, vmax=vmax, vmin=vmin, cmap=colormap, values=values)

    plot_window.vectors[plot_window.fig].extend(plot)



def elflux(ex,ey,ez,vec,ed=None,scale=.1, colormap = 'jet', unit='W/m^2'):
    """
    Routine for flux (or other) vectors.

    :param array ex: node x-coordinates
    :param array ey: node y-coordinates
    :param array ez: node z-coordinates
    :param array vec: vectors
    :param array ed: element displacements
    :param float scale: vector scale
    :param string colormap: Colormap for vector scalar
    :param string unit: vector unit
    """
    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window

    nel = np.size(ex,0)
    upd_scale = (1/np.max(vec))*scale

    if ed is None:
        ed = np.zeros((nel,3))

    points = []
    vectors = []
    text = []
    vmin = np.min(vec)
    vmax = np.max(vec)
    values = []
    for i in range(nel):
        flux_tot = np.sqrt(vec[i,0]**2 + vec[i,1]**2 + vec[i,2]**2)
        points.append([np.average(ex[i])+ed[i,0], np.average(ey[i])+ed[i,1], np.average(ez[i])+ed[i,2]])
        text.append(f'Flux at El. {i+1}: {np.round(flux_tot,3)} {unit}')
        values.append([flux_tot,flux_tot,flux_tot,flux_tot,flux_tot,flux_tot])
        vectors.append(vec[i]*upd_scale)

    plot = vdu.vectors(points, vectors, c='k', alpha=1, shaftLength=0.8, shaftWidth=0.05, headLength=0.25, headWidth=0.2, fill=True, text=text, vmax=vmax, vmin=vmin, cmap=colormap, values=values)

    plot_window.vectors[plot_window.fig].extend(plot)





### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
# Functions for importing/exporting


def import_mat(file,list=None):
    """
    Routine for importing from MATLAB.

    :param str file: File name for import (with or without .mat-file type)
    :param list list: List of variables to import, and the order of them

    :yield x: matlab variables
    """
    data = {} # dict to be returned by loadmat
    loadmat(file, data, variable_names=list)

    if list == None: # remove unnecessary entries
        keys = ['__header__', '__version__', '__globals__']
        for key in keys:
            data.pop(key)
        return data # returns the data, random ordering
    else: # supplying a 'list' is recommended
        ret = dict.fromkeys(list,None) # returns the data, ordering by user
        for key,val in ret.items():
            x = data[key] 
            if key == 'edof' or key == 'Edof' or key == 'EDOF':
                x = np.delete(x,0,1) # auto convert Edof from Matlab to Python
            yield x # returns the data one-by-one, ordering by user



def export_vtk(file, meshes):
    """
    Routine for exporting to VTK (ParaView).

    :param str file: File name for export (without .vtk-file type)
    :param list meshes: List of mesh actors or single mesh actor
    """
    if isinstance(meshes, list):
        mesh = v.merge(meshes)
    else:
        mesh = meshes
    
    v.io.write(mesh, file+".vtk")





### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
# Functions for handling rendering


def figure(fig,bg='white',flat=False,hover=False):
    """
    Routine for choosing what figure to add objects to, and creating necessary lists etc.

    :param int fig: Figure number
    :param str bg: Background color of figure. If set to 'black', some objects in figure will automatically be white
    :param bool flat: 2D mode
    :param bool hover: Hover mode (highlights actors and outputs node/element info by hovering instead of clicking)
    """
    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window

    if fig < 1:
        print("figure: Please give a positive integer (> 0)")
        sys.exit()
    else:
        plot_window.fig = fig - 1

    #print("Selecting figure ",fig)
    if fig > 1:
        while len(plot_window.meshes) < plot_window.fig + 1:
            plot_window.geometries.append([])
            plot_window.meshes.append([])
            plot_window.nodes.append([])

            plot_window.msg.append([])
            plot_window.proj.append([])
            plot_window.rulers.append([])
            plot_window.vectors.append([])

            plot_window.dia_lines.append([])
            plot_window.dia_points.append([])
            plot_window.dia_axes.append([])

            plot_window.keyframes.append([])
            plot_window.keyframe_dict.append([])

        while len(plot_window.mode_2D) < plot_window.fig + 1:
            plot_window.bg.append(bg)
            plot_window.mode_2D.append(False)
            plot_window.mode_hover.append(False)
            plot_window.mode_anim.append(False)

    elif fig == 1:
        plot_window.bg.append(bg)
        plot_window.mode_2D.append(False)
        plot_window.mode_hover.append(False)
        plot_window.mode_anim.append(False)

    plot_window.bg[fig-1] = bg

    if flat == True:
        plot_window.mode_2D[fig-1] = True

    if hover == True:
        plot_window.mode_hover[fig-1] = True



def show_and_wait():
    """
    Routine for starting visualization & showing figure(s).
    """
    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window
    plot_window.render()
