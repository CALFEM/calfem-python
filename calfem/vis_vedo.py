# -*- coding: utf-8 -*-
"""
CALFEM VTK

Contains all the functions for 3D visualization in CALFEM using VTK

@author: Andreas Åmand
"""

import numpy as np
import vedo as v
#import polyscope as p
#import dolfin
from vedo.dolfin import dataurl, download, plot
import sys
from PyQt5 import Qt
from PyQt5 import uic
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import webbrowser


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
# Tools

class tools:
    def get_coord_from_edof(edof_row,dof,element_type):
        if element_type == 1 or element_type == 2 or element_type == 5:
            edof_row1,edof_row2 = np.split(edof_row,2)
            coord1 = int(np.where(np.all(edof_row1==dof,axis=1))[0])
            coord2 = int(np.where(np.all(edof_row2==dof,axis=1))[0])
            return coord1, coord2
        elif element_type == 3 or element_type == 4:
            edof_row1,edof_row2,edof_row3,edof_row4,edof_row5,edof_row6,edof_row7,edof_row8 = np.split(edof_row,8)
            coord1 = int(np.where(np.all(edof_row1==dof,axis=1))[0])
            coord2 = int(np.where(np.all(edof_row2==dof,axis=1))[0])
            coord3 = int(np.where(np.all(edof_row3==dof,axis=1))[0])
            coord4 = int(np.where(np.all(edof_row4==dof,axis=1))[0])
            coord5 = int(np.where(np.all(edof_row5==dof,axis=1))[0])
            coord6 = int(np.where(np.all(edof_row6==dof,axis=1))[0])
            coord7 = int(np.where(np.all(edof_row7==dof,axis=1))[0])
            coord8 = int(np.where(np.all(edof_row8==dof,axis=1))[0])
            coords = np.array([coord1, coord2, coord3, coord4, coord5, coord6, coord7, coord8])
            return coords
        elif element_type == 6:
            print('lol plattor stöds inte än')

    def get_a_from_coord(coord_row_num,num_of_deformations,a,scale):
        dx = a[coord_row_num*num_of_deformations]*scale
        dy = a[coord_row_num*num_of_deformations+1]*scale
        dz = a[coord_row_num*num_of_deformations+2]*scale
        return dx, dy, dz

    def get_node_elements(coord,scale,alpha):
        nnode = np.size(coord, axis = 0)
        ncoord = np.size(coord, axis = 1)
        nodes = []
        for i in range(nnode):
            if ncoord == 3:
                node = v.Sphere(c='white').scale(1.5*scale).pos([coord[i,0],coord[i,1],coord[i,2]]).alpha(alpha)
                node.name = f"Node nr.{i}"
                nodes.append(node)
            elif ncoord == 2:
                node = v.Sphere(c='white').scale(1.5*scale).pos([coord[i,0],coord[i,1],0]).alpha(alpha)
                node.name = f"Node nr.{i}"
                nodes.append(node)
            elif ncoord == 1:
                node = v.Sphere(c='white').scale(1.5*scale).pos([coord[i,0],0,0]).alpha(alpha)
                node.name = f"Node nr.{i}"
                nodes.append(node)
        return nodes


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
# Initialization

def init_app():
    global vedo_app
    #global anim
    vedo_app = Qt.QApplication.instance()
    #anim = animations()

    if vedo_app is None:
        print("No QApplication instance found. Creating one.")
        # if it does not exist then a QApplication is created
        vedo_app = Qt.QApplication(sys.argv)
    else:
        print("QApplication instance found.")

    return vedo_app

class VedoPlotWindow:
    __instance = None
    @staticmethod
    def intance():
        """ Static access method. """
        if VedoPlotWindow.__instance == None:
            VedoPlotWindow()

        return VedoPlotWindow.__instance

    def __init__(self):
        """ Virtually private constructor. """
        if VedoPlotWindow.__instance != None:
            raise Exception("This class is a singleton!")
        else:
            VedoPlotWindow.__instance = self
            self.plot_window = VedoMainWindow()


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
# Main window

class VedoMainWindow(Qt.QMainWindow):
#class VedoMainWindow():
    
    def __init__(self):
        Qt.QMainWindow.__init__(self)
        #self.animations = animations()
        #self.initialize()
        
    #def initialize(self):
            
        #uic.loadUi("../QtVTKMainWindow.ui", self)
        #uic.loadUi("QtVTKMainWindow.ui", self)
        self.frame = Qt.QFrame()
        self.layout = Qt.QVBoxLayout()
        self.vtkWidget = QVTKRenderWindowInteractor(self.frame)

        # Create renderer and add the vedo objects and callbacks
        #self.plotter = v.Plotter(qtWidget=self.vtkWidget,bg="black",bg2="blackboard",axes=4).show(viewup=[0,0,1], rate=30)
        self.plotter = v.Plotter(qtWidget=self.vtkWidget,bg="black")#.show(viewup=[0,0,1],rate=30,interactive=True)
        #self.plotter = v.Plotter(qtWidget=self.vtkWidget)


        
        #Test click

        self.plotter.addCallback('mouse click', self.click)
        #self.plotter.addCallback('MouseClick', self.click)

        self.silcont = [None]
        self.click_msg = v.Text2D("", pos="bottom-center", c='k', bg='r9', alpha=0.8)

        #self.plotter.show(self.click_msg).close()
        #self.plotter.show(self.click_msg)
        #self.plotter += self.click_msg
        self.plotter.add(self.click_msg)

        
        #self.plotter.show(viewup=[0,0,1],rate=30,mode=10)
        #self.plotter.show(mode=10,axes=4)
        # show(*actors, at=None, axes=None, resetcam=None, zoom=False, interactive=None, viewup='', azimuth=0, elevation=0, roll=0, camera=None, interactorStyle=0, mode=None, rate=None, bg=None, bg2=None, size=None, title=None, q=False)
        #self.plotter.addLegendBox()
        #axes=4 -> koordinataxlar
        #axes=4 -> koordinataxlar i hörnet
        #axes=7 -> mått
        #axes=12 -> gradskiva
        self.layout.addWidget(self.vtkWidget)
        self.frame.setLayout(self.layout)
        self.setCentralWidget(self.frame)

        #self.plotter.show(mode=1)

        self.show()
        
        #self.plotter.show(axes=4, viewup=[1,1,1], rate=30, resetcam=True)

        #v.interactive().show()
        #self.show()


        # Animation
        #self.anim_start = 0.0
        #self.anim_end = 1
        self.anim_steps = 101.0
        self.anim_coord = None
        self.anim_def_coord = None
        self.anim_element_type = None
        """
        # QtButtons
        self.actionShow_Axis.triggered.connect(lambda: self.show_axis()) #Show Axis
        self.actionCALFEM_for_Python_documentation.triggered.connect(lambda: webbrowser.open('https://calfem-for-python.readthedocs.io/en/latest/')) #Documentation
        self.actionShow_browser.triggered.connect(lambda: self.obj_browser())
        self.actionShow_Grid.triggered.connect(lambda: self.show_grid())

        # Camera
        self.actionReset_Camera.triggered.connect(lambda: self.reset_camera()) #Reset camera
        self.actionActor.triggered.connect(lambda: self.actor_mode())
        self.actionDefault.triggered.connect(lambda: self.default_mode())
        self.actionTrackball.triggered.connect(lambda: self.trackball_mode())
        self.actionSelection.triggered.connect(lambda: self.select_mode())

        # Animation
        self.actionStart_Animation.triggered.connect(lambda: self.start_animation())
        self.actionEdit_Parameters.triggered.connect(lambda: self.edit_parameters())
        """






        # View parameters
        self.axis = True
        self.grid = False
        self.mode = 'default' # def=10, trackball=0, actor=1




        """
        # Click a sphere to highlight it
        from vedo import Text2D, Sphere, Plotter
        import numpy as np

        spheres = []
        for i in range(25):
            p = np.random.rand(2)
            s = Sphere(r=0.05).pos(p).color('k5')
            s.name = f"sphere nr.{i} at {p}"
            spheres.append(s)

        def func(evt):
            if not evt.actor: return
            sil = evt.actor.silhouette().lineWidth(6).c('red5')
            msg.text("You clicked: "+evt.actor.name)
            plt.remove(silcont.pop()).add(sil)
            silcont.append(sil)

        silcont = [None]
        msg = Text2D("", pos="bottom-center", c='k', bg='r9', alpha=0.8)

        plt = Plotter(axes=1, bg='black')
        plt.addCallback('mouse click', func)
        plt.show(spheres, msg, __doc__, zoom=1.2).close()
        """




        """
        # Generate the silhouette of a mesh
        # as seen along a specified direction
        
        from vedo import *

        s = Hyperboloid().rotateX(20)

        sx = s.clone().projectOnPlane('x').c('r').x(-3) # sx is 2d
        sy = s.clone().projectOnPlane('y').c('g').y(-3)
        sz = s.clone().projectOnPlane('z').c('b').z(-3)

        show(s,
             sx, sx.silhouette('2d'), # 2d objects dont need a direction
             sy, sy.silhouette('2d'),
             sz, sz.silhouette('2d'),
             __doc__,
             axes={'zxGrid':True, 'yzGrid':True},
             viewup='z',
        ).close()
        """








    






    def click(self,evt):
        if not evt.actor:
            print('not actor')
            return
        sil = evt.actor.silhouette().lineWidth(6).c('red5')
        self.click_msg.text(evt.actor.name)
        self.plotter.remove(self.silcont.pop()).add(sil)
        self.silcont.append(sil)







    def render_geometry(self,meshes,nodes=None,merge=False):
        
        

        # Mesh/elements plotted, possibly merged for correct numbering
        if merge == True:
            mesh = v.merge(meshes,flag=True)
            self.plotter += mesh
            #self.plotter.add(mesh)

            mesh.clean()
            #mesh.computeNormals().clean().lw(0.1)
            pids = mesh.boundaries(returnPointIds=True)
            bpts = mesh.points()[pids]

            pts = v.Points(bpts, r=1, c='red')
            labels = mesh.labels('id', scale=0.02).c('w')

            self.plotter += pts
            self.plotter += labels
        else:
            nel = np.size(meshes, axis = 0)
            for i in range(nel):
                self.plotter += meshes[i]
                #self.plotter.add(meshes[i])

        # Optional nodes
        if nodes is not None:
            nnode = np.size(nodes, axis = 0)
            for i in range(nnode):
                self.plotter += nodes[i]
                #self.plotter.add(nodes[i])

        #self.plotter.render(resetcam=True)
        self.plotter.addHoverLegend().show(resetcam=True,axes=4)
        #self.plotter.interactive()
        #self.show()

        


    def render_addon(self,addon):
        self.plotter += addon
        self.plotter.show(resetcam=True,axes=4)


    """

    def render_beam_geometry(self,nodes,elements):
        
        nnode = np.size(nodes, axis = 0)
        for i in range(nnode):
            self.plotter += nodes[i]
            #self.plotter.add(nodes[i])
        #self.plotter.add(nodes,resetcam=True)
        #self.plotter.add(elements,resetcam=True)

        nel = np.size(elements, axis = 0)
        for i in range(nel):
            self.plotter += elements[i]
            #self.plotter.add(elements[i])

        #self.plotter.render(resetcam=True)
        self.plotter.show(resetcam=True,mode=10,axes=4)
        #self.plotter.interactive()

        #self.plotter.addButton(
        #    self.show_hide(elements,nodes),
        #    pos=(0.9,0.95),
        #    states=['Hide undeformed mesh','Show undeformed mesh'],
        #    c=['w','w'],
        #    bc=['blackboard','blackboard'],
        #    font='courier',
        #    size=8
        #    )

    def render_solid_geometry(self,meshes):

        #self.plotter.add(mesh,resetcam=True)

        nel = np.size(meshes, axis = 0)
        mesh = v.merge(meshes,flag=True)
        self.plotter += mesh

        mesh.computeNormals().clean().lw(0.1)
        pids = mesh.boundaries(returnPointIds=True)
        bpts = mesh.points()[pids]

        pts = v.Points(bpts, r=1, c='red')
        labels = mesh.labels('id', scale=0.02).c('w')

        self.plotter += pts
        self.plotter += labels

        self.plotter.show(resetcam=True,mode=10,axes=4)

        #for i in range(nel):
        #    self.plotter += mesh[i]
        #    #self.plotter.add(mesh[i])
        #    #ps_mesh = p.register_surface_mesh("My vedo mesh",solid_mesh.points(),solid_mesh.faces(),color=[0.5,0,0],smooth_shade=False)
        #    #mesh[i].computeNormals().clean().lw(0.1)
        #
        #    pids = mesh[i].boundaries(returnPointIds=True)
        #    bpts = mesh[i].points()[pids]
        #
        #    pts = v.Points(bpts, r=1, c='red')
        #    labels = mesh[i].labels('id', scale=0.02).c('w')
        #    #print(labels)
        #    #print(pts)
        #    #print(bpts)
        #    #print(pids)
        #
        #    self.plotter += pts
        #    self.plotter += labels
        #
        #    self.plotter.show(resetcam=True,mode=10,axes=4)
        #
        #    #self.plotter.show(mesh[i],pts,labels,resetcam=True,mode=10,axes=4)
        #    #self.plotter.show(mesh[i],resetcam=True,mode=10,axes=4)

        #self.plotter.show(resetcam=True,mode=10,axes=4)
        #self.plotter.show(resetcam=True,viewup=[0,0,1],rate=30,interactive=True,mode=10)
        #self.plotter.interactive()
        #ps_mesh.add_scalar_quantity("heights", solid_mesh.points()[:,2], defined_on='vertices')
        #p.show()
        """














    # Functions run from QtButtons

    def reset_camera(self):
        print("reset camera")
        # återställer zoom & position
        #self.plotter.resetCamera()
        self.plotter.show(resetcam=True)
        #self.plotter.show(resetcam=True,viewup=[0,0,1],rate=30,interactive=True,mode=10, at=0)
        #self.plotter.resetCamera()
        #self.iren.Render()

    def show_axis(self):
        print("show axis")
        #self.plotter(axes=13)
        #self.plotter.addGlobalAxes(axtype=4)
        #self.plotter.render()
        #self.plotter.axes(0)
        if self.axis == True:
            self.axis = False
            self.plotter.show(axes=0)
            #self.plotter.addGlobalAxes(0)
        else:
            self.axis = True
            self.plotter.show(axes=4)
            #self.plotter.addGlobalAxes(4)
        #self.plotter.render()
        #self.show()

        print(self.axis)

    def show_grid(self):
        grid_parameters = dict(
            #xtitle='My variable \Omega^\lowerxi_lm  in units of \mum^3', # latex-style syntax
            #ytitle='This is my highly\ncustomized y-axis',
            #ztitle='z in units of Å', # many unicode chars are supported (type: vedo -r fonts)
            yValuesAndLabels=[(-3.2,'Mark^a_-3.2'), (-1.2,'Carmen^b_-1.2'), (3,'John^c_3')],
            textScale=1.3,       # make all text 30% bigger
            numberOfDivisions=5, # approximate number of divisions on longest axis
            axesLineWidth= 2,
            gridLineWidth= 1,
            zxGrid2=True,        # show zx plane on opposite side of the bounding box
            yzGrid2=True,        # show yz plane on opposite side of the bounding box
            xyPlaneColor='green7',
            xyGridColor='dg',    # darkgreen line color
            xyAlpha=0.2,         # grid opacity
            xTitlePosition=0.5,  # title fractional positions along axis
            xTitleJustify="top-center", # align title wrt to its axis
            yTitleSize=0.02,
            yTitleBox=True,
            yTitleOffset=0.05,
            yLabelOffset=0.4,
            yHighlightZero=True, # draw a line highlighting zero position if in range
            yHighlightZeroColor='red',
            zLineColor='blue',
            zTitleColor='blue',
            zTitleBackfaceColor='v', # violet color of axis title backface
            labelFont="Quikhand",
            yLabelSize=0.025,    # size of the numeric labels along Y axis
            yLabelColor='dg',    # color of the numeric labels along Y axis
            )
        #grid = v.base.BaseActor()
        #self.plotter += grid
        #self.plotter.add(grid, at=0)

        #self.plotter.show(resetcam=True,viewup=[0,0,1],rate=30,interactive=True,mode=10,axes=dict(zLabelSize=.04, numberOfDivisions=10))
        #self.plotter.show(axes=9, at=0)
        self.plotter.addGlobalAxes(grid_parameters)
        self.plotter.render()

    def obj_browser(self):
        meshes = self.plotter.getMeshes()
        #meshes = self.plotter.getVolumes()
        v.applications.Browser(meshes, prefix='Mesh No. ')

    def start_animation(self):
        print('start anim')
        animate(self.anim_steps,self.anim_coord,self.anim_def_coord,self.anim_element_type)

    def edit_parameters(self):
        print('edit parameters')
        #edit_anim()

    def trackball_mode(self):
        print('actor')
        self.plotter.show(mode=0)

    def default_mode(self):
        print('actor')
        self.plotter.show(mode=10)

    def actor_mode(self):
        print('actor')
        self.plotter.show(mode=1)

    def select_mode(self):
        print('actor')
        self.plotter.show(mode=6)

    def show_hide(elements=None,nodes=None,mesh=None):
        if mesh == None:
            elements.alpha(0)
            self.plotter.render()




### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###


def draw_geometry(
    edof,
    coord,
    dof,
    element_type,
    el_values=None,
    label=None,
    colormap='jet',
    scale=0.02,
    alpha=1,
    nseg=2,
    render_nodes=True,
    color=None
    #export=None
    ):



    # Element types: 1: Spring, 2: Bar, 3: Flow, 4: Solid, 5: Beam, 6: Plate

    # 2 node: 1,2,5 (1-3D)
    # 3 node: 3,4 (Triangular 2D)
    # 4 node: 3,4,6 (Quadratic/rectangular/isoparametric 2D)
    # 8 node: 3,4 (Isoparametric 2D or 3D)

    app = init_app()
    plot_window = VedoPlotWindow.intance().plot_window
    #animations = animations().instance()

    if element_type == 1:

        #ndof = np.size(dof, axis = 0)*np.size(dof, axis = 1)
        #nnode = np.size(coord, axis = 0)
        #ncoord = np.size(coord, axis = 1) # Använd för att visualisera 1-3 dim.
        #nodes = []

        #for i in range(nnode):
        #    nodes.append(v.Sphere(c='white').scale(1.5*scale).pos([coord[i,0],0,0]).alpha(alpha))            

        nel = np.size(edof, axis = 0)
        elements = []

        for i in range(nel):
            coord1,coord2 = tools.get_coord_from_edof(edof[i,:],dof,element_type)

            #print(coord[coord1,0])
            #print(coord[coord2,0])
            #spring = v.Spring([coord[coord1,0],0,0],[coord[coord2,0],0,0])
            spring = v.Spring([coord[coord1,0],0,0],[coord[coord2,0],0,0],r=1.5*scale).alpha(alpha)
            elements.append(spring)

        #print(elements,nodes)
        if render_nodes == True:
            nodes = tools.get_node_elements(coord,scale,alpha)
            plot_window.render_geometry(elements,nodes)
        else:
            plot_window.render_geometry(elements)


    elif element_type == 2:
        print("bar")
        if render_nodes == True:
            nodes = tools.get_node_elements(coord,scale,alpha)

        nel = np.size(edof, axis = 0)
        elements = []

        res = 4

        vmin, vmax = np.min(el_values), np.max(el_values)

        el_values_array = np.zeros((1,4*res))[0,:]

        for i in range(nel):
            coord1,coord2 = tools.get_coord_from_edof(edof[i,:],dof,element_type)

            tube = v.Cylinder([[coord[coord1,0],coord[coord1,1],coord[coord1,2]],[coord[coord2,0],coord[coord2,1],coord[coord2,2]]],r=scale,res=res,c=color).alpha(alpha)
            tube.name = f"Bar nr.{i}"
            elements.append(tube)

            if el_values is not None:
                #el_value1 = el_values[2*i,:]
                #el_value2 = el_values[2*i+1,:]
                #el_value = el_values[i]

                el_values_array[1] = el_values[i]
                el_values_array[3] = el_values[i]
                el_values_array[5] = el_values[i]
                el_values_array[7] = el_values[i]
                el_values_array[12] = el_values[i]
                el_values_array[13] = el_values[i]
                el_values_array[14] = el_values[i]
                el_values_array[15] = el_values[i]

                el_values_array[0] = el_values[i]
                el_values_array[2] = el_values[i]
                el_values_array[4] = el_values[i]
                el_values_array[6] = el_values[i]
                el_values_array[8] = el_values[i]
                el_values_array[9] = el_values[i]
                el_values_array[10] = el_values[i]
                el_values_array[11] = el_values[i]

                tube.cmap(colormap, el_values_array, on="points", vmin=vmin, vmax=vmax)
                #tube.addScalarBar(label)

            #print(coord[coord1,0])
            #print(coord[coord2,0])
            #spring = v.Spring([coord[coord1,0],0,0],[coord[coord2,0],0,0])
            #tube = v.Cylinder([[coord[coord1,0],coord[coord1,1],coord[coord1,2]],[coord[coord2,0],coord[coord2,1],coord[coord2,2]]],r=scale,res=res,c='white').alpha(alpha)
            #elements.append(tube)

        #print(elements,nodes)
        if render_nodes == True:
            plot_window.render_geometry(elements,nodes)
        else:
            plot_window.render_geometry(elements)


        return elements

    elif element_type == 3:
        print("flow")

    elif element_type == 4:

        #nnode = np.size(coord, axis = 0)
        #nodes = []

        #for i in range(nnode):
        #    nodes.append(v.Sphere(c='white').scale(1.5*scale).pos([coord[i,0],coord[i,1],coord[i,2]]).alpha(alpha))

        if render_nodes == True:
            nodes = tools.get_node_elements(coord,scale,alpha)

        meshes = []
        nel = np.size(edof, axis = 0)
        #print(coord)
        for i in range(nel):
            #print(coord[i])
            coords = tools.get_coord_from_edof(edof[i,:],dof,4)

            #print(coords)
            #print(coord[coords,:])

            mesh = v.Mesh([coord[coords,:],[[0,1,2,3],[4,5,6,7],[0,3,7,4],[1,2,6,5],[0,1,5,4],[2,3,7,6]]],alpha=alpha)
            meshes.append(mesh)

        #print(meshes[1])

        #if export is not None:
        #    v.io.write(mesh, export+".vtk")

        #plot_window.render_solid_geometry(meshes)
        if render_nodes == True:
            plot_window.render_geometry(meshes,nodes,merge=True)
        else:
            plot_window.render_geometry(meshes,merge=True)
        #plot_window.anim_element_type = 4

    #plot_window.anim_coord = coord

    elif element_type == 5:
        #scale = 0.02
        #nnode = np.size(coord, axis = 0)
        #nodes = []
        #ndof = np.size(a, axis = 0)
        ncoord = np.size(coord, axis = 0)
        #def_nodes = []
        #def_coord = np.zeros([ncoord,3])

        #for i in range(nnode):
        #    nodes.append(v.Sphere(c='white').scale(1.5*scale).pos([coord[i,0],coord[i,1],coord[i,2]]).alpha(alpha))
        if render_nodes == True:
            nodes = tools.get_node_elements(coord,scale,alpha)


        nel = np.size(edof, axis = 0)
        elements = []
        
        res = 4

        vmin, vmax = np.min(el_values), np.max(el_values)

        el_values_array = np.zeros((1,4*res))[0,:]

        for i in range(nel):
            coord1,coord2 = tools.get_coord_from_edof(edof[i,:],dof,5)
            #print(coord1,coord2)
            #tube = v.Cylinder([[coord[coord1,0],coord[coord1,1],coord[coord1,2]],[coord[coord2,0],coord[coord2,1],coord[coord2,2]]],r=scale,res=res,c='white').alpha(alpha)
            #elements.append(tube)

            if nseg > 2:
                steps = np.float32(1/(nseg-1))
                #print(steps)

                dx = (coord[coord2,0]-coord[coord1,0])*steps
                dy = (coord[coord2,1]-coord[coord1,1])*steps
                dz = (coord[coord2,2]-coord[coord1,2])*steps
                #print(coord2,0,def_coord[coord1,0])

                for j in range(nseg-1):
                    x1 = coord[coord1,0]+dx*j
                    y1 = coord[coord1,1]+dy*j
                    z1 = coord[coord1,2]+dz*j

                    x2 = coord[coord1,0]+dx*(j+1)
                    y2 = coord[coord1,1]+dy*(j+1)
                    z2 = coord[coord1,2]+dz*(j+1)

                    tube = v.Cylinder([[x1,y1,z1],[x2,y2,z2]],r=scale,res=res,c=color).alpha(alpha)
                    elements.append(tube)

                    if el_values is not None:
                        #el_value1 = el_values[nseg*i,:]
                        #el_value2 = el_values[nseg*i+1,:]
                        el_value1 = el_values[nseg*i]
                        el_value2 = el_values[nseg*i+1]

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

                        tube.cmap(colormap, el_values_array, on="points", vmin=vmin, vmax=vmax)
                        #tube.addScalarBar(label)


                #for k in range(nseg):
                #    el_value[k] = el_values[nseg*i,:]



            else:
                tube = v.Cylinder([[coord[coord1,0],coord[coord1,1],coord[coord1,2]],[coord[coord2,0],coord[coord2,1],coord[coord2,2]]],r=scale,res=res,c=color).alpha(alpha)
                elements.append(tube)

                if el_values is not None:
                    #el_value1 = el_values[2*i,:]
                    #el_value2 = el_values[2*i+1,:]
                    el_value1 = el_values[2*i]
                    el_value2 = el_values[2*i+1]

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

                    tube.cmap(colormap, el_values_array, on="points", vmin=vmin, vmax=vmax)
                    #tube.addScalarBar(label)


        #if export is not None:
        #    v.io.write(mesh, export+".vtk")

        #plot_window.render_beam_geometry(nodes,elements)
        if render_nodes == True:
            plot_window.render_geometry(elements,nodes)
        else:
            plot_window.render_geometry(elements)
        #plot_window.anim_element_type = 5
        #anim.elements = elements
        #anim.nodes = nodes

        return elements
    
        #anim.type = 1
        #anim.coords = coords

    elif element_type == 6:
        print("plate")

    else:
        print("Invalid element type, please declare 'element_type'. The element types are; 1: Spring, 2: Bar, 3: Flow, 4: Solid, 5: Beam, 6: Plate")
        sys.exit()
        #try:
        #    raise ValueError("Invalid element type, please declare 'element_type'. The element types are; 1: Spring, 2: Bar, 3: Flow, 4: Solid, 5: Beam, 6: Plate")
        #except ValueError:
            #return ValueError("Invalid element type, please declare 'element_type'. The element types are; 1: Spring, 2: Bar, 3: Flow, 4: Solid, 5: Beam, 6: Plate")
        #    return
            #print('test')

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

def add_scalar_bar(meshes,label,pos=[0.8,0.05]):
    app = init_app()
    plot_window = VedoPlotWindow.intance().plot_window

    #mesh = v.merge(meshes,flag=True)

    #scalar_bar = v.addons.addScalarBar(mesh,title=label,pos=pos)

    nel = np.size(meshes, axis = 0)
    for i in range(nel):
        scalar_bar = v.addons.addScalarBar(meshes[i],title=label,pos=pos)
    plot_window.render_addon(scalar_bar)

def add_legend(meshes):
    app = init_app()
    plot_window = VedoPlotWindow.intance().plot_window

    #mesh = v.merge(meshes,flag=True)

    legend_box = v.addons.LegendBox(meshes,pos='top-left')

    #nel = np.size(meshes, axis = 0)
    #for i in range(nel):
    #    scalar_bar = v.addons.LegendBox(meshes[i],pos='top-left')
    plot_window.render_addon(legend_box)

def add_text(text,color,pos=[0,1]):
    app = init_app()
    plot_window = VedoPlotWindow.intance().plot_window

    #mesh = v.merge(meshes,flag=True)

    #legend_box = v.addons.LegendBox(meshes,pos='top-left')
    #msg = v.Text2D(text, pos=pos, c=color, bg='r9', alpha=0.8)
    msg = v.Text2D(text, pos=pos, c=color)

    #nel = np.size(meshes, axis = 0)
    #for i in range(nel):
    #    scalar_bar = v.addons.LegendBox(meshes[i],pos='top-left')
    plot_window.render_addon(msg)

def add_silouette(text,color,pos=[0,1]):
    app = init_app()
    plot_window = VedoPlotWindow.intance().plot_window

    #mesh = v.merge(meshes,flag=True)

    #legend_box = v.addons.LegendBox(meshes,pos='top-left')
    #msg = v.Text2D(text, pos=pos, c=color, bg='r9', alpha=0.8)
    msg = v.Text2D(text, pos=pos, c=color)

    #nel = np.size(meshes, axis = 0)
    #for i in range(nel):
    #    scalar_bar = v.addons.LegendBox(meshes[i],pos='top-left')
    plot_window.render_addon(msg)

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###


def draw_displaced_geometry(
    edof,
    coord,
    dof,
    element_type,
    a,
    el_values=None,
    label=None,
    colormap='jet',
    scale=0.02,
    alpha=1,
    def_scale=1,
    nseg=2,
    render_nodes=True,
    color='white',
    offset_x = 0.0,
    offset_y = 0.0,
    offset_z = 0.0
    #export=None
    ):
    # Element types: 1: Spring, 2: Bar, 3: Flow, 4: Solid, 5: Beam, 6: Plate

    # 2 node: 1,2,5 (1-3D)
    # 3 node: 3,4 (Triangular 2D)
    # 4 node: 3,4,6 (Quadratic/rectangular/isoparametric 2D)
    # 8 node: 3,4 (Isoparametric 2D or 3D)

    app = init_app()
    plot_window = VedoPlotWindow.intance().plot_window
    #animations = animations().instance()

    if element_type == 1:
        print("spring")
        ndof = np.size(dof, axis = 0)*np.size(dof, axis = 1)

        #ndof = np.size(dof, axis = 0)*np.size(dof, axis = 1)
        #nnode = np.size(coord, axis = 0)
        #ncoord = np.size(coord, axis = 1) # Använd för att visualisera 1-3 dim.
        #nodes = []

        #for i in range(nnode):
        #    nodes.append(v.Sphere(c='white').scale(1.5*scale).pos([coord[i,0],0,0]).alpha(alpha))
        if render_nodes == True:
            nodes = tools.get_node_elements(coord,scale,alpha)

        nel = np.size(edof, axis = 0)
        elements = []

        for i in range(nel):
            coord1,coord2 = tools.get_coord_from_edof(edof[i,:],dof,element_type)

            #print(coord[coord1,0])
            #print(coord[coord2,0])
            #spring = v.Spring([coord[coord1,0],0,0],[coord[coord2,0],0,0])
            spring = v.Spring([coord[coord1,0]+offset_x,0+offset_y,0+offset_z],[coord[coord2,0]+offset_x,0+offset_y,0+offset_z],r=scale).alpha(alpha)
            elements.append(spring)

        #print(elements,nodes)
        if render_nodes == True:
            plot_window.render_geometry(elements,nodes)
        else:
            plot_window.render_geometry(elements)

    elif element_type == 2:
        #scale = 0.02
        #ndof = np.size(a, axis = 0)
        ndof = np.size(dof, axis = 0)*np.size(dof, axis = 1)
        ncoord = np.size(coord, axis = 0)
        def_nodes = []
        def_coord = np.zeros([ncoord,3])

        for i in range(0, ncoord):
            #if a.any() == None:
            #    x = coord[i,0]
            #    y = coord[i,1]
            #    z = coord[i,2]
            #else:
            a_dx, a_dy, a_dz = tools.get_a_from_coord(i,6,a,def_scale)

            x = coord[i,0]+a_dx
            y = coord[i,1]+a_dy
            z = coord[i,2]+a_dz

            def_coord[i] = [x,y,z]

            def_nodes.append(v.Sphere(c='white').scale(1.5*scale).pos([x,y,z]).alpha(alpha))

        #beamdata.def_coord = def_coord

        nel = np.size(edof, axis = 0)
        def_elements = []
        #scale = 0.02
        res = 4

        #def_coords = []

        vmin, vmax = np.min(el_values), np.max(el_values)

        #if nseg > 2:
        #    el_values = np.zeros((1,nseg))[0,:]

        el_values_array = np.zeros((1,4*res))[0,:]

        for i in range(nel):
            coord1,coord2 = tools.get_coord_from_edof(edof[i,:],dof,element_type)

            #tube = v.Cylinder(height=h/scale,res=res).scale(scale).pos([x,y,z]).orientation([dx,dy,dz])


    


            tube = v.Cylinder([[def_coord[coord1,0],def_coord[coord1,1],def_coord[coord1,2]],[def_coord[coord2,0],def_coord[coord2,1],def_coord[coord2,2]]],r=scale,res=res,c=color).alpha(alpha)
            def_elements.append(tube)
            if el_values is not None:
                #el_value1 = el_values[2*i,:]
                #el_value2 = el_values[2*i+1,:]

                el_values_array[1] = el_values[i]
                el_values_array[3] = el_values[i]
                el_values_array[5] = el_values[i]
                el_values_array[7] = el_values[i]
                el_values_array[12] = el_values[i]
                el_values_array[13] = el_values[i]
                el_values_array[14] = el_values[i]
                el_values_array[15] = el_values[i]

                el_values_array[0] = el_values[i]
                el_values_array[2] = el_values[i]
                el_values_array[4] = el_values[i]
                el_values_array[6] = el_values[i]
                el_values_array[8] = el_values[i]
                el_values_array[9] = el_values[i]
                el_values_array[10] = el_values[i]
                el_values_array[11] = el_values[i]

                tube.cmap(colormap, el_values_array, on="points", vmin=vmin, vmax=vmax)
                #tube.addScalarBar(label)

        #plot_window.render_beam_geometry(def_nodes,def_elements)

        if render_nodes == True:
            plot_window.render_geometry(def_elements,def_nodes)
        else:
            plot_window.render_geometry(def_elements)

        return def_elements

    elif element_type == 3:
        print("flow")

    elif element_type == 4:
        print("solid")
        #plot_window.anim_element_type = 2

    elif element_type == 5:
        #scale = 0.02
        #ndof = np.size(a, axis = 0)
        ndof = np.size(dof, axis = 0)*np.size(dof, axis = 1)
        ncoord = np.size(coord, axis = 0)
        def_nodes = []
        def_coord = np.zeros([ncoord,3])

        for i in range(0, ncoord):
            #if a.any() == None:
            #    x = coord[i,0]
            #    y = coord[i,1]
            #    z = coord[i,2]
            #else:
            a_dx, a_dy, a_dz = tools.get_a_from_coord(i,6,a,def_scale)

            x = coord[i,0]+a_dx
            y = coord[i,1]+a_dy
            z = coord[i,2]+a_dz

            def_coord[i] = [x,y,z]

            def_nodes.append(v.Sphere(c='white').scale(1.5*scale).pos([x,y,z]).alpha(alpha))

        #beamdata.def_coord = def_coord

        nel = np.size(edof, axis = 0)
        def_elements = []
        #scale = 0.02
        res = 4

        #def_coords = []

        vmin, vmax = np.min(el_values), np.max(el_values)

        #if nseg > 2:
        #    el_values = np.zeros((1,nseg))[0,:]

        el_values_array = np.zeros((1,4*res))[0,:]

        for i in range(nel):
            coord1,coord2 = tools.get_coord_from_edof(edof[i,:],dof,5)

            #tube = v.Cylinder(height=h/scale,res=res).scale(scale).pos([x,y,z]).orientation([dx,dy,dz])


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
                    tube = v.Cylinder([[x1,y1,z1],[x2,y2,z2]],r=scale,res=res,c=color).alpha(alpha)
                    def_elements.append(tube)
                    if el_values is not None:
                        #el_value1 = el_values[nseg*i,:]
                        #el_value2 = el_values[nseg*i+1,:]
                        el_value1 = el_values[nseg*i]
                        el_value2 = el_values[nseg*i+1]

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

                        tube.cmap(colormap, el_values_array, on="points", vmin=vmin, vmax=vmax)
                        #tube.addScalarBar(label)


                #for k in range(nseg):
                #    el_value[k] = el_values[nseg*i,:]



            else:
                tube = v.Cylinder([[def_coord[coord1,0],def_coord[coord1,1],def_coord[coord1,2]],[def_coord[coord2,0],def_coord[coord2,1],def_coord[coord2,2]]],r=scale,res=res,c=color).alpha(alpha)
                def_elements.append(tube)

                if el_values is not None:
                    #el_value1 = el_values[2*i,:]
                    #el_value2 = el_values[2*i+1,:]
                    el_value1 = el_values[2*i]
                    el_value2 = el_values[2*i+1]

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

                    tube.cmap(colormap, el_values_array, on="points", vmin=vmin, vmax=vmax)
                    #tube.addScalarBar(label)

        #plot_window.render_beam_geometry(def_nodes,def_elements)

        if render_nodes == True:
            plot_window.render_geometry(def_elements,def_nodes)
        else:
            plot_window.render_geometry(def_elements)

        return def_elements
        
        #plot_window.anim_element_type = 1
        #anim.type = 1
        #anim.def_elements = def_elements
        #anim.def_nodes = def_nodes
        #anim.def_coords = def_coords

    elif element_type == 6:
        print("plate")
        #plot_window.anim_element_type = 2
    else:
        print("Invalid element type, please declare 'element_type'. The element types are; 1: Spring, 2: Bar, 3: Flow, 4: Solid, 5: Beam, 6: Plate")
        sys.exit()
        #try:
        #    raise ValueError("Invalid element type, please declare 'element_type'. The element types are; 1: Spring, 2: Bar, 3: Flow, 4: Solid, 5: Beam, 6: Plate")
        #except ValueError:
            #return ValueError("Invalid element type, please declare 'element_type'. The element types are; 1: Spring, 2: Bar, 3: Flow, 4: Solid, 5: Beam, 6: Plate")
        #    return
            #print('test')
    #plot_window.anim_def_coord = coord


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
"""
class animations():
    
    #def __init__(self):
    #    #Qt.QMainWindow.__init__(self)
    #    #self.initialize()
    #    #self.type = None

    #    # Type 1
    #    self.elements = None
    #   self.nodes = None
    #    self.def_elements = None
    #    self.def_nodes = None

    #    self.mesh = None
    #    self.def_mesh = None
    

    def __init__(self):

        self.start = 0
        self.end = 1
        self.steps = 101

    def __call__(self):
        #print('test')
        #self.elements = None
        #self.nodes = None
        #self.def_elements = None
        #self.def_nodes = None

        #self.mesh = None
        #self.def_mesh = None

        self.type = None
        self.coords = None
        self.def_coords = None

    def edit(self):
        print('edit_parameters')


    #def animate(coord_start,coord_end,el_values,element_type):
    def animate(self):
        app = init_app()
        plot_window = VedoPlotWindow.intance().plot_window

        t = np.arange(self.start, self.end, self.steps)
        print(t)
        
        #if self.elements != None:
        #    element_type = 1
        #elif self.mesh != None:
        #    element_type = 2
        
        if element_type == 1:
            ncoord = np.size(self.coord, axis = 0)
            def_coord = np.zeros([ncoord,3])
"""

def animate(steps,coord,def_coord,element_type):
    app = init_app()
    plot_window = VedoPlotWindow.intance().plot_window

    t = np.arange(0, 1+1/(steps-1), 1/(steps-1))
    print(t)
    #print(start)
    #print(end)
    #print(steps)
    
    #if self.elements != None:
    #    element_type = 1
    #elif self.mesh != None:
    #    element_type = 2
    
    if element_type == 5:
        ncoord = np.size(coord, axis = 0)
        def_coord = np.zeros([ncoord,3])


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###


#Start Calfem-vedo visualization
def show_and_wait():
    #app = Qt.QApplication(sys.argv)
    #window = MainWindow()
    app = init_app()
    app.exec_()

