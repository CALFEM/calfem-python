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
#import webbrowser
from scipy.io import loadmat
import calfem.core as cfc
import vedo_utils as vdu
import vtk

# Examples using this module:
    # exv1: Spring
    # exv2: 3D Truss (beams & bars)
    # exv3: 3D Flow
    # exv4a: 3D Solid (using multiple renderers)
    # exv4b: 3D Solid w/ animation
    # exv5: 2D Plate visualized in 3D

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

# --- Tillfälliga anteckningar ---


    # 2D mode: self.plotter.parallelProjection(value=True, at=0) + annat 'mode' för kameran
    # Balkdiagram: vedo.pyplot.plot()
    # Klicka på element, se nedan
            
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

    # Siluetter + mått, se nedan

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

    # self.plotter.show(mode=?), 0 -> trackball, 1 -> actor, 10 -> terrain
    # self.plotter.show(axes=?), 4 -> koordinataxlar, 4 -> koordinataxlar i hörnet, 7 -> mått, 12 -> gradskiva






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
        #self.n = 0
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
        #self.plotter[self.n].addCallback('mouse click', self.click)
        self.silcont = [None]
        self.click_msg = v.Text2D("", pos="bottom-center", bg='auto', alpha=0.1, font='Calco',c='black')
        #self.plotter[self.n].add(self.click_msg)

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
        #if evt.isAssembly: # endast för testning ifall en assembly skapas
        #    print('assembly')
        #    self.click_msg.text(evt.actor.info)
        if evt.actor:
            #print(evt.actor.mapper())
            #if self.silcont != [] or self.silcont != [None]:
            # Silouette
            sil = evt.actor.silhouette().lineWidth(5).c('red5')
            self.plt[evt.title].remove(self.silcont.pop()).add(sil)
            self.silcont.append(sil)

            #self.silcont[0].c('black')

            # Color
            #evt.actor.c('red5')

            #sil = evt.actor.silhouette().lineWidth(5).c('red5')
            #self.plt[evt.title].remove(self.silcont.pop()).add(evt.actor)
            #self.silcont.append(evt.actor)
            #print('Title: ',evt.title,', Plotter: ',self.plt[evt.title])
            self.click_msg.text(evt.actor.name)
            
            

            #self.plt[evt.title].remove(self.silcont.pop()).add(sil)
            #evt.interactor.add(sil)
            #evt.interactor.pop()
            #.remove(self.silcont.pop())
            #.add(sil)
            #self.plotter.remove(self.silcont.pop()).add(sil)
            #self.silcont.append(sil)
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

        #print('keyframes',self.keyframes)
        #print('keyframe_dict',self.keyframe_dict)

        for plot in range(self.rendered, self.fig+1):
            print('anim',self.mode_anim)
            #print('2D',self.mode_2D)
            #print(f'Figure {i}: 2D mode is {self.mode_2D[i]}')
            if self.mode_2D[plot] is True:
                pp = True
                #mode = 12
                #mode = 6
                mode = 5
            else:
                pp = False
                mode = 0
            
            #type(dof) is list
            #if self.keyframes[self.fig] != []:
            if self.mode_anim[plot] == True:
                v.settings.immediateRendering = True
                print('animation',plot+1)
                #if self.bg[self.fig] != 'white':
                #    opts = dict(bg=self.bg[self.fig], axes=4, interactive=False, title=f'Figure {i+1} - CALFEM vedo visualization tool')
                #else:
                #    opts = dict(axes=4, interactive=False, title=f'Figure {i+1} - CALFEM vedo visualization tool')
                opts = dict(axes=4, interactive=0, title=f'Figure {plot+1} - CALFEM vedo visualization tool')
                keyframes = self.keyframes[self.fig]
                keyframe_dict = self.keyframe_dict[self.fig]
                print(keyframes)
                print(keyframe_dict)
                nmesh = len(keyframes)
                #print(len(keyframes))
                #print(len(keyframe_dict))
                #plt = v.Plotter(**opts).show()
                plt = v.Plotter(**opts)

                #plt.enableErase()
                
                #plt.parallelProjection(value=pp)

                #plt.addCallback('mouse click', self.click)
                #self.plt[f'Figure {i+1} - CALFEM vedo visualization tool'] = plt

                plt += v.Text2D('Press ESC to exit', pos='bottom-middle')

                #msg = v.Text2D(text, pos=pos, alpha=1, c=color)
                #plot_window.msg[plot_window.fig] += [msg]

                #print(self.keyframes)
                #print(keyframe_dict[0])
                #print(keyframes[0][0])
                #sys.exit()






                for j in range(len(self.msg[plot])):
                    if self.bg[self.fig] == 'black':
                        self.msg[plot][j].c('w')
                    plt.add(self.msg[plot][j])
                        #plt += self.msg[i][j]
                        #plt.add(self.click_msg)

                #if self.proj[i]:
                for j in range(len(self.proj[plot])):
                    if self.bg[self.fig] == 'black':
                        self.proj[plot][j].c('w')
                    plt.add(self.proj[plot][j])

                #if self.rulers[i]:
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
                    
                    #it = 0
                    #plt += keyframes[0][keyframe_dict[0][1]]
                    #plt.resetCamera()
                    for key, val in keyframe_dict[0].items():
                        #v.show(keyframes[0][val], **opts)
                        #print(key,val,keyframes[0][val])
                        #if key == 1:
                            #plt.clear(keyframes[0][keyframe_dict[0][val]])
                        #    time.sleep(self.dt/10000)
                            #plt.clear()
                        #    plt.pop()
                            #plt.remove()
                        #    plt += keyframes[0][val]
                        #    plt.show(resetcam=False)
                        #else:
                        #    plt += keyframes[0][val]
                        #    plt.show(resetcam=True)

                        time.sleep(self.dt/1000)
                        for i in range(nmesh):
                            plt.pop()
                        

                        for i in range(nmesh):
                            plt += keyframes[i][val]
                        #if key == 1:
                        #    plt.show(resetcam=True)
                        #else:
                        #    plt.show(resetcam=False)
                        plt.show(resetcam=False)

                        

                        #plt.show(keyframes[0][val]).close()
                        
                        #plt.add(keyframes[0][val],resetcam=True)
                        #plt.allowInteraction()
                        #plt.enableRenderer(0)
                        #plt.render(resetcam=True)
                        #time.sleep(self.dt/1000)
                        
                        #it += 1



                    #plt.show()
                    if plt.escaped: break  # if ESC is hit during the loop

                plt.close()
                v.settings.immediateRendering = False

                #sys.exit()
                #plt.close()
                '''
                it = 0
                for j in zip(keyframes):
                    print(j)
                    sys.exit()
                #for j in range(1,len(self.keyframes[i])):
                    if it ==0:
                        #self.plt[f'Figure {i+1} - CALFEM vedo visualization tool'] += j
                        #self.plt[f'Figure {i+1} - CALFEM vedo visualization tool'].render(resetcam=True)
                        plt.show(j)
                        #plt += j
                        #plt.render(resetcam=True)
                        #self.plt[f'Figure {i+1} - CALFEM vedo visualization tool'].show(j)
                    else:
                        #time.sleep(self.dt/1000)
                        #self.plt[f'Figure {i+1} - CALFEM vedo visualization tool'].show(j)
                        #plt.clear()
                        #plt += mesh
                        #plt.render(resetcam=True)
                        #plt.show(mesh).interactive().close()
                        #print(j)
                        #self.plt[f'Figure {i+1} - CALFEM vedo visualization tool'].clear(keyframes[it-1])
                        #self.plt[f'Figure {i+1} - CALFEM vedo visualization tool'] += j
                        plt.clear(keyframes[it-1])
                        #plt.clear()
                        plt += j
                        #plt.add(j)
                        #plt.show(j).close()
                        plt.render(resetcam=True)
                        #plt.interactive()
                        #v.show(axes=4,interactive=True)
                        #plt.show(j)
                        #self.plt[f'Figure {i+1} - CALFEM vedo visualization tool'].render(resetcam=True)
                        #v.interactive().close()
                    it += 1

                    plt.close()
                '''
                #plt.show(mesh).interactive().close()

                #plt.show(mesh, title='Animation - CALFEM vedo visualization tool')
            
            else:
                print('non-animation',plot+1)
                if self.bg[self.fig] != 'white':
                    opts = dict(bg=self.bg[self.fig], mode=mode, axes=4, interactive=False, new=True, title=f'Figure {plot+1} - CALFEM vedo visualization tool')
                    if self.bg[self.fig] == 'black':
                        self.click_msg.c('w')
                        #for j in range(len(self.nodes[i])):
                        #    if self.nodes[i][j].c('black'):
                        #        self.nodes[i][j].c('w')
                        for j in range(len(self.geometries[plot])):
                            for k in range(len(self.geometries[plot][j])):
                                #print(self.geometries[i][j])
                                self.geometries[plot][j][k].c('w')
                else:
                    opts = dict(mode=mode, axes=4, interactive=False, new=True, title=f'Figure {plot+1} - CALFEM vedo visualization tool')
                plt = v.show(self.geometries[plot], self.meshes[plot], self.nodes[plot], self.click_msg, **opts)
                plt.parallelProjection(value=pp)

            #plt.addGlobalAxes(11)#
            #plt.addShadows()
            #plt.addScaleIndicator(pos=(0.7, 0.05), s=0.02, length=2, lw=4, c='k1', alpha=1, units='', gap=0.05)
            if self.mode_hover[plot] == True or self.mode_2D[plot] == True:
                plt.addCallback('MouseMove', self.click)
            else:
                plt.addCallback('mouse click', self.click)
            
            #plt += self.click_msg#.addHoverLegend(useInfo=True,s=1.25,maxlength=96)
            #print('Figure text: ',self.msg[i])
            #print('Projections: ',self.proj[i])
            #self.plt.append(plt)
            self.plt[f'Figure {plot+1} - CALFEM vedo visualization tool'] = plt

            #if self.msg[i]:
            for j in range(len(self.msg[plot])):
                if self.bg[self.fig] == 'black':
                    self.msg[plot][j].c('w')
                plt.add(self.msg[plot][j])
                    #plt += self.msg[i][j]
                    #plt.add(self.click_msg)

            #if self.proj[i]:
            for j in range(len(self.proj[plot])):
                if self.bg[self.fig] == 'black':
                    self.proj[plot][j].c('w')
                plt.add(self.proj[plot][j])

            #if self.rulers[i]:
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
        '''
        def animate(self):
            for i in range(self.rendered, self.fig+1):
                opts = dict(axes=4, interactive=False, new=True, title=f'Figure {i+1} - CALFEM vedo visualization tool')
                plt = v.show(self.keyframes[i] **opts)
        '''
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
# Plotting functions, for adding things to a rendering window
    
# Add scalar bar to a renderer
def add_scalar_bar(
    label,
    pos=[0.8,0.05],
    text_pos="bottom-right",
    font_size=24,
    color='black',
    on = 'mesh'
    #sx=1000,
    #sy=1000
    #size=(100, 50)
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

# Add text to a render window
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

    #while len(plot_window.msg) < plot_window.fig + 1:
        

    #plot_window.msg.append([])

    msg = v.Text2D(text, pos=pos, alpha=1, c=color)
    plot_window.msg[plot_window.fig] += [msg]

# Add 3D text to a scene
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

# Add axes to a mesh
def add_mesh_axes(mesh):
    """
    Routine for automatically adding axes for a mesh to a figure.

    :param list mesh: mesh for adding axes
    """
    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window

    axes = mesh.buildAxes()
    plot_window.dia_axes[plot_window.fig].append(axes)

# Add grid & axes
#def add_axes(obj=None, xtitle=None, ytitle=None, ztitle=None, xrange=None, yrange=None, zrange=None, c=None, numberOfDivisions=None, digits=None, limitRatio=0.04, axesLineWidth=1, gridLineWidth=1, htitle='', hTitleSize=0.03, hTitleFont=None, hTitleItalic=True, hTitleColor=None, hTitleJustify='bottom-center', hTitleRotation=0, hTitleOffset=(0, 0.01, 0), titleDepth=0, titleFont='', textScale=1.0, xTitlePosition=0.95, yTitlePosition=0.95, zTitlePosition=0.95, xTitleOffset=0.025, yTitleOffset=0.0275, zTitleOffset=0.02, xTitleJustify=None, yTitleJustify=None, zTitleJustify=None, xTitleRotation=0, yTitleRotation=0, zTitleRotation=0, xTitleBox=False, yTitleBox=False, xTitleSize=0.025, yTitleSize=0.025, zTitleSize=0.025, xTitleColor=None, yTitleColor=None, zTitleColor=None, xTitleBackfaceColor=None, yTitleBackfaceColor=None, zTitleBackfaceColor=None, xTitleItalic=0, yTitleItalic=0, zTitleItalic=0, xyGrid=True, yzGrid=False, zxGrid=False, xyGrid2=False, yzGrid2=False, zxGrid2=False, xyGridTransparent=False, yzGridTransparent=False, zxGridTransparent=False, xyGrid2Transparent=False, yzGrid2Transparent=False, zxGrid2Transparent=False, xyPlaneColor=None, yzPlaneColor=None, zxPlaneColor=None, xyGridColor=None, yzGridColor=None, zxGridColor=None, xyAlpha=0.075, yzAlpha=0.075, zxAlpha=0.075, xyFrameLine=None, yzFrameLine=None, zxFrameLine=None, xyFrameColor=None, yzFrameColor=None, zxFrameColor=None, xLineColor=None, yLineColor=None, zLineColor=None, xHighlightZero=False, yHighlightZero=False, zHighlightZero=False, xHighlightZeroColor='r', yHighlightZeroColor='g', zHighlightZeroColor='b', showTicks=True, xTickLength=0.015, yTickLength=0.015, zTickLength=0.015, xTickThickness=0.0025, yTickThickness=0.0025, zTickThickness=0.0025, xMinorTicks=1, yMinorTicks=1, zMinorTicks=1, tipSize=None, labelFont='', xLabelColor=None, yLabelColor=None, zLabelColor=None, xLabelSize=0.016, yLabelSize=0.016, zLabelSize=0.016, xLabelOffset=0.8, yLabelOffset=0.8, zLabelOffset=0.8, xLabelJustify=None, yLabelJustify=None, zLabelJustify=None, xLabelRotation=0, yLabelRotation=0, zLabelRotation=0, xAxisRotation=0, yAxisRotation=0, zAxisRotation=0, xValuesAndLabels=None, yValuesAndLabels=None, zValuesAndLabels=None, xyShift=0, yzShift=0, zxShift=0, xShiftAlongY=0, xShiftAlongZ=0, yShiftAlongX=0, yShiftAlongZ=0, zShiftAlongX=0, zShiftAlongY=0, xUseBounds=True, yUseBounds=False, zUseBounds=False, xInverted=False, yInverted=False, zInverted=False, useGlobal=False, tol=0.0001):
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
            #hTitleFont='Kanopus',
            #hTitleJustify='bottom-right',
            #hTitleColor='red2',
            #hTitleSize=0.035,
            #hTitleOffset=(0,0.075,0),
            #hTitleRotation=45,
            #zHighlightZero=True,
            #xyFrameLine=2, yzFrameLine=1, zxFrameLine=1,
            #xyFrameColor='red3',
            #xyShift=1.05, # move xy 5% above the top of z-range
            xyGrid=xyGrid,
            yzGrid=yzGrid,
            zxGrid=zxGrid,
            xyGrid2=xyGrid,
            yzGrid2=yzGrid,
            zxGrid2=zxGrid,
            xyGridTransparent=xyGridTransparent, yzGridTransparent=yzGridTransparent, zxGridTransparent=zxGridTransparent, xyGrid2Transparent=xyGrid2Transparent, yzGrid2Transparent=yzGrid2Transparent, zxGrid2Transparent=zxGrid2Transparent,
            #zxShift=1.0,
            #xTitleJustify='bottom-right',
            #xTitleOffset=-1.175,
            #xLabelOffset=-1.75,
            #yLabelRotation=90,
            #zInverted=True,
            #tipSize=0.25,
            xrange=xrange,
            yrange=yrange,
            zrange=zrange
          )

    #axes = v.Axes(obj=obj, xtitle=xtitle, ytitle=ytitle, ztitle=ztitle, xrange=xrange, yrange=yrange, zrange=zrange, c=c, numberOfDivisions=numberOfDivisions, digits=digits, limitRatio=limitRatio, axesLineWidth=axesLineWidth, gridLineWidth=gridLineWidth, htitle=htitle, hTitleSize=hTitleSize, hTitleFont=hTitleFont, hTitleItalic=hTitleItalic, hTitleColor=hTitleColor, hTitleJustify=hTitleJustify, hTitleRotation=hTitleRotation, hTitleOffset=hTitleOffset, titleDepth=titleDepth, titleFont=titleFont, textScale=textScale, xTitlePosition=xTitlePosition, yTitlePosition=yTitlePosition, zTitlePosition=zTitlePosition, xTitleOffset=xTitleOffset, yTitleOffset=yTitleOffset, zTitleOffset=zTitleOffset, xTitleJustify=xTitleJustify, yTitleJustify=yTitleJustify, zTitleJustify=zTitleJustify, xTitleRotation=xTitleRotation, yTitleRotation=yTitleRotation, zTitleRotation=zTitleRotation, xTitleBox=xTitleBox, yTitleBox=yTitleBox, xTitleSize=xTitleSize, yTitleSize=yTitleSize, zTitleSize=zTitleSize, xTitleColor=xTitleColor, yTitleColor=yTitleColor, zTitleColor=zTitleColor, xTitleBackfaceColor=xTitleBackfaceColor, yTitleBackfaceColor=yTitleBackfaceColor, zTitleBackfaceColor=zTitleBackfaceColor, xTitleItalic=xTitleItalic, yTitleItalic=yTitleItalic, zTitleItalic=zTitleItalic, xyGrid=xyGrid, yzGrid=yzGrid, zxGrid=zxGrid, xyGrid2=xyGrid2, yzGrid2=yzGrid2, zxGrid2=zxGrid2, xyGridTransparent=xyGridTransparent, yzGridTransparent=yzGridTransparent, zxGridTransparent=zxGridTransparent, xyGrid2Transparent=xyGrid2Transparent, yzGrid2Transparent=yzGrid2Transparent, zxGrid2Transparent=zxGrid2Transparent, xyPlaneColor=xyPlaneColor, yzPlaneColor=yzPlaneColor, zxPlaneColor=zxPlaneColor, xyGridColor=xyGridColor, yzGridColor=yzGridColor, zxGridColor=zxGridColor, xyAlpha=xyAlpha, yzAlpha=yzAlpha, zxAlpha=zxAlpha, xyFrameLine=xyFrameLine, yzFrameLine=yzFrameLine, zxFrameLine=zxFrameLine, xyFrameColor=xyFrameColor, yzFrameColor=yzFrameColor, zxFrameColor=zxFrameColor, xLineColor=xLineColor, yLineColor=yLineColor, zLineColor=zLineColor, xHighlightZero=xHighlightZero, yHighlightZero=yHighlightZero, zHighlightZero=zHighlightZero, xHighlightZeroColor=xHighlightZeroColor, yHighlightZeroColor=yHighlightZeroColor, zHighlightZeroColor=zHighlightZeroColor, showTicks=showTicks, xTickLength=xTickLength, yTickLength=yTickLength, zTickLength=zTickLength, xTickThickness=xTickThickness, yTickThickness=yTickThickness, zTickThickness=zTickThickness, xMinorTicks=xMinorTicks, yMinorTicks=yMinorTicks, zMinorTicks=zMinorTicks, tipSize=tipSize, labelFont=labelFont, xLabelColor=xLabelColor, yLabelColor=yLabelColor, zLabelColor=zLabelColor, xLabelSize=xLabelSize, yLabelSize=yLabelSize, zLabelSize=zLabelSize, xLabelOffset=xLabelOffset, yLabelOffset=yLabelOffset, zLabelOffset=zLabelOffset, xLabelJustify=xLabelJustify, yLabelJustify=yLabelJustify, zLabelJustify=zLabelJustify, xLabelRotation=xLabelRotation, yLabelRotation=yLabelRotation, zLabelRotation=zLabelRotation, xAxisRotation=xAxisRotation, yAxisRotation=yAxisRotation, zAxisRotation=zAxisRotation, xValuesAndLabels=xValuesAndLabels, yValuesAndLabels=yValuesAndLabels, zValuesAndLabels=zValuesAndLabels, xyShift=xyShift, yzShift=yzShift, zxShift=zxShift, xShiftAlongY=xShiftAlongY, xShiftAlongZ=xShiftAlongZ, yShiftAlongX=yShiftAlongX, yShiftAlongZ=yShiftAlongZ, zShiftAlongX=zShiftAlongX, zShiftAlongY=zShiftAlongY, xUseBounds=xUseBounds, yUseBounds=yUseBounds, zUseBounds=zUseBounds, xInverted=xInverted, yInverted=yInverted, zInverted=zInverted, useGlobal=useGlobal, tol=tol)
    plot_window.dia_axes[plot_window.fig].append(axes)

# Add silhouette with or without measurements to a renderer
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

    

    #print("Size of projections: ",len(plot_window.proj))
    #while len(plot_window.proj) < plot_window.fig + 1:
        
    #if rulers == True:
        #while len(plot_window.rulers) < plot_window.fig + 1:
            
    #print("Size of projections after loop: ",len(plot_window.proj))

    if plane == 'xy':
        #assem = v.merge(plot_window.meshes[plot_window.fig])
        assem = plot_window.meshes[plot_window.fig]
        plot_window.meshes[plot_window.fig]
        proj = assem.projectOnPlane('z').z(offset).silhouette('2d').c(color)
        plot_window.proj[plot_window.fig] += [proj]
        if rulers == True:
            ruler = v.addons.RulerAxes(proj, xtitle='', ytitle='', ztitle='', xlabel='', ylabel='', zlabel='', xpad=1, ypad=1, zpad=1, font='Normografo', s=None, italic=0, units='m', c=color, alpha=1, lw=1, precision=3, labelRotation=0, axisRotation=0, xycross=True)
            plot_window.rulers[plot_window.fig] += [ruler]
        """
        for i in range(len(plot_window.meshes[plot_window.fig])):
            proj = plot_window.meshes[plot_window.fig][i].clone().projectOnPlane('z').c(color).z(offset)
            plot_window.proj[plot_window.fig] += [proj]
        """
        #plot_window.meshes[plot_window.fig][0].addShadow(plane='x')
    elif plane == 'xz':
        meshes = []
        for i in range(len(plot_window.meshes[plot_window.fig])):
            #print(plot_window.meshes[plot_window.fig,i])
            meshes.append(plot_window.meshes[plot_window.fig][i].clone())
        #assem = v.merge(plot_window.meshes[plot_window.fig])
        assem = v.merge(meshes)
        #assem = v.Assembly(plot_window.meshes[plot_window.fig])
        #assem = plot_window.meshes[plot_window.fig][0]
        proj = assem.projectOnPlane('y').y(offset).silhouette('2d').c(color)
        plot_window.proj[plot_window.fig] += [proj]
        if rulers == True:
            ruler = v.addons.RulerAxes(proj, xtitle='', ytitle='', ztitle='', xlabel='', ylabel='', zlabel='', xpad=1, ypad=1, zpad=1, font='Normografo', s=None, italic=0, units='m', c=color, alpha=1, lw=1, precision=3, labelRotation=0, axisRotation=0, xycross=True)
            plot_window.rulers[plot_window.fig] += [ruler]
    elif plane == 'yz':
        #assem = v.merge(plot_window.meshes[plot_window.fig])
        assem = v.Assembly(plot_window.meshes[plot_window.fig])
        proj = assem.projectOnPlane('x').x(offset).silhouette('2d').c(color)
        plot_window.proj[plot_window.fig] += [proj]
        if rulers == True:
            ruler = v.addons.RulerAxes(proj, xtitle='', ytitle='', ztitle='', xlabel='', ylabel='', zlabel='', xpad=0, ypad=0.1, zpad=0.1, font='Normografo', s=None, italic=0, units='m', c=color, alpha=1, lw=1, precision=3, labelRotation=0, axisRotation=0, xycross=True)
            plot_window.rulers[plot_window.fig] += [ruler]
    else:
        print("Please choose a plane to project to. Set plane to 'xy', 'xz' or 'yz'")
        sys.exit()

    

    #plot_window.proj.append([])

    #plot_window.proj[plot_window.fig] += [proj]

#Add measurements to a renderer
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

    #while len(plot_window.rulers) < plot_window.fig + 1:
        #plot_window.rulers.append([])

    assem = v.merge(plot_window.meshes[plot_window.fig])

    ruler = v.addons.RulerAxes(assem, xtitle=xtitle, ytitle=ytitle, ztitle=ztitle, xlabel=xlabel, ylabel=ylabel, zlabel=zlabel, xpad=0.1, ypad=0.1, zpad=0.1, font='Normografo', s=None, italic=0, units='m', c=(0,0,0), alpha=alpha, lw=1, precision=3, labelRotation=0, axisRotation=0, xycross=False)

    plot_window.rulers[plot_window.fig] += [ruler]

# Beam diagrams
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

    #print('ex',ex)
    #print('ey',ey)
    #print('ez',ez)
    print('es',es)
    #print('eci',eci)

    nel = np.size(ex,0)
    nseg = np.size(eci,1)
    upd_scale = scale*(1/np.max(np.absolute(es)))

    for i in range(nel):
        #l = np.sqrt( (ex[i,1]-ex[i,0])**2 + (ey[i,1]-ey[i,0])**2 + (ez[i,1]-ez[i,0])**2 )
        #print('l',l)
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
        #plot_window.dia_points[plot_window.fig].append(graph)
        for j in range(len(pts)):
            plot_window.dia_points[plot_window.fig].append(pts[j])
        
        if invert == True:
            ticks = -np.linspace(np.max(es[i])*upd_scale, np.min(es[i])*upd_scale, 10)
            #ticks = np.flip(ticks)
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






        #axes.yrange(np.min(es[i]),np.max(es[i]))
        #print(axes.unpack('yValuesAndLabels'))
        #axes.unpack('yAxis').scale(.00005)
        plot_window.dia_axes[plot_window.fig].append(axes)


        #text = []

        #for key,val in points.items():
        #    pts.append(v.Point(pos=val[0], r=12, c='black', alpha=1))
        #    text.append(v.Text3D(key, pos=val[0], s=scale, font='Normografo', hspacing=1.15, vspacing=2.15, depth=0, italic=False, justify='bottom-left', c='black', alpha=1, literal=False))

            
            #plot_window.meshes[plot_window.fig].append(text)
    
    
    #print('x',x)

# Principal stresses
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
    #print(n)

    upd_scale = (1/np.max(val))*scale
    #upd_scale = (1/np.max(vec))*scale

    if ed is None:
        ed = np.zeros((nel,3))

    #x_comp = eigenvectors[0, :]
    #y_comp = eigenvectors[1, :]
    #labels = ['$v_1$', '$v_2$']
    #plot(x_comp, y_comp, ['r','b'], [-1,1], [-1,1], '$x_1$', '$x_2$', 'Plot of eigenvectors', labels, offsets)

    #print('Eigenvalues el. 1')
    #print(val[0,:])
    #print('Eigenvectors el. 1')
    #print(vec[0,:,:])
    #print('---')
    pts = []
    points = []
    vectors = []
    text = []
    vmin = np.min(val)
    vmax = np.max(val)
    #vmin = np.min(vec)
    #vmax = np.max(vec)
    values = []
    for i in range(nel):
        x = np.average(ex[i])+ed[i,0]
        y = np.average(ey[i])+ed[i,1]
        z = np.average(ez[i])+ed[i,2]
        for j in range(n):
            points.append([x,y,z])
            text.append(f'Principal stress {j+1} at El. {i+1}: {np.round(val[i,j],3)} {unit}')
            #value = val[i,j]
            values.append([val[i,j],val[i,j],val[i,j],val[i,j],val[i,j],val[i,j]])
            vector = vec[i, :, j]*val[i,j]*upd_scale
            vectors.append(vector)

        #point = v.Point([x,y,z])
        #pts.append(point)

        

        #x_comp = vec[i, 0, :]
        #y_comp = vec[i, 1, :]
        #z_comp = vec[i, 2, :]

        #print('components')
        #print(x_comp)
        #print(y_comp)
        #print(z_comp)
        #print('---')
        #sys.exit()
        '''
        for j in range(n):
            x_comp = vec[i, j, :]
            x_comp = vec[i, :]
        '''

    #pointcloud = v.merge(pts)
    #plot_window.meshes[plot_window.fig].append(pointcloud)

    #quiver = v.pyplot.quiver(points, vectors, c='k', alpha=1, shaftLength=0.8, shaftWidth=0.05, headLength=0.25, headWidth=0.2, fill=True)

    plot = vdu.vectors(points, vectors, c='k', alpha=1, shaftLength=0.8, shaftWidth=0.05, headLength=0.25, headWidth=0.2, fill=True, text=text, vmax=vmax, vmin=vmin, cmap=colormap, values=values)

    plot_window.vectors[plot_window.fig].extend(plot)

# Flux (or other) vectors
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
        #text.append('')
        text.append(f'Flux at El. {i+1}: {np.round(flux_tot,3)} {unit}')
        #value = val[i,j]
        values.append([flux_tot,flux_tot,flux_tot,flux_tot,flux_tot,flux_tot])
        #vector = vec[i, :, j]*val[i,j]*upd_scale
        vectors.append(vec[i]*upd_scale)

    print('Points & vectors')
    print(points[0])
    print(vectors[0])
    print('---')

    plot = vdu.vectors(points, vectors, c='k', alpha=1, shaftLength=0.8, shaftWidth=0.05, headLength=0.25, headWidth=0.2, fill=True, text=text, vmax=vmax, vmin=vmin, cmap=colormap, values=values)

    plot_window.vectors[plot_window.fig].extend(plot)




"""
#def add_vectors(edof,coord,dof,v,element_type):
def add_vectors(ex,ey,ez):
    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window

    nel = np.size(ex,0)

    #v.Tensors([ex,ey,ez],source='Arrow')

    #for i in range(nel)


    '''
    nel = np.size(edof, axis = 0)
    x = np.zeros((nel,1))
    y = np.zeros((nel,1))
    z = np.zeros((nel,1))
    for i in range(nel):
        coords = vdu.get_coord_from_edof(edof[i],dof,element_type)
        x[i] = np.average(coord[coords,0])
        y[i] = np.average(coord[coords,1])
        z[i] = np.average(coord[coords,2])

        mag = np.sqrt(v[i][0]**2 + v[i][1]**2 + v[i][2]**2)
        alpha = np.arccos(v[i][0]/mag)
        beta = np.arccos(v[i][1]/mag)
        gamma = np.arccos(v[i][2]/mag)
        #print('mag',mag)
        #print('alpha',alpha)
        #print('beta',beta)
        #print('gamma',gamma)
        #p1 = [-np.cos(alpha)]
        #for j in zip(coords):
        #    coord[j]
    #arrow = v.Arrow().scale(0.04)
    #field = v.Glyph([x,y,z],arrow,orientationArray=vectors)
    #field = v.Tensors([x,y,z],source='Arrow')
    plot_window.vectors[plot_window.fig] += [field]
    '''

def tensors(ex,ey,ez,array, ed = None,scale=1):
    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window

    nel = np.size(ex,0)

    x = np.zeros((nel,1))
    y = np.zeros((nel,1))
    z = np.zeros((nel,1))

    pts = []
    for i in range(nel):
        x = np.average(ex[i])
        y = np.average(ey[i])
        z = np.average(ez[i])
        point = v.Point([x,y,z])
        #point.SetInputData(tensors[:,:,i])
        #print(point.inputdata())
        #point.inputdata().Tensors(tensor[i])
        pts.append(point)

    pointcloud = v.merge(pts)

    upd_scale = (1/np.max(array))*scale

    arrow = v.Cone().scale(upd_scale)

    #print(array[:,0]*0.00000000000000001)

    

    gl = v.Glyph(pointcloud,arrow,np.transpose(array),scaleByVectorSize=True, colorByVectorSize=True)

    #ag = vtk.vtkRandomAttributeGenerator()
    #ag.SetInputData(pointcloud.polydata())
    #ag.GenerateAllDataOn()
    #ag.Update()

    #print(ag.GetOutput())

    #ts = v.Tensors(ag.GetOutput(), source='cube', scale=0.1)

    #v.show(pointcloud, ts, interactive=True).close()

    #plot_window.meshes[plot_window.fig] += [pointcloud]
    #plot_window.meshes[plot_window.fig] += [gl]
    plot_window.meshes[plot_window.fig].append(gl)


    #pointcloud.SetInputData(tensors)

    #pts = v.pointcloud.Points([x,y,z])

    #v.Tensors(pointcloud)
"""





# Test of isolines/contour plots
'''
def eliso(mesh):
    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window

    isol = mesh.isolines(n=10).color('b')

    plot_window.meshes[plot_window.fig].append(isol)


def elcont(mesh):
    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window

    isob = mesh.isobands(n=5)

    plot_window.meshes[plot_window.fig].append(isob)
'''




### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
# Main functions for plotting geometries, meshes & results from CALFEM calculations

# NOTE: Only rectangular surfaces supported
def draw_geometry(points=None,lines=None,surfaces=None,scale=0.05,points_alpha=1,lines_alpha=1,surfaces_alpha=1):
    """
    Routine for geometry for spring, bar, flow, solid, beam or plate elements.

    :param array edof: element topology by degrees of freedom [nel x (n_dofs_per_node)|(n_dofs_per_node+1)*n_nodes ]

    :return array el_values: Global scalar values for element
    """
    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window

    if surfaces == None and lines == None and points == None:
        print("Please input either: (points), (points, lines) or (points, lines, surfaces) from geometry module")
        sys.exit()
    else:
        if surfaces is not None:
            print('-----------------')
            print('points',points)
            print('-----------------')
            print('lines',lines)
            print('-----------------')
            print('surfaces',surfaces)
            print('-----------------')


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
            #s_text = []
            for key,val in surfaces.items():
                surface_lines = val[1]
                

                point_key_dict = {}
                print(' ')
                npts = len(val[1])
                print('Surface lines:',surface_lines)
                #print(f'Surface line {lns[i]}, points:',lines[i][1])
                #sys.exit()
                #print('-----------------')
                #print('Surface lines:',val[1])
                indx = 0
                for line in surface_lines:
                    #point_keys.append(lines[line][1])
                    point_key_dict[indx] = lines[line][1]
                    indx += 1

                point_list = []
                for i in range(npts):
                    l0 = point_key_dict[i]
                    if i == npts - 1:
                        l1 = point_key_dict[0]
                    else:
                        l1 = point_key_dict[i+1]

                    print(f'First point for line {i}',l0)
                    print(f'Second point for line {i}',l1)

                    if l0[1] == l1[0] or l0[1] == l1[1]:
                        print('Ok, adding...')
                        point_list.append(l0[0])
                    elif l0[0] == l1[0] or l0[0] == l1[1]:
                        print('Flipping, adding...')
                        point_list.append(l0[1])
                    else:
                        print('ERROR')
                        #point_list.append(l0[1])

                    #print('-----------------')

                '''
                for line in range(len(point_keys)-1):
                    print(point_keys[line])
                    first_point = point_keys[line]
                    second_point = point_keys[line+1]
                    if first_point[]
                '''


                #sys.exit()

                print('Point keys:',point_list)
                print('-----------------')
                print(' ')
                #point_keys = list(dict.fromkeys(point_keys))
                
                
                coords = []
                for i in point_list:
                    #print(f'Point {i}, coordinates:',points[i][0])
                    coords.append(points[i][0])



                #print('Point keys, dupl. removed:',point_keys)
                print('Surface coordinates:',coords)

                print(f'Surface {key} lines:',val[1])

                pln = v.Mesh([coords, [[0,1,2,3]]], alpha=surfaces_alpha, c='grey6')
                
                #if key == 0 or key == 1 or key == 5:
                    
                #    pln = v.Mesh([coords, [[0,1,2,3]]], alpha=surfaces_alpha, c='grey6')
                #elif key == 2 or key == 3 or key == 4:
                #    pln = v.Mesh([coords, [[0,1,3,2]]], alpha=surfaces_alpha, c='grey6')
                surf.append(pln)
                pln.name = f"Surface {key}"

                print('-----------------')

                """
                    if indx == 0:
                        point_keys.append(lines[i][1][0])
                        '''
                    elif indx == npts - 1:
                        #print('Last point should trigger')
                        if point_keys[indx-1] != lines[i][1][0]:
                            point_keys.append(lines[i][1][1])
                        else:
                            point_keys.append(lines[i][1][0])
                        '''
                    else:
                        print('Index:',indx)
                        #print(val[1])
                        #print(point_keys)
                        #print(point_keys[indx])
                        #print(lines[i][1][0])
                        #if point_keys[indx-1] != lines[i][1][0]:
                        if lines[i][1][0] in point_keys:
                            if lines[i][1][1] in point_keys:
                                print('ERROR')
                                print(lines[i][1][0])
                                print(lines[i][1][1])
                            else:
                                point_keys.append(lines[i][1][1])
                        else:
                            point_keys.append(lines[i][1][0])

                    #point_keys.append(lines[i][1][1])
                    indx += 1
                """
                #print('Point keys:',point_keys)
                #point_keys = list(dict.fromkeys(point_keys))
                
                '''
                coords = []
                for i in point_keys:
                    #print(f'Point {i}, coordinates:',points[i][0])
                    coords.append(points[i][0])



                print('Point keys, dupl. removed:',point_keys)
                print('Surface coordinates:',coords)
                print('-----------------')


                pln = v.Mesh([coords, [[0,1,2,3]]], alpha=surfaces_alpha, c='grey6')
                surf.append(pln)
                pln.name = f"Surface {key}"
                '''
                '''
                ### NOTE: only 4 point surfaces implemented
                l12 = lines[val[1][0]][1]
                l34 = lines[val[1][2]][1]

                p1=points[l12[0]][0]
                p2=points[l12[1]][0]
                p3=points[l34[0]][0]
                p4=points[l34[1]][0]

                x = np.average([p1[0],p2[0],p3[0],p4[0]])
                y = np.average([p1[1],p2[1],p3[1],p4[1]])
                z = np.average([p1[2],p2[2],p3[2],p4[2]])

                #print(p1)
                normal = np.cross(np.asarray(p2) - np.asarray(p1), np.asarray(p3) - np.asarray(p2))

                dx = np.asarray(p1) + 0.5 * (np.asarray(p2) - np.asarray(p1))
                dy = np.asarray(p1) + 0.5 * (np.asarray(p3) - np.asarray(p1))



                print('normal',normal)

                print('dx',dx)
                print('dy',dy)

                print('------------')

                if normal[0] != 0.0:
                    sy = dx[1]*2
                    sx = dy[2]*2
                elif normal[1] != 0.0:
                    sy = dx[0]*2
                    sx = dy[2]*2
                elif normal[2] != 0.0:
                    sx = dx[0]*2
                    sy = dy[1]*2


                #print(x)
                #print(y)
                #print(z)
                #print(p4)
                '''
                #sys.exit()
                
                
                #pln = v.Plane(pos=(x, y, z), normal=normal, sx=sx, sy=sy, c='gray6', alpha=1)
                #surf.append(pln)
                #pln.name = f"Surface {key}"
                #s_text.append(v.Text3D(key, pos=[x, y, z], s=scale, font='Normografo', hspacing=1.15, vspacing=2.15, depth=0, italic=False, justify='bottom-left', c='black', alpha=1, literal=False))

            plot_window.geometries[plot_window.fig].append(surf)
            #plot_window.geometries[plot_window.fig].append(s_text)

            
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
            print("Please provide point coordinates along with lines")
            sys.exit()

        

    #print(points)

    #lines = []

    #for i in range(len(points)):
        #print(i)
    #    for j in range(len(points)):
    #        lines.append([i,j])

    
    #geometry
    #pts = v.Points(points).ps(10)#.renderPointsAsSpheres()
    #geometry = v.utils.geometry(pts.tomesh())
    #dly = v.delaunay2D(pts, mode='fit').lw(1)

    #if vol == True:
        #volume = v.mesh2Volume(geometry).mode(4)
        #plot_window.meshes[plot_window.fig].append(volume)
    #    plot_window.meshes[plot_window.fig].append(geometry)
    #else:
    
'''
Element types: 1: Spring, 2: Bar, 3: Flow, 4: Solid, 5: Beam, 6: Plate

    2 node: 1,2,5 (1-3D)
    3 node: 3,4 (Triangular 2D)
    4 node: 3,4,6 (Quadratic/rectangular/isoparametric 2D)
    8 node: 3,4 (Isoparametric 2D or 3D)
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

    :param array edof: element topology by degrees of freedom [nel x (n_dofs_per_node)|(n_dofs_per_node+1)*n_nodes ]

    :return array el_values: Global scalar values for element
    """

    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window

    #eq_els = np.transpose(eq_els)
    #eq = np.transpose(eq)

    if np.size(coord, axis = 1) == 1:
        print('only x')
        coord = np.append(coord, np.zeros((np.size(coord, axis = 0),1)), axis=1)
        coord = np.append(coord, np.zeros((np.size(coord, axis = 0),1)), axis=1)
        #np.c_[ coord, np.zeros((np.size(coord, axis = 0),1)), np.zeros((np.size(coord, axis = 0),1)) ]
    elif np.size(coord, axis = 1) == 2:
        print('only x & y')
        coord = np.append(coord, np.zeros((np.size(coord, axis = 0),1)), axis=1)
        #np.c_[ coord, np.zeros((np.size(coord, axis = 0),1)) ]

    print('coord',coord)


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
                        #print(eq_dict)
                        print('j',j)
                        #print(eq)
                        #print(indx)
                        print(eq_els)
                        print(eq[indx])
                        eq_dict[j[0]] = eq[indx][0]
                        indx += 1

        for i in range(nel):
            coord1,coord2 = vdu.get_coord_from_edof(edof[i,:],dof,element_type)
            '''
            if element_type == 1 and spring == True:
                element = v.Spring([def_coord[coord1,0],def_coord[coord1,1],def_coord[coord1,2]],[def_coord[coord2,0],def_coord[coord2,1],def_coord[coord2,2]],r=scale*1.5,c=color).alpha(alpha)
                element.name = f"Spring element {i+1}"
                elements.append(element)
            '''
            if element_type == 1 and spring == True:
                #spring = v.Spring([coord[coord1,0],0,0],[coord[coord2,0],0,0],r=1.5*scale).alpha(alpha)
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
                #if element_type == 3 or element_type == 4:
                if np.any(np.isin(eq_els, i, assume_unique=True)) == True:
                    beam = v.Cylinder([[coord[coord1,0],coord[coord1,1],coord[coord1,2]],[coord[coord2,0],coord[coord2,1],coord[coord2,2]]],r=scale,res=4,c=f_color).alpha(alpha)
                    #mesh = v.Mesh([coord[coords,:],[[0,1,2,3],[4,5,6,7],[0,3,7,4],[1,2,6,5],[0,1,5,4],[2,3,7,6]]],alpha=alpha,c=f_color).lw(1)
                else:
                    #mesh = v.Mesh([coord[coords,:],[[0,1,2,3],[4,5,6,7],[0,3,7,4],[1,2,6,5],[0,1,5,4],[2,3,7,6]]],alpha=alpha).lw(1)
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
                    #print(eq_dict)
                    print('j',j)
                    #print(eq)
                    #print(indx)
                    print(eq_els)
                    print(eq[indx])
                    eq_dict[j[0]] = eq[indx][0]
                    indx += 1

            #coord[i] += offset

            if element_type == 3:
                coords = vdu.get_coord_from_edof(edof[i,:],dof,3)
            elif element_type == 4:
                coords = vdu.get_coord_from_edof(edof[i,:],dof,4)
            elif element_type == 6:
            #for i in range(nel):
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
                #print('eq',eq)
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

        #if export is not None:
        #    v.io.write(mesh, export+".vtk")

        if render_nodes == True:
            if element_type == 3:
                nodes = vdu.get_node_elements(coord,scale,alpha,dof,bcPrescr,bc,bc_color,fPrescr,f,f_color,dofs_per_node=1)
            elif element_type == 4 or element_type == 6:
                nodes = vdu.get_node_elements(coord,scale,alpha,dof,bcPrescr,bc,bc_color,fPrescr,f,f_color,dofs_per_node=3)

            plot_window.meshes[plot_window.fig].extend(meshes)
            plot_window.nodes[plot_window.fig].extend(nodes)
            print("Adding mesh to figure ",plot_window.fig+1)
        else:
            plot_window.meshes[plot_window.fig].extend(meshes)
            print("Adding mesh to figure ",plot_window.fig+1)
        return meshes

    else:
        print("Invalid element type, please declare 'element_type'. The element types are:\n    1 - Spring\n    2 - Bar\n    3 - Flow\n    4 - Solid\n    5 - Beam\n    6 - Plate")
        sys.exit()

# Creates a deformed mesh for rendering, see element types above
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

    :param array edof: element topology by degrees of freedom [nel x (n_dofs_per_node)|(n_dofs_per_node+1)*n_nodes ]

    :return array el_values: Global scalar values for element
    """
    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window

    values = scalars

    if np.size(coord, axis = 1) == 1:
        print('only x')
        coord = np.append(coord, np.zeros((np.size(coord, axis = 0),1)), axis=1)
        coord = np.append(coord, np.zeros((np.size(coord, axis = 0),1)), axis=1)
        #np.c_[ coord, np.zeros((np.size(coord, axis = 0),1)), np.zeros((np.size(coord, axis = 0),1)) ]
    elif np.size(coord, axis = 1) == 2:
        print('only x & y')
        coord = np.append(coord, np.zeros((np.size(coord, axis = 0),1)), axis=1)
        #np.c_[ coord, np.zeros((np.size(coord, axis = 0),1)) ]

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
        print("Invalid element type, please declare 'element_type'. The element types are:\n    1 - Spring\n    2 - Bar\n    3 - Flow\n    4 - Solid\n    5 - Beam\n    6 - Plate")
        sys.exit()

    if a is None:
        a = np.zeros((nnode*ndof_per_n,1))
    elif element_type == 3:
        print('NOTE: Element type is flow, but deformation matrix given. Deformation was set to 0 for all DOFs')
        a = np.zeros((nnode*ndof_per_n,1))

    #if element_type == 1 and values is not None:
    #    print('NOTE: Colormapping for spring elements is not supported')

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
        ndof = np.size(dof, axis = 0)*np.size(dof, axis = 1)
        ncoord = np.size(coord, axis = 0)

        #nodes = []

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
                #def_coord[i] += offset
            elif element_type == 2 or element_type == 5:
                a_dx, a_dy, a_dz = vdu.get_a_from_coord(i,6,a,def_scale)

                x = coord[i,0]+a_dx
                y = coord[i,1]+a_dy
                z = coord[i,2]+a_dz

                def_coord[i] = [x,y,z]
                #def_coord[i] += offset

                #def_nodes.append(v.Sphere(c='white').scale(1.5*scale).pos([x,y,z]).alpha(alpha))

        nel = np.size(edof, axis = 0)
        elements = []
        res = 4

        el_values_array = np.zeros((1,4*4))[0,:]

        for i in range(nel):
            coord1,coord2 = vdu.get_coord_from_edof(edof[i,:],dof,element_type)
            '''
            if element_type == 1 and spring == True:
                #spring = v.Spring([coord[coord1,0],0,0],[coord[coord2,0],0,0],r=1.5*scale).alpha(alpha)
                spring = v.Spring([coord[coord1,0],coord[coord1,1],coord[coord1,2]],[coord[coord2,0],coord[coord2,1],coord[coord2,2]],r=1.5*scale,c=color).alpha(alpha)
                spring.name = f"Spring element {i+1}"
                elements.append(spring)
            elif element_type == 1 and spring == False:
                spring = v.Cylinder([[coord[coord1,0],coord[coord1,1],coord[coord1,2]],[coord[coord2,0],coord[coord2,1],coord[coord2,2]]],r=scale,res=4,c=color).alpha(alpha)
                spring.name = f"Spring element {i+1}"
                elements.append(spring)
            '''
            if element_type == 1 and spring == True:
                element = v.Spring([def_coord[coord1,0],def_coord[coord1,1],def_coord[coord1,2]],[def_coord[coord2,0],def_coord[coord2,1],def_coord[coord2,2]],r=scale*1.5,c=color).alpha(alpha)
                element.name = f"Spring element {i+1}"
                elements.append(element)

                #if values is not None:
                #    print('Colormapping not supported for springs')
                
                if values is not None:
                    el_values_array = []
                    for j in range(14):
                        el_values_array.append(values[i])
                    element.celldata[scalar_title] = el_values_array
                    element.cmap(colormap, scalar_title, on="cells", vmin=vmin, vmax=vmax)
                
            #if element_type == 1 and spring == False:
            #    element = v.Cylinder([[def_coord[coord1,0],def_coord[coord1,1],def_coord[coord1,2]],[def_coord[coord2,0],def_coord[coord2,1],def_coord[coord2,2]]],r=scale,res=4,c=color).alpha(alpha)
            #    element.name = f"Spring element {i+1}"
            #    elements.append(element)
            elif element_type == 2 or element_type == 1 and spring == False:
                element = v.Cylinder([[def_coord[coord1,0],def_coord[coord1,1],def_coord[coord1,2]],[def_coord[coord2,0],def_coord[coord2,1],def_coord[coord2,2]]],r=scale,res=4,c=color).alpha(alpha)
                #bar.info = f"Bar nr. {i}, at [{def_coord[coord1,0]+0.5*(def_coord[coord2,0]-def_coord[coord1,0])},{def_coord[coord1,1]+0.5*(def_coord[coord2,1]-def_coord[coord1,1])},{def_coord[coord1,2]+0.5*(def_coord[coord2,2]-def_coord[coord1,2])}]"
                if element_type == 2:
                    element.name = f"Bar element {i+1}"
                elif element_type == 1 and spring == False:
                    element.name = f"Spring element {i+1}"

                elements.append(element)

                if values is not None:
                    #bar.info = bar.info + f", max el. value {values[i]}"
                    #el_values_array[1] = values[i]
                    #el_values_array[3] = values[i]
                    #el_values_array[5] = values[i]
                    #el_values_array[7] = values[i]
                    #el_values_array[12] = values[i]
                    #el_values_array[13] = values[i]
                    #el_values_array[14] = values[i]
                    #el_values_array[15] = values[i]

                    #el_values_array[0] = values[i]
                    #el_values_array[2] = values[i]
                    #el_values_array[4] = values[i]
                    #el_values_array[6] = values[i]
                    #el_values_array[8] = values[i]
                    #el_values_array[9] = values[i]
                    #el_values_array[10] = values[i]
                    #el_values_array[11] = values[i]

                    #bar.cmap(colormap, el_values_array, on="points", vmin=vmin, vmax=vmax)
                    el_values_array = []
                    for j in range(6):
                        el_values_array.append(values[i])
                    element.celldata[scalar_title] = el_values_array
                    #element.cmap(colormap, [values[i],values[i],values[i],values[i],values[i],values[i]], on="cells", vmin=vmin, vmax=vmax)
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
                            #print('Beam: scalars',values)
                            el_value1 = values[nseg*i+j]
                            #beam.name = f"Beam nr. {i}, seg. {j}, element value: {np.round(values[nseg*i],2)}"
                            #beam.celldata["val"] = [values[nseg*i],values[nseg*i],values[nseg*i],values[nseg*i],values[nseg*i],values[nseg*i]]
                            el_value2 = values[nseg*i+j+1]

                            #print('beam',i+1,'segment',j+1,f'val {el_value1} and {el_value2}')
                            
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
                    #nodes.append(v.Sphere(c='white').scale(1.5*scale).pos([x,y,z]).alpha(alpha))
                plot_window.meshes[plot_window.fig].extend(elements)
                plot_window.nodes[plot_window.fig].extend(nodes)
            else:
                plot_window.meshes[plot_window.fig].extend(elements)

        return elements



    # Elements w/ a volume/surface (flow, solid & plate)
    elif element_type == 3 or element_type == 4 or element_type == 6:
        ex,ey,ez = cfc.coordxtr(edof,coord,dof)
        #if element_type == 3 or element_type == 4:
        #    ex,ey,ez = cfc.coordxtr(edof,coord,dof)
        #elif element_type == 6:
        #    ex,ey,ez = cfc.coordxtr(edof,coord,dof)
            #ez = np.zeros((nel,4))

        coord[:] += offset

        ed = cfc.extractEldisp(edof,a)
        if element_type == 3:
            if val != 'nodal_values_by_el':
                coord2, topo, node_dofs, a_node, node_scalars = vdu.convert_to_node_topo(edof,ex,ey,ez,ed,ignore_first=False,dofs_per_node=1)
            else:
                coord2, topo, node_dofs, a_node, node_scalars = vdu.convert_to_node_topo(edof,ex,ey,ez,ed,values,ignore_first=False,dofs_per_node=1)
            def_coord = coord2 + a_node*def_scale
            #for i in def_coord:
            #    i += offset
        elif element_type == 4:
            if val != 'nodal_values_by_el':
                coord2, topo, node_dofs, a_node, node_scalars = vdu.convert_to_node_topo(edof,ex,ey,ez,ed,ignore_first=False)
            else:
                coord2, topo, node_dofs, a_node, node_scalars = vdu.convert_to_node_topo(edof,ex,ey,ez,ed,values,ignore_first=False)
            def_coord = coord2 + a_node*def_scale
            #for i in def_coord:
            #    i += offset
        elif element_type == 6:
            coord2, topo, node_dofs, a_node, test = vdu.convert_to_node_topo(edof,ex,ey,ez,ed,ignore_first=False)
            def_coord = coord2
            def_coord[:,2] = a_node[:,0]*def_scale
            #for i in def_coord:
            #    i += offset

        #def_coord = coord2 + a_node*def_scale

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
                #print(val)
                #vmin, vmax = np.min(values), np.max(values)
                
                el_values = vdu.convert_el_values(edof,values)
                mesh.celldata[scalar_title] = el_values

                mesh.cmap(colormap, scalar_title, on="cells", n=colors, vmax=vmax, vmin=vmin)
            
            elif val and val == 'nodal_values_by_el':
                #print(val)
                #vmin, vmax = np.min(values), np.max(values)
                nodal_values = vdu.convert_nodal_values(edof,topo,dof,values)
                #nodal_values = vdu.convert_a(coord2,coord,nodal_values,1)
                mesh.pointdata[scalar_title] = nodal_values
                #mesh.pointdata["val"] = node_scalars
                print(ug.celldata.keys())
                #nodal_values = vdu.convert_nodal_values(edof,dof,coord,coord2,values)
                mesh.cmap(colormap, scalar_title, on="points", n=colors, vmax=vmax, vmin=vmin)


            elif val and val == 'nodal_values':
                print(val)
                #values = vdu.convert_nodal_values(edof,topo,dof,values)
                #vmin, vmax = np.min(values), np.max(values)
                mesh.pointdata[scalar_title] = values
                mesh.cmap(colormap, scalar_title, on="points", n=colors, vmax=vmax, vmin=vmin)
                #ug.pointdata["val"] = values
                #nodal_values = vdu.convert_nodal_values(edof,dof,coord,coord2,values)
                #mesh.cmap(colormap, values, on="points", vmin=vmin, vmax=vmax)

        elif element_type == 6:
            if val and val == 'el_values':
                #print(val)
                #vmin, vmax = np.min(values), np.max(values)
                
                #el_values = vdu.convert_el_values(edof,values)
                mesh.celldata[scalar_title] = values

                mesh.cmap(colormap, scalar_title, on="cells", n=colors, vmax=vmax, vmin=vmin)
        
        '''
        print('Number of topo cells: ',np.size(topo, axis=0))
        print('Number of mesh cells: ',np.size(mesh.faces(), axis=0))

        print('Number of topo points: ',np.size(coord2, axis=0))
        print('Number of mesh points: ',np.size(mesh.points(), axis=0))

        print('Number of Vedo mesh coordinates: ',mesh.N())
        '''

        if only_ret == False:

            if render_nodes == True:
                #nodes = vdu.get_node_elements(def_coord,scale,alpha)
                if element_type == 3:
                    nodes = vdu.get_node_elements(def_coord,scale,alpha,dof,dofs_per_node=1)
                elif element_type == 4 or element_type == 6:
                    nodes = vdu.get_node_elements(def_coord,scale,alpha,dof,dofs_per_node=3)
                    #nodes.append(v.Sphere(c='white').scale(1.5*scale).pos([x,y,z]).alpha(alpha))
                #plot_window.meshes[plot_window.fig] += mesh
                #plot_window.nodes[plot_window.fig] += nodes
                plot_window.meshes[plot_window.fig].extend(elements)
                plot_window.nodes[plot_window.fig].extend(nodes)
                #print("Meshes are ",np.size(plot_window.meshes, axis=0),"X",np.size(plot_window.meshes, axis=1))
                #print("Adding mesh to figure ",plot_window.fig+1)
            else:
                plot_window.meshes[plot_window.fig].append(mesh)
                #print("Meshes are ",np.size(plot_window.meshes, axis=0),"X",np.size(plot_window.meshes, axis=1))
                #print("Adding mesh to figure ",plot_window.fig+1)

        return mesh



'''
def test(edof,
    ex,
    ey,
    ez,
    a=None,
    el_values=None,
    colormap='jet',
    scale=0.02,
    alpha=1,
    def_scale=1,
    nseg=2,
    color='white',
    offset = [0, 0, 0],
    merge=False,
    t=None):

    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window


    coord, topo, node_dofs = convert_to_node_topo(edof,ex,ey,ez,ignore_first=False)
    #return coords, topo, node_dofs

    if a is None:
        a = np.zeros([np.size(coord, axis = 0)*np.size(coord, axis = 1),1])

    nnode = np.size(coord, axis = 0)
    nel = np.size(topo, axis = 0)
    ndof = np.size(topo, axis = 1)
    print(np.size(topo, axis = 0))
    print(np.size(topo, axis = 1))
    print(np.size(a, axis = 0))
    #print(np.size(a, axis = 1))
    #print(coord[0][0])
    #print(coord[0][1])
    #print(coord[0][2])

    """
    if ct == vtk.VTK_HEXAHEDRON:
                    cell = vtk.vtkHexahedron()
                elif ct == vtk.VTK_TETRA:
                    cell = vtk.vtkTetra()
                elif ct == vtk.VTK_VOXEL:
                    cell = vtk.vtkVoxel()
                elif ct == vtk.VTK_WEDGE:
                    cell = vtk.vtkWedge()
                elif ct == vtk.VTK_PYRAMID:
                    cell = vtk.vtkPyramid()
                elif ct == vtk.VTK_HEXAGONAL_PRISM:
                    cell = vtk.vtkHexagonalPrism()
                elif ct == vtk.VTK_PENTAGONAL_PRISM:
                    cell = vtk.vtkPentagonalPrism()
    """

    #def_nodes = []
    def_coord = np.zeros([nnode,3])
    #celltype = [vtk.VTK_HEXAHEDRON] * nel

    #pdata = np.zeros((nnode), dtype=float)

    for i in range(nnode):
        #a_dx, a_dy, a_dz = get_a_from_coord(i,3,a,def_scale)
        #x = coord[i][0]+a_dx
        #y = coord[i][1]+a_dy
        #z = coord[i][2]+a_dz
        #def_coord[i] = [x,y,z]
        def_coord[i,0] = a[i*3]
        def_coord[i,1] = a[i*3+1]
        def_coord[i,2] = a[i*3+2]


    #meshes = []
    #nel = np.size(edof, axis = 0)

    #for i, dofs in enumerate(node_dofs):
        #v = u0[dofs-1]
        #pdata[i] = np.linalg.norm(v)

    #ugrid = v.UGrid([coord, topo, celltype])
    print(coord[0])
    print(topo)
    mesh = v.Mesh([coord, topo[0]]).lw(1)
    #ugrid.pointdata["mag"] = pdata

    #mesh = ugrid.tomesh()
    
    #for i in range(nel):
    #coords = get_coord_from_edof(edof[i,:],dof,4)

    #print(topo)

    #mesh = v.UGrid([def_coord, topo, [10]]).lw(10)
    #[[0,1,2,3],[4,5,6,7],[0,3,7,4],[1,2,6,5],[0,1,5,4],[2,3,7,6]],
    mesh.color(c='red')
    #mesh.name = f"Mesh nr. {i+1}"
    #meshes.append(mesh)

            


    #plot_window.meshes[plot_window.fig].extend(mesh)
    plot_window.meshes[plot_window.fig] += [mesh]
    print("Adding mesh to figure ",plot_window.fig+1)
'''




# Creates a deformed mesh for rendering, see element types above
#def draw_vectors(edof,coord,dof,element_type,vect):




### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
# Functions regarding animations

'''
def draw_displaced_mesh(
    # Main input
    edof,
    coord,
    dof,
    element_type,
    a=None,
    values=None,
    
    # Other parameters
    scale=0.02,
    alpha=1,
    def_scale=1,
    render_nodes=False,
    color='white',
    offset = [0, 0, 0],
    only_ret=False,
    lines=False,

    # Input for colormapping
    colormap='jet',
    colors = 256,
    vmax=None,
    vmin=None,

    # Element-specific input
    spring = True,
    nseg=2
    ):
'''

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

    :param array edof: element topology by degrees of freedom [nel x (n_dofs_per_node)|(n_dofs_per_node+1)*n_nodes ]

    :return array el_values: Global scalar values for element
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

    print(values)

    if element_type == 1 or element_type == 3:
        ndof_per_n = 1
    elif element_type == 2 or element_type == 4 or element_type == 6:
        ndof_per_n = 3
    elif element_type == 5:
        ndof_per_n = 6

    if a is None:
        a = np.zeros((nnode*ndof_per_n,1))
    elif element_type == 3:
        print('NOTE: Element type is flow, but deformation matrix given. Deformation was set to 0 for all DOFs')
        a = np.zeros((nnode*ndof_per_n,1))

    #v.settings.immediateRendering = True

    #camera = dict(viewAngle=30)
    if negative == True:
        timesteps = np.arange(0, 1+1/(steps), 1/(steps))
        unique_timesteps = timesteps
        timesteps = np.append(timesteps, np.flip(np.arange(0, 1, 1/(steps))))
        timesteps = np.append(timesteps, np.flip(np.arange(-1, 0, 1/(steps))))
        unique_timesteps = np.append(unique_timesteps, np.flip(np.arange(-1, 0, 1/(steps))))
        timesteps = np.append(timesteps, np.arange(-1+1/(steps), 0, 1/(steps)))
        #print('To be appended:',np.arange(-1+1/(steps), 0, 1/(steps)))
        #print('Looping incl. negative deformation, timesteps:',timesteps)
        #sys.exit()
    elif loop == True:
        timesteps = np.arange(0, 1+1/(steps), 1/(steps))
        unique_timesteps = timesteps
        timesteps = np.append(timesteps, np.flip(np.arange(0+1/(steps), 1, 1/(steps))))
        #print('No looping, timesteps:',timesteps)
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
        print('Creating keyframe', it+1)
        print('t',t)
        mesh = None
        #mesh = draw_displaced_mesh(edof,coord,dof,element_type,a*t,values*t,scale=scale,alpha=alpha,def_scale=def_scale,colormap=colormap,colors=colors,vmax=vmax,vmin=vmin,only_ret=True)
        
        if values is not None and animate_colormap == True:
            
            if t >= 0:
                print('animate_colormap & t>=0')
                mesh = draw_displaced_mesh(edof,coord,dof,element_type,a*t,values*t,scale=scale,alpha=alpha,def_scale=def_scale,colormap=colormap,colors=colors,vmax=vmax,vmin=vmin,scalar_title=scalar_title,nseg=nseg,spring=spring,only_ret=True)
            
            else:
                print('animate_colormap & t<0')
                mesh = draw_displaced_mesh(edof,coord,dof,element_type,a*t,-values*t,scale=scale,alpha=alpha,def_scale=def_scale,colormap=colormap,colors=colors,vmax=vmax,vmin=vmin,scalar_title=scalar_title,nseg=nseg,spring=spring,only_ret=True)
        
        elif values is not None and animate_colormap == False:
            print('no animate_colormap')
            mesh = draw_displaced_mesh(edof,coord,dof,element_type,a*t,values,scale=scale,alpha=alpha,def_scale=def_scale,colormap=colormap,colors=colors,vmax=vmax,vmin=vmin,scalar_title=scalar_title,nseg=nseg,spring=spring,only_ret=True)
        
        else:
            print('no scalars')
            mesh = draw_displaced_mesh(edof,coord,dof,element_type,a*t,scale=scale,alpha=alpha,def_scale=def_scale,nseg=nseg,spring=spring,only_ret=True)
        
        if element_type == 1 or element_type == 2 or element_type == 5:
            
            mesh = v.merge(mesh)
            #mesh = v.Assembly(mesh)

        keyframes[np.round(t,3)] = mesh
        
        it += 1
    if only_export == False:
        plot_window.keyframes[plot_window.fig].append(keyframes)
        plot_window.keyframe_dict[plot_window.fig].append(keyframe_dict)

    if export == True:
        for key, val in keyframe_dict.items():
            #if element_type == 1 or element_type == 2 or element_type == 5:
            #    mesh = keyframes[val].unpack()
            #else:
            mesh = keyframes[val]
            output = file+f'_{int(key)}'
            export_vtk(output,mesh)
    
    #print('keyframe_dict',keyframe_dict)
    #res = list(sorted({ele for val in keyframe_dict.values() for ele in val}))
    #print(unique_timesteps)
    #print(timesteps)
    #print(start)
    #print(end)
    #print(steps)
    
    #if self.elements != None:
    #    element_type = 1
    #elif self.mesh != None:
    #    element_type = 2

    #dt = 0.1

    #t = np.arange(0.0, 10.0, dt)

    #pb = v.ProgressBar(0, len(t), c="b")

    #nsteps = np.size(timesteps,0)
    """
    it = 0
    if element_type == 4:
        #ncoord = np.size(coord, axis = 0)
        #nel = np.size(edof, axis = 0)


        if values is not None:
            if vmax is None and vmin is None:
                vmin, vmax = np.min(values), np.max(values)
            elif vmax is None:
                vmax = np.max(values)
            elif vmin is None:
                vmin = np.min(values)

        #keyframes = []

        for t in unique_timesteps:
            '''
            def draw_displaced_mesh(
                # Main input
                edof,
                coord,
                dof,
                element_type,
                a=None,
                values=None,
                
                # Other parameters
                scale=0.02,
                alpha=1,
                def_scale=1,
                render_nodes=False,
                color='white',
                offset = [0, 0, 0],
                only_ret=False,
                lines=False,

                # Input for colormapping
                colormap='jet',
                colors = 256,
                vmax=None,
                vmin=None,

                # Element-specific input
                spring = True,
                nseg=2
                ):
            '''

            mesh = draw_displaced_mesh(edof,coord,dof,element_type,a*t,values*t,scale=scale,alpha=alpha,def_scale=def_scale,colormap=colormap,colors=colors,vmax=vmax,vmin=vmin,only_ret=True)
            #keyframes.append(mesh)
            keyframes[t] = mesh

            if export == True:
                output = file+f'_{int(it)}'
                export_vtk(output,mesh)
            #plt.render(resetcam=True)
            #plt.show(mesh).interactive().close()
            it += 1
        plot_window.keyframes[plot_window.fig].extend(keyframes)
        
        

        #print('keyframes',keyframes)
        
        #plot_window.
        
        #def_coord = np.zeros([ncoord,3])
        #def_coord = np.zeros((ncoord,3,steps))


        #meshes = []


        #for i in range(nel):
        #    coords = vdu.get_coord_from_edof(edof[i,:],dof,4)


        #    mesh = v.Mesh([def_coord[coords,:],[[0,1,2,3],[4,5,6,7],[0,3,7,4],[1,2,6,5],[0,1,5,4],[2,3,7,6]]],alpha=alpha).lw(1)
            #mesh.info = f"Mesh nr. {i}"
            #mesh.name = f"Mesh nr. {i+1}"
        #    meshes.append(mesh)

        #mesh = v.merge(meshes)

        #plt = v.show(mesh, axes=4, interactive=0)

        
        #plt += mesh


        

        
        #for t in timesteps:

            
        
        for i in range(0, ncoord):
            a_dx, a_dy, a_dz = vdu.get_a_from_coord(i,3,a,def_scale)

            a_dx = t*a_dx
            a_dy = t*a_dy
            a_dz = t*a_dz

            x_step = a_dx/(steps)
            y_step = a_dy/(steps)
            z_step = a_dz/(steps)

            for j in range(0, steps):

                x = coord[i,0]+x_step*j
                y = coord[i,1]+y_step*j
                z = coord[i,2]+z_step*j

                def_coord[i,:] = [x,y,z]
        

        '''
        ex,ey,ez = cfc.coordxtr(edof,coord,dof)
        ed = cfc.extractEldisp(edof,a)
        coord2, topo, node_dofs, a_node, node_scalars = vdu.convert_to_node_topo(edof,ex,ey,ez,ed,ignore_first=False,dofs_per_node=3)
        
        new_coord = coord2 + a_node*def_scale*t

        ct = vtk.VTK_HEXAHEDRON

        celltypes = [ct] * nel

        ug=v.UGrid([new_coord, topo, celltypes])
        ug.points(new_coord)
        
        mesh = ug.tomesh().lw(1).alpha(alpha)

        # Element values
        el_values = vdu.convert_el_values(edof,values)
        mesh.celldata["val"] = el_values

        mesh.cmap(colormap, "val", on="cells", vmin=vmin*t, vmax=vmax*t)
        #meshes = []
        '''
        
        #vmin, vmax = np.min(el_values), np.max(el_values)
        for i in range(nel):
            coords = vdu.get_coord_from_edof(edof[i,:],dof,4)

            mesh = v.Mesh([def_coord[coords,:],[[0,1,2,3],[4,5,6,7],[0,3,7,4],[1,2,6,5],[0,1,5,4],[2,3,7,6]]],alpha=alpha).lw(1)
            meshes.append(mesh)

            if values is not None and np.size(values, axis = 1) == 1:
                el_values_array = np.zeros((1,6))[0,:]
                el_values_array[0] = values[i]*t
                el_values_array[1] = values[i]*t
                el_values_array[2] = values[i]*t
                el_values_array[3] = values[i]*t
                el_values_array[4] = values[i]*t
                el_values_array[5] = values[i]*t
                #if title is not None:
                #    mesh.cmap(colormap, el_values_array, on="cells", vmin=vmin, vmax=vmax).addScalarBar(title=title,horizontal=True,useAlpha=False,titleFontSize=16)
                #else:
                mesh.cmap(colormap, el_values_array, on="cells", vmin=vmin, vmax=vmax)
            #for j in range(steps):
                #mesh = v.Mesh([def_coord[coords,:,j],[[0,1,2,3],[4,5,6,7],[0,3,7,4],[1,2,6,5],[0,1,5,4],[2,3,7,6]]],alpha=alpha)
                #meshes[i,j] = mesh
            #mesh = v.Mesh([def_coord[coords,:],[[0,1,2,3],[4,5,6,7],[0,3,7,4],[1,2,6,5],[0,1,5,4],[2,3,7,6]]],alpha=alpha)
            #mesh.info = f"Mesh nr. {i}"
            #meshes.append(mesh)
            mesh = v.merge(meshes)
        """

            #plt.clear()
            #plt += mesh
            #if export == True:
            #    output = file+f'_{int(10*t)}'
            #    export_vtk(output,mesh)
            #plt.render(resetcam=True)
            #plt.show(mesh).interactive().close()
            


        #v.interactive().close()
        #v.interactive().close()

        #v.settings.immediateRendering = False

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
# Functions for importing/exporting

def import_mat(file,list=None):
    """
    Routine for importing from MATLAB.

    :param array edof: element topology by degrees of freedom [nel x (n_dofs_per_node)|(n_dofs_per_node+1)*n_nodes ]

    :return array el_values: Global scalar values for element
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

    :param array edof: element topology by degrees of freedom [nel x (n_dofs_per_node)|(n_dofs_per_node+1)*n_nodes ]

    :return array el_values: Global scalar values for element
    """
    if isinstance(meshes, list):
        mesh = v.merge(meshes)
    else:
        mesh = meshes
    
    #print('meshes')
    #print(meshes)
    #print('mesh')
    #print(mesh)
    #mesh = v.merge(meshes)
    #for i in range(len(meshes)):
    v.io.write(mesh, file+".vtk")

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
# Functions for handling rendering


#def figure(fig=None):
def figure(fig,bg='white',flat=False,hover=False,axes=False):
    """
    Routine for choosing what figure to add objects to, and creating necessary lists etc.

    :param array edof: element topology by degrees of freedom [nel x (n_dofs_per_node)|(n_dofs_per_node+1)*n_nodes ]

    :return array el_values: Global scalar values for element
    """
    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window

    #if fig == None:
    #    fig = plot_window.fig + 1

    if fig < 1:
        print("Please give a positive integer (> 0)")
        sys.exit()
    else:
        plot_window.fig = fig - 1

    print("Selecting figure ",fig)
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

            #print(plot_window.mode_2D)
            
        #plot_window.mode_2D.append(False)
    elif fig == 1:
        plot_window.bg.append(bg)
        plot_window.mode_2D.append(False)
        plot_window.mode_hover.append(False)
        plot_window.mode_anim.append(False)

    plot_window.bg[fig-1] = bg

    if flat == True:
        plot_window.mode_2D[fig-1] = True
        #print(plot_window.mode_2D)
        #plot_window.mode_2D.append(True)
        #print(plot_window.mode_2D)
        #sys.exit()
    #else:
    #    plot_window.mode_2D.append(False)

    if hover == True:
        #print(plot_window.mode_2D)
        #plot_window.mode_hover.append(True)
        plot_window.mode_hover[fig-1] = True


# Lägg till figurnummer här???
# Start Calfem-vedo visualization
def show_and_wait():
    """
    Routine for showing figure(s).

    :param array edof: element topology by degrees of freedom [nel x (n_dofs_per_node)|(n_dofs_per_node+1)*n_nodes ]

    :return array el_values: Global scalar values for element
    """
    app = init_app()
    plot_window = VedoPlotWindow.instance().plot_window
    plot_window.render()
