# -*- coding: utf-8 -*-

import os, sys

from PyQt5.QtCore import pyqtSlot, pyqtSignal, QThread
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QDialog, QWidget, QMainWindow, QFileDialog
from PyQt5.QtGui import QIcon
from PyQt5 import uic

#from scipy.interpolate import spline
import calfem.geometry as cfg
import calfem.vis_mpl as cfv
import numpy as np

from matplotlib.lines import Line2D
from matplotlib.artist import Artist
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt

def dist(x, y):
    """
    Return the distance between two points.
    """
    d = x - y
    return np.sqrt(np.dot(d, d))


def dist_point_to_segment(p, s0, s1):
    """
    Get the distance of a point to a segment.
      *p*, *s0*, *s1* are *xy* sequences
    This algorithm from
    http://geomalgorithms.com/a02-_lines.html
    """
    v = s1 - s0
    w = p - s0
    c1 = np.dot(w, v)
    if c1 <= 0:
        return dist(p, s0)
    c2 = np.dot(v, v)
    if c2 <= c1:
        return dist(p, s1)
    b = c1 / c2
    pb = s0 + b * v
    return dist(p, pb)

class GeomInteracter:
    """
    An path editor.

    Key-bindings

      't' toggle vertex markers on and off.  When vertex markers are on,
          you can move them, delete them


    """

    showverts = True
    epsilon = 5  # max pixel distance to count as a vertex hit

    def __init__(self, geometry):

        self.geometry = geometry

        self.ax = pathpatch.axes
        canvas = self.ax.figure.canvas
        self.pathpatch = pathpatch
        self.pathpatch.set_animated(True)

        self.lines = []

        for curve in self.geometry.curves:
            c = self.geometry.curves[curve]
            x = []
            y = []
            if c[0] == "Spline":
                if len(c[1])==2:
                    x.append(self.geometry.points[c[1][0]][0][0])
                    x.append(self.geometry.points[c[1][1]][0][0])
                    y.append(self.geometry.points[c[1][0]][0][1])
                    y.append(self.geometry.points[c[1][1]][0][1])
                    l, = ax.plot(x, y, marker='o', markerfacecolor='r', animated=True)
                    self.lines.append(l)
                elif len(c[1])==3:
                    x.append(self.geometry.points[c[1][0]][0][0])
                    x.append(self.geometry.points[c[1][1]][0][0])
                    x.append(self.geometry.points[c[1][2]][0][0])
                    y.append(self.geometry.points[c[1][0]][0][1])
                    y.append(self.geometry.points[c[1][1]][0][1])
                    y.append(self.geometry.points[c[1][2]][0][1])
                    l, = ax.plot(x, y, marker='o', markerfacecolor='r', animated=True)
                    self.lines.append(l)
                else:
                    pass


        self._ind = None  # the active vert

        canvas.mpl_connect('draw_event', self.draw_callback)
        canvas.mpl_connect('button_press_event', self.button_press_callback)
        canvas.mpl_connect('key_press_event', self.key_press_callback)
        canvas.mpl_connect('button_release_event', self.button_release_callback)
        canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        self.canvas = canvas

    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        #self.ax.draw_artist(self.pathpatch)
        for line in self.lines:
            self.ax.draw_artist(line)
        self.canvas.blit(self.ax.bbox)

    def pathpatch_changed(self, pathpatch):
        """This method is called whenever the pathpatch object is called."""
        # only copy the artist props to the line (except visibility)
        vis = self.line.get_visible()
        plt.Artist.update_from(self.line, pathpatch)
        self.line.set_visible(vis)  # don't use the pathpatch visibility state

    def get_ind_under_point(self, event):
        """
        Return the index of the point closest to the event position or *None*
        if no point is within ``self.epsilon`` to the event position.
        """
        # display coords
        xy = np.asarray(self.pathpatch.get_path().vertices)
        xyt = self.pathpatch.get_transform().transform(xy)
        xt, yt = xyt[:, 0], xyt[:, 1]
        d = np.sqrt((xt - event.x)**2 + (yt - event.y)**2)
        ind = d.argmin()

        if d[ind] >= self.epsilon:
            ind = None

        return ind

    def button_press_callback(self, event):
        """Callback for mouse button presses."""
        if not self.showverts:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        self._ind = self.get_ind_under_point(event)

    def button_release_callback(self, event):
        """Callback for mouse button releases."""
        if not self.showverts:
            return
        if event.button != 1:
            return
        self._ind = None

    def key_press_callback(self, event):
        """Callback for key presses."""
        if not event.inaxes:
            return
        if event.key == 't':
            self.showverts = not self.showverts
            self.line.set_visible(self.showverts)
            if not self.showverts:
                self._ind = None

        self.canvas.draw()

    def motion_notify_callback(self, event):
        """Callback for mouse movements."""
        if not self.showverts:
            return
        if self._ind is None:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        x, y = event.xdata, event.ydata

        vertices = self.pathpatch.get_path().vertices

        vertices[self._ind] = x, y
        self.line.set_data(zip(*vertices))

        self.canvas.restore_region(self.background)
        #self.ax.draw_artist(self.pathpatch)
        self.ax.draw_artist(self.line)
        self.canvas.blit(self.ax.bbox)

class PolygonInteractor(object):
    """
    A polygon editor.

    Key-bindings

      't' toggle vertex markers on and off.  When vertex markers are on,
          you can move them, delete them

      'd' delete the vertex under point

      'i' insert a vertex at point.  You must be within epsilon of the
          line connecting two existing vertices

    """

    showverts = True
    epsilon = 5  # max pixel distance to count as a vertex hit

    def __init__(self, ax, poly):
        if poly.figure is None:
            raise RuntimeError('You must first add the polygon to a figure '
                               'or canvas before defining the interactor')
        self.ax = ax
        canvas = poly.figure.canvas
        self.poly = poly

        x, y = zip(*self.poly.xy)
        self.line = Line2D(x, y, marker='o', markerfacecolor='r', animated=True)
        self.line.set_color('g')
        self.ax.add_line(self.line)

        self.cid = self.poly.add_callback(self.poly_changed)
        self._ind = None  # the active vert

        canvas.mpl_connect('draw_event', self.draw_callback)
        canvas.mpl_connect('button_press_event', self.button_press_callback)
        canvas.mpl_connect('key_press_event', self.key_press_callback)
        canvas.mpl_connect('button_release_event', self.button_release_callback)
        canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        self.canvas = canvas

    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.ax.draw_artist(self.poly)
        self.ax.draw_artist(self.line) 
        # do not need to blit here, this will fire before the screen is
        # updated

    def poly_changed(self, poly):
        'this method is called whenever the polygon object is called'
        # only copy the artist props to the line (except visibility)
        vis = self.line.get_visible()
        Artist.update_from(self.line, poly)
        self.line.set_visible(vis)  # don't use the poly visibility state

    def get_ind_under_point(self, event):
        'get the index of the vertex under point if within epsilon tolerance'

        # display coords
        xy = np.asarray(self.poly.xy)
        xyt = self.poly.get_transform().transform(xy)
        xt, yt = xyt[:, 0], xyt[:, 1]
        d = np.hypot(xt - event.x, yt - event.y)
        indseq, = np.nonzero(d == d.min())
        ind = indseq[0]

        if d[ind] >= self.epsilon:
            ind = None

        return ind

    def button_press_callback(self, event):
        'whenever a mouse button is pressed'
        if not self.showverts:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        self._ind = self.get_ind_under_point(event)

    def button_release_callback(self, event):
        'whenever a mouse button is released'
        if not self.showverts:
            return
        if event.button != 1:
            return
        self._ind = None

    def key_press_callback(self, event):
        'whenever a key is pressed'
        print("key_press_callback")

        if not event.inaxes:
            return
        if event.key == 't':
            self.showverts = not self.showverts
            self.line.set_visible(self.showverts)
            if not self.showverts:
                self._ind = None
        elif event.key == 'd':
            ind = self.get_ind_under_point(event)
            if ind is not None:
                self.poly.xy = np.delete(self.poly.xy,
                                         ind, axis=0)
                self.line.set_data(zip(*self.poly.xy))
        elif event.key == 'i':
            xys = self.poly.get_transform().transform(self.poly.xy)
            p = event.x, event.y  # display coords
            for i in range(len(xys) - 1):
                s0 = xys[i]
                s1 = xys[i + 1]
                d = dist_point_to_segment(p, s0, s1)
                if d <= self.epsilon:
                    self.poly.xy = np.insert(
                        self.poly.xy, i+1,
                        [event.xdata, event.ydata],
                        axis=0)
                    self.line.set_data(zip(*self.poly.xy))
                    break
        if self.line.stale:
            self.canvas.draw_idle()

    def motion_notify_callback(self, event):
        'on mouse movement'
        if not self.showverts:
            return
        if self._ind is None:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        x, y = event.xdata, event.ydata

        self.poly.xy[self._ind] = x, y
        if self._ind == 0:
            self.poly.xy[-1] = x, y
        elif self._ind == len(self.poly.xy) - 1:
            self.poly.xy[0] = x, y
        self.line.set_data(zip(*self.poly.xy))

        self.canvas.restore_region(self.background)
        self.ax.draw_artist(self.poly)
        self.ax.draw_artist(self.line)
        self.canvas.blit(self.ax.bbox)

class MainWindow(QMainWindow):
    """MainWindow-klass som hanterar vårt huvudfönster"""

    def __init__(self, app, geometry):
        """Class constructor"""
        
        super().__init__()

        # --- Lagra en referens till applikationsinstansen i klassen
        
        self.app = app

        # --- Geometri

        self.geometry = geometry
                       
        # --- Läs in gränssnitt från fil
        
        # Läs in gränssnitt från fil

        uic.loadUi("cfeditor.ui", self)
              
        # --- Koppla kontroller till händelsemetoder

        self.fig = cfv.figure()

        self.fw = cfv.figure_widget(self.fig)
        self.fw.setFocusPolicy(Qt.StrongFocus)

        self.setCentralWidget(self.fw)

        self.axes = self.fig.gca()
        self.axes.set_xlim(0, 300)
        self.axes.set_ylim(0, 300)

        #print(self.geometry.points)
        #print(self.geometry.curves)
        


        # --- Se till att visa fönstret

       
        self.fw.show()
        self.show()
        self.raise_()
        
        # --- Se till att vi har en initierad modell  

    def on_plt_click(self, event):
        print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
            ('double' if event.dblclick else 'single', event.button,
            event.x, event.y, event.xdata, event.ydata))

    def on_plt_keypress(self, event):
        print("keypres")

def cfedit(geometry):
    app = QApplication(sys.argv)
    window = MainWindow(app, geometry)
    sys.exit(app.exec_())     
       

if __name__ == '__main__':

    # ---- Create Geometry ------------------------------------------------------

    g = cfg.geometry()

    # Add Points:

    points = [
        [0,0], 
        [0,100], 
        [0,150], 
        [100,0], 
        [150,0], 
        [100,-100], 
        [150,-100]
    ]

    for p in points:
        g.point(p)

    # Add Splines:

    g.spline([1,2], marker=2, el_on_curve=4)
    g.spline([3,4], el_on_curve=4)
    g.circle([1,0,3], el_on_curve = 10)
    g.circle([2,0,4], el_on_curve = 10)
    g.spline([3,5], el_on_curve = 6)
    g.spline([5,6], marker=3, el_on_curve = 4)
    g.spline([6,4], el_on_curve = 6)

    # Add Surfaces:
    #
    # When we set markers for surfaces, and have 2D elements, we can find which 
    # region an element is in via the list 'elementmarkers', which is returned by 
    # GmshMesher.create()

    g.structuredSurface([0,2,1,3], marker = 10)
    g.structuredSurface([1,4,5,6], marker = 11)

    cfedit(g)

