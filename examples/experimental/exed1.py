import sys
from math import sqrt

from PyQt5.QtCore import pyqtSlot, pyqtSignal, QThread
from PyQt5.QtWidgets import QApplication, QDialog, QWidget, QMainWindow, QFileDialog, QGraphicsScene, QGraphicsView, QGraphicsItem, QGraphicsEllipseItem, QGraphicsLineItem
from PyQt5.QtGui import QPen, QBrush
from PyQt5.Qt import Qt
from PyQt5 import uic

import numpy as np

import calfem.geometry as cfg
#import calfem.mesh as cfm
#import calfem.vis as cfv
#import calfem.utils as cfu
#import calfem.core as cfc

class NodeItem(QGraphicsEllipseItem):
    def __init__(self, x, y, width, height, parent = None):
        super().__init__(x, y, width, height, parent)

    def sceneEvent(self, event):
        print(event.type())
        return False


class MainWindow(QMainWindow):
    """Main window class of our UI"""
    def __init__(self, g):
        """Constructor"""
        super(MainWindow, self).__init__()

        self.g = g

        # Load user interface from UI-file

        uic.loadUi('exed1.ui', self)

        self.init_scene()

    def calc_bounding_box(self):
        """Calc geometry bounding box"""

        self.max_x = -1e300
        self.max_y = -1e300
        self.min_x = 1e300
        self.min_y = 1e300

        for p in self.g.points.values():
            if p[0][0]>self.max_x:
                self.max_x = p[0][0]
            if p[0][0]<self.min_x:
                self.min_x = p[0][0]
            if p[0][1]>self.max_y:
                self.max_y = p[0][1]
            if p[0][1]<self.min_y:
                self.min_y = p[0][1]

        self.bw = self.max_x - self.min_x
        self.bh = self.max_y - self.min_y
        self.hs = (self.bw+self.bh)*0.5*0.075


    def init_scene(self):
        """Setup scene"""
        self.scene = QGraphicsScene()
        self.green_brush = QBrush(Qt.green)
        self.gray_brush = QBrush(Qt.gray)
        self.white_brush = QBrush(Qt.white)
        self.pen = QPen(Qt.gray)
        self.pen.setWidth(4)
        self.pen.setCosmetic(True)
        self.line_pen = QPen(Qt.red)
        self.line_pen.setWidth(4)
        self.line_pen.setCosmetic(True)

        self.calc_bounding_box()

        for c in self.g.curves.values():
            if c[0] == 'Spline':
                i0 = c[1][0]
                i1 = c[1][1]
                p0 = self.g.points[i0][0]
                p1 = self.g.points[i1][0]
                gi = self.scene.addLine(p0[0], p0[1], p1[0], p1[1], self.pen)

        for p in self.g.points.values():
            gi = QGraphicsEllipseItem(-self.hs/2, -self.hs/2, self.hs, self.hs)
            gi.setBrush(self.white_brush)
            gi.setPen(self.pen)
            gi.setPos(p[0][0], p[0][1])
            gi.setFlag(QGraphicsItem.ItemIsMovable)
            gi.setFlag(QGraphicsItem.ItemIsSelectable)
            self.scene.addItem(gi)
            #gi = self.scene.addEllipse(-self.hs/2, -self.hs/2, self.hs, self.hs, self.pen, self.white_brush)


        print(self.scene.sceneRect())

        #self.scene.setSceneRect(self.min_x, self.min_y, self.max_x-self.min_x, self.max_y-self.min_y)

        self.graphics_view.setScene(self.scene)
        #self.graphics_view.fitInView(self.min_x, self.min_y, self.max_x-self.min_x, self.max_y-self.min_y)
        self.graphics_view.scale(200.0, 200.0)


        #ellipse = self.scene.addEllipse(20,20, 200,200, self.pen, self.greenBrush)
        #rect = self.scene.addRect(self.min_x, self.min_y, self.max_x-self.min_x, self.max_y-self.min_y, self.pen)

        #ellipse.setFlag(QGraphicsItem.ItemIsMovable)
        #rect.setFlag(QGraphicsItem.ItemIsMovable)
        #ellipse.setFlag(QGraphicsItem.ItemIsSelectable)

        self.scene.focusItemChanged.connect(self.on_focus_item_changed)
        self.scene.selectionChanged.connect(self.on_selection_changed)

    def on_focus_item_changed(self, new_focus_item, old_focus_item, reason):
        print("on_focus_item_changed", new_focus_item, reason)

    def on_selection_changed(self):
        if len(self.scene.selectedItems())>0:
            print("on_selection_changed", self.scene.selectedItems()[0].pos())

    def on_item_mouse_move(self):
        print("on_item_mouse_move")

def edit_geometry(g):
    app = QApplication(sys.argv)
    widget = MainWindow(g)
    widget.show()
    sys.exit(app.exec_())


if __name__ == "__main__":

    g = cfg.Geometry()  # Create a GeoData object that holds the geometry.

    # Add points:
    #  The first parameter is the coordinates. These can be in 2D or 3D.
    #  The other parameters are not defined in this example. These parameters are
    #  ID, marker, and elSize.
    #  Since we do not specify an ID the points are automatically assigned IDs,
    #  starting from 0.

    g.point([0, 0])
    g.point([2, 0])
    g.point([2, 1])
    g.point([0, 1])
    g.point([0.5, 0.3])
    g.point([0.3, 0.7])
    g.point([0.7, 0.7])
    g.point([0.8, 0.5])
    g.point([1.7, 0.5])
    g.point([1.5, 0.5])
    g.point([1.7, 0.7])

    g.ellipse([7, 8, 9, 10], marker=50)
    g.spline([0, 1], marker=80)
    g.spline([2, 1])                      # 2
    g.spline([3, 2])                      # 3
    g.spline([0, 3])                      # 4
    g.spline([7, 9], marker=50)           # 5
    g.spline([10, 9])                     # 6
    g.spline([4, 5, 6, 4])                # 7 - This is a closed spline.

    g.surface([4, 3, 2, 1], [[7], [5, 6, 0]])    

    edit_geometry(g)


