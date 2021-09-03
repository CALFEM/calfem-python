# -*- coding: iso-8859-15 -*-
"""
CALFEM Editor Module

Contains functions for implementing a interactive geometry editor.

Written by Karl Eriksson
"""

import os, sys, itertools
from functools import cmp_to_key


import PyQt5

from PyQt5.QtCore import QPoint, QPointF, QLineF, QRectF, QSize, QLocale, Qt

from PyQt5.QtGui import QPen, QColor, QBrush, QPolygon, QPolygonF, QPainter, QFont, QIntValidator, \
    QDoubleValidator, QCloseEvent, QCursor

from PyQt5.QtWidgets import QApplication, QDialog, QWidget, QMainWindow, QFileDialog, QGraphicsView, QGraphicsScene, \
    QGraphicsSceneMouseEvent, QGraphicsItem, QLabel, QGraphicsPolygonItem, QButtonGroup, QToolButton, \
    QGraphicsEllipseItem, QMessageBox, QLineEdit, QFormLayout, QGraphicsLineItem, QTabWidget, QGraphicsTextItem, QFrame, \
    QVBoxLayout, QGridLayout, QScrollArea, QPushButton, QWidgetItem, QGraphicsRectItem, QStackedWidget, QComboBox, \
    QScrollBar, QSpinBox, QTextBrowser

from PyQt5.uic import loadUi

import numpy as np

import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_mpl as cfv
import calfem.utils as cfu
import calfem.editor_resources

setattr(QGraphicsEllipseItem, "marker", None)
setattr(QGraphicsEllipseItem, "localIndex", None)
setattr(QGraphicsLineItem, "localIndex", None)

app = None


class EditorWindow(QMainWindow):
    """MainWindow-klass som hanterar vårt huvudfönster"""

    def __init__(self):
        """Constructor"""
        super(QMainWindow, self).__init__()
        self.app = app

        root = os.path.dirname(os.path.realpath(__file__))
        loadUi(os.path.join(root, 'editor.ui'), self)

        # loadUi('editor.ui', self)   loadUI kan ladd ui-fil och lägga till objekt direkt i klassen.
        self.setWindowTitle("CALFEM Geometry Editor")

        scene = EditorScene(self)
        self.scene = scene
        self.graphicsView.setScene(scene)
        self.graphicsView.setRenderHint(QPainter.Antialiasing)
        self.graphicsView.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.graphicsView.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        self.setMouseTracking(True)

        self.graphicsView.setVisible(True)
        self.return_g = False

        self.loadGeometryButton.clicked.connect(scene.load_from_g)
        self.saveGeomButton.clicked.connect(scene.on_save_geometry_action)
        self.saveMeshButton.clicked.connect(scene.on_save_mesh_action)
        self.saveArraysButton.clicked.connect(scene.on_save_arrays_action)
        self.saveArraysMatlabButton.clicked.connect(scene.on_save_arrays_matlab_action)

        self.gridButtonSurface.pressed.connect(self.toggle_grid)
        self.gridButtonBorder.pressed.connect(self.toggle_grid)
        self.gridSpacingButton.clicked.connect(self.resize_grid_spacing)
        self.gridSnapButtonSurface.clicked.connect(self.scene.toggle_gridsnap)
        self.gridSnapButtonBorder.clicked.connect(self.scene.toggle_gridsnap)
        self.mergeButton.clicked.connect(scene.on_action_merge)
        self.tabWidget.currentChanged.connect(self.on_action_tab)

        #  Disable being able to show mesh or geometry before any geometry has been added
        self.tabWidget.setTabEnabled(2, False)
        self.tabWidget.setTabEnabled(3, False)

        self.arrowButtonSurface.setDown(True)
        self.selectGroupSurface = QButtonGroup(self)
        self.selectGroupSurface.addButton(self.arrowButtonSurface)
        self.selectGroupSurface.addButton(self.rectButton)
        self.selectGroupSurface.addButton(self.polyButton)
        self.selectGroupSurface.addButton(self.addPolyHoleButton)
        self.selectGroupSurface.addButton(self.addRectHoleButton)
        self.selectGroupSurface.addButton(self.panningButtonSurface)

        self.selectGroupBorder = QButtonGroup(self)
        self.selectGroupBorder.addButton(self.arrowButtonBorder)
        self.selectGroupBorder.addButton(self.panningButtonBorder)
        self.selectGroupBorder.addButton(self.setMarkerButton)
        self.selectGroupBorder.addButton(self.splitEdgeButton)

        #self.arrowButton_2 = self.arrowButton

        self.arrowButtonSurface.clicked.connect(self.on_action_arrow)
        self.arrowButtonBorder.clicked.connect(self.on_action_arrow)
        self.panningButtonSurface.clicked.connect(self.on_action_panning)
        self.panningButtonBorder.clicked.connect(self.on_action_panning)
        self.zoomInButtonSurface.clicked.connect(self.zoom_in)
        self.zoomInButtonBorder.clicked.connect(self.zoom_in)
        self.zoomOutButtonSurface.clicked.connect(self.zoom_out)
        self.zoomOutButtonBorder.clicked.connect(self.zoom_out)
        self.rectButton.clicked.connect(self.on_action_rect)
        self.polyButton.clicked.connect(self.on_action_poly)
        self.addPolyHoleButton.clicked.connect(self.on_action_poly_hole)
        self.addRectHoleButton.clicked.connect(self.on_action_rect_hole)
        self.splitEdgeButton.clicked.connect(self.on_action_split)
        self.setMarkerButton.clicked.connect(self.on_action_set_marker)

        # Create a layout behind the drawing canvas to be able to shift to when showing geometry or mesh
        layout = self.verticalLayout
        self.fig = cfv.figure()
        self.canvas = cfv.figure_widget(self.fig)
        layout.addWidget(self.canvas)
        scene.figure_canvas = self.canvas

        scene.scroll_area = self.scrollArea

        self.labelX.setText("")
        self.labelY.setText("")
        self.gridSpacingSpinBox.setValue(self.scene.grid_spacing)
        self.scene.update_labels = self.update_labels
        self.scene.toggle_tab_enabled = self.toggle_tab_enabled
        self.scene.get_graphics_view_size = self.get_graphics_view_size
        self.scene.overlap_warning = self.overlap_warning
        self.scene.set_tooltip = self.set_tooltip
        self.scene.set_view_scale = self.set_view_scale
        self.scene.append_text_browser = self.append_text_browser

        self.radioButtonDrawPointTrue.toggled.connect(self.toggle_draw_points)
        self.radioButtonDisplayPointLabelsTrue.toggled.connect(self.toggle_display_point_labels)
        self.radioButtonDisplayEdgeLabelsTrue.toggled.connect(self.toggle_display_edge_labels)

        self.radioButtonTriangle.toggled.connect(self.toggle_el_type)
        self.refreshMeshButton.clicked.connect(self.scene.show_mesh)
        self.DOFSpinBox.valueChanged.connect(self.update_dofs_per_node)
        self.elSizeSpinBox.valueChanged.connect(self.update_el_size)

        self.clearTextBrowserButton.clicked.connect(self.textBrowser.clear)
        self.hideTextBrowserButton.clicked.connect(self.toggle_text_browser)
        self.showTextBrowserButton.clicked.connect(self.toggle_text_browser)

        self.stackedWidget.setCurrentIndex(1)
        self.tabWidget.setCurrentIndex(0)
        self.tabWidget.setTabText(0, "Surface Mode")
        self.tabWidget.setTabText(1, "Border Mode")
        self.tabWidget.setTabText(2, "Show Geometry")
        self.tabWidget.setTabText(3, "Show Mesh")

        self.overlap_warning_choice = None
        self.grid_snap_on = False
        self.panning_on = False

        self.max_view_scale = 8
        self.min_view_scale = 0.25

        self.graphicsView.wheelEvent = self.override_scroll

        self.text_browser_active = True
        self.showTextBrowserButton.setVisible(False)

        self.show()
        self.raise_()

    def toggle_text_browser(self):
        if self.text_browser_active:
            self.clearTextBrowserButton.setVisible(False)
            self.hideTextBrowserButton.setVisible(False)
            self.showTextBrowserButton.setVisible(True)
            self.textBrowser.setMaximumHeight(4) #.setVisible(False)

            self.text_browser_active = False

        else:
            self.clearTextBrowserButton.setVisible(True)
            self.hideTextBrowserButton.setVisible(True)
            self.showTextBrowserButton.setVisible(False)
            self.textBrowser.setMaximumHeight(200)
            self.text_browser_active = True

    def append_text_browser(self, text):
        self.textBrowser.append(text)

    def update_dofs_per_node(self):
        self.scene.dofs_per_node = self.DOFSpinBox.value()

    def update_el_size(self):
        self.scene.el_size_factor = self.elSizeSpinBox.value()

    def toggle_el_type(self):
        if self.scene.el_type == 2:
            self.scene.el_type = 3
        else:
            self.scene.el_type = 2

    def toggle_draw_points(self):
        self.scene.draw_points = not self.scene.draw_points
        self.scene.show_geom()

    def toggle_display_edge_labels(self):
        self.scene.label_curves = not self.scene.label_curves
        self.scene.show_geom()

    def toggle_display_point_labels(self):
        self.scene.label_points = not self.scene.label_points
        self.scene.show_geom()

    def toggle_grid(self):
        self.scene.toggle_grid()

    def override_scroll(self, event):
        delta = event.angleDelta().y()
        factor = 1
        if delta > 0:
            factor = 2
        elif delta < 0:
            factor = 0.5
        self.scene.set_view_scale(factor)



    def zoom_in(self):
        self.set_view_scale(2)

    def zoom_out(self):
        self.set_view_scale(0.5)

    def set_view_scale(self, factor):

        new_scale = self.graphicsView.transform().m11() * factor

        if self.min_view_scale <= new_scale <= self.max_view_scale:
            self.graphicsView.scale(factor, factor)
        print(self.graphicsView.transform().m11())

    def set_tooltip(self, string):
        self.labelTooltip.setText(string)

    def toggle_grid_snap(self):
        self.grid_snap_on = not self.grid_snap_on
        self.scene.toggle_gridsnap()

    def overlap_warning(self):
        msg = QMessageBox()
        msg.setWindowTitle("Warning!")
        msg.setText("There exists overlapping geometries")
        msg.setIcon(QMessageBox.Warning)
        msg.setStandardButtons(QMessageBox.Cancel | QMessageBox.Ignore)

        msg.setDefaultButton(QMessageBox.Retry)
        msg.setInformativeText("")

        msg.setDetailedText("details")

        msg.buttonClicked.connect(self.popup_button)
        msg.exec_()
        return self.overlap_warning_choice

    def popup_button(self, i):
        self.overlap_warning_choice = i.text()

    def resize_grid_spacing(self):
        self.scene.grid_spacing = self.gridSpacingSpinBox.value()
        if self.scene.grid_on:
            self.scene.toggle_grid("off")
            self.scene.create_grid()
            self.scene.toggle_grid("on")
        else:
            self.scene.create_grid()
            self.scene.toggle_grid("off")

    def get_graphics_view_size(self):
        return self.graphicsView.size().width(), self.graphicsView.size().height()

    def on_action_tab(self):
        self.set_tooltip("")
        index = self.tabWidget.currentIndex()
        self.scene.set_view(index)
        self.on_action_arrow()
        if self.scene.canceled:
            index = 0
            self.scene.canceled = False
            self.tabWidget.setCurrentIndex(0)
        if index == 0:
            self.stackedWidgetRight.setCurrentIndex(1)
            self.stackedWidget.setCurrentIndex(1)
            self.graphicsView.setVisible(True)

            if self.scene.grid_snap:
                self.gridSnapButtonSurface.setChecked(True)
            else:
                self.gridSnapButtonSurface.setChecked(False)

            if self.scene.grid_on:
                self.gridButtonSurface.setChecked(True)
            else:
                self.gridButtonSurface.setChecked(False)

        elif index == 1:
            self.stackedWidget.setCurrentIndex(1)
            self.stackedWidgetRight.setCurrentIndex(1)
            self.graphicsView.setVisible(True)

            if self.scene.grid_snap:
                self.gridSnapButtonBorder.setChecked(True)
            else:
                self.gridSnapButtonBorder.setChecked(False)

            if self.scene.grid_on:
                self.gridButtonBorder.setChecked(True)
            else:
                self.gridButtonBorder.setChecked(False)

        elif index == 2:
            self.stackedWidgetRight.setCurrentIndex(2)
            self.stackedWidget.setCurrentIndex(0)
        elif index == 3:
            self.stackedWidgetRight.setCurrentIndex(0)
            self.stackedWidget.setCurrentIndex(0)

    def update_labels(self, x, y):
        self.labelX.setText(str(x))
        self.labelY.setText(str(-y))

    def toggle_tab_enabled(self, index, boolean):
        self.tabWidget.setTabEnabled(index, boolean)

    def closeEvent(self, a0: QCloseEvent) -> None:
        if self.return_g:
            self.scene.build_gmsh()

    def load_scene_from_g(self, g):
        self.scene.load_from_g(g)

    def on_action_panning(self):
        if not self.panning_on:
            self.toggle_panning()
        for button in self.selectGroupSurface.buttons():
            button.setDown(False)
        for button in self.selectGroupBorder.buttons():
            button.setDown(False)
        self.panningButtonSurface.setDown(True)
        self.panningButtonBorder.setDown(True)

    def toggle_panning(self):
        if self.panning_on:
            self.graphicsView.setDragMode(QGraphicsView.NoDrag)
            self.panning_on = False
        else:
            self.graphicsView.setDragMode(QGraphicsView.ScrollHandDrag)
            self.panning_on = True
            self.scene.mode = "Panning"

    def on_action_poly(self):
        """
        De-selects all other toolButtons and sets the current mode to drawing polygon
        """
        self.graphicsView.setCursor(Qt.CrossCursor)

        self.set_tooltip("Press left mouse button to select starting point of polygon")
        for button in self.selectGroupSurface.buttons():
            button.setDown(False)
        self.polyButton.setDown(True)
        self.scene.mode = "Draw Poly"
        self.scene.hole_mode = False

        if self.panning_on:
            self.toggle_panning()

    def on_action_poly_hole(self):
        """
        De-selects all other toolButtons and sets the current mode to drawing polygon
        """
        self.graphicsView.setCursor(Qt.CrossCursor)
        self.set_tooltip("Press left mouse button to select starting point of polygon hole")
        for button in self.selectGroupSurface.buttons():
            button.setDown(False)
        self.addPolyHoleButton.setDown(True)
        self.scene.mode = "Draw Poly"
        self.scene.hole_mode = True

        if self.panning_on:
            self.toggle_panning()

    def on_action_rect(self):
        """
        De-selects all other toolButtons and sets the current mode to drawing rectangle
        """
        self.graphicsView.setCursor(Qt.CrossCursor)
        self.set_tooltip("Press left mouse button to select starting point of rectangle")
        for button in self.selectGroupSurface.buttons():
            button.setDown(False)
        self.rectButton.setDown(True)
        self.scene.mode = "Draw Rect"
        self.scene.hole_mode = False

        if self.panning_on:
            self.toggle_panning()

    def on_action_rect_hole(self):
        """
        De-selects all other toolButtons and sets the current mode to drawing polygon
        """
        self.graphicsView.setCursor(Qt.CrossCursor)
        self.set_tooltip("Press left mouse button to select starting point of rectangular hole")
        for button in self.selectGroupSurface.buttons():
            button.setDown(False)
        self.addRectHoleButton.setDown(True)
        self.scene.mode = "Draw Rect"
        self.scene.hole_mode = True

        if self.panning_on:
            self.toggle_panning()

    def on_action_split(self):
        self.set_tooltip("Click on existing edge to split by adding an extra node")
        self.scene.mode = "Split Line"
        self.graphicsView.setCursor(Qt.ArrowCursor)

        for button in self.selectGroupBorder.buttons():
            button.setDown(False)
        self.splitEdgeButton.setDown(True)

        if self.panning_on:
            self.toggle_panning()

    def on_action_set_marker(self):
        self.set_tooltip("Select edge to add marker")
        self.scene.mode = "Set Marker"
        self.graphicsView.setCursor(Qt.PointingHandCursor)

        for button in self.selectGroupBorder.buttons():
            button.setDown(False)
        self.setMarkerButton.setDown(True)

        if self.panning_on:
            self.toggle_panning()

    def on_action_arrow(self):
        """
        De-selects all other toolButtons and sets the current mode Arrow (normal mode, mouse acting as pointer)
        """
        self.set_tooltip("")
        self.graphicsView.setCursor(Qt.ArrowCursor)
        for button in self.selectGroupSurface.buttons():
            button.setDown(False)
        for button in self.selectGroupBorder.buttons():
            button.setDown(False)
        self.arrowButtonSurface.setDown(True)
        self.arrowButtonBorder.setDown(True)
        self.scene.mode = "Arrow"

        if self.panning_on:
            self.toggle_panning()


class EditorScene(QGraphicsScene, QMainWindow):

    def __init__(self, parent=None):
        QGraphicsScene.__init__(self, parent)
        self.x_prev = None
        self.y_prev = None
        self.g = cfg.Geometry()
        self.g_index = 0

        self.LUBlue = QColor(0, 0, 128)
        self.LUBronze = QColor(156, 97, 20)

        self.poly_list = []
        self.hole_list = []
        self.edge_list = []

        self.point_coord_list = np.zeros((1, 2))
        self.point_marker_list = []
        self.line_marker_list = []

        self.grid_on = True
        self.grid_built = False
        self.grid = []
        self.grid_snap = False
        self.grid_snap_last_x = 0
        self.grid_snap_last_y = 0
        self.grid_spacing = 20
        self.grid_max_scale = 10
        self.grid_min_scale = 0.2
        self.create_grid()
        self.toggle_grid(True)
        self.toggle_gridsnap()

        self.first_draw = True
        self.prev_point = None
        self.connecting_rect = None
        self.connecting_line = None
        self.connecting_line_list = []
        self.drawing_poly = QPolygonF()
        self.drawing_points = []
        self.drawing_rect = QPolygonF()
        self.hole_mode = False

        self.g = None
        self.mesh = None

        self.line_marker_dict = {}
        self.line_marker_index = 1

        self.numbers_shown = False
        self.line_number_list = []

        self.potential_edge_splitters = []

        # Standard Mesh settings
        self.el_type = 2
        self.dofs_per_node = 1
        self.el_size_factor = 25

        # Standard Geom view settings
        self.draw_points = True
        self.label_points = True
        self.label_curves = True

        self.mode = "Arrow"
        self.view = "Surface View"

        self.update_labels = None
        self.toggle_tab_enabled = None
        self.get_graphics_view_size = None
        self.overlap_warning = None
        self.set_tooltip = None
        self.figure_canvas = None
        self.scroll_area = None
        self.append_text_browser = None

        self.node_splitter = self.addEllipse(-3, -3, 6, 6)
        self.node_splitter.setVisible(False)
        self.split_edge = None

        self.canceled = False

        self.prev_selected_point = None

    def set_view(self, index):
        """
        Cycle between the diffrent tabs in the tabwidget and execute the corresponding action of the view
        """
        view = index
        self.reset_scroll_area()
        if view == 0:
            self.toggle_surface_mode()
            if self.numbers_shown:
                self.toggle_display_numbers()
            self.view = "Surface View"
        elif view == 1:
            self.toggle_border_mode()
            self.view = "Border View"
        elif view == 2:
            if self.show_geom() == "Canceled":
                self.canceled = True
            else:
                self.view = "Geom View"
        elif view == 3:
            self.show_mesh()
            self.mesh_scroll_area()
            self.view = "Mesh View"

    def move_corner(self, circ, poly, new_x, new_y):
        """
        Update the corner of a polygon when dragging a corner circle of the polygon
        """
        poly_list = self.poly_to_list(poly, "Local")  # Extract the position of the polygon points before moving the point

        index = poly_list.index(circ.pos())  # Get the selected circles index in the polygon

        circ.setPos(new_x - poly.scenePos().x(), new_y - poly.scenePos().y())  # Move the selected circle

        poly_list[index] = circ.pos()  # Update the coords of the point in the polygon list
        poly.setPolygon(QPolygonF(poly_list))  # Update the polygon with the new list

        # Loop through all the edges of the polygon to determine which two lines are connected to the moved point,
        # update these edges with the new point to match the movement,also move any edge labeltext connected to the edge
        # based on that we know that the connected lines are one with the same index as the corner point and one with
        # index-1
        # exception case 1:  when selecting corner with index 0, then we use that the last edge is indexed with a
        # negative sign
        # exception case 2: when selecting corner with highest index last edge is indexed with a negative sign,
        # catch this with a special if statement for the highest index
        for item in poly.childItems():
            if isinstance(item, PyQt5.QtWidgets.QGraphicsLineItem):
                if circ.localIndex == 0:
                    if item.localIndex < 0:
                        line = item.line()
                        line.setP2(circ.pos())
                        item.setLine(line)
                        if item.childItems():
                            text = item.childItems()[0]
                            text.setPos((item.line().x1() + item.line().x2()) / 2,
                                        (item.line().y1() + item.line().y2()) / 2)
                if item.localIndex == circ.localIndex:
                    line = item.line()
                    line.setP2(circ.pos())
                    item.setLine(line)
                    if item.childItems():
                        text = item.childItems()[0]
                        text.setPos((item.line().x1() + item.line().x2()) / 2,
                                    (item.line().y1() + item.line().y2()) / 2)
                if item.localIndex == circ.localIndex + 1:
                    line = item.line()
                    line.setP1(circ.pos())
                    item.setLine(line)
                    if item.childItems():
                        text = item.childItems()[0]
                        text.setPos((item.line().x1() + item.line().x2()) / 2,
                                    (item.line().y1() + item.line().y2()) / 2)
                if circ.localIndex == poly.polygon().size() - 1:
                    if item.localIndex < 0:
                        line = item.line()
                        line.setP1(circ.pos())
                        item.setLine(line)
                        if item.childItems():
                            text = item.childItems()[0]
                            text.setPos((item.line().x1() + item.line().x2()) / 2,
                                        (item.line().y1() + item.line().y2()) / 2)

    def keyPressEvent(self, event):
        # If delete is pressed and a polygon is selected: remove the selected polygon
        if event.key() == 16777223:
            if self.selectedItems():
                if isinstance(self.selectedItems()[0], PyQt5.QtWidgets.QGraphicsPolygonItem):
                    self.delete_polygon(self.selectedItems()[0])

    def mouseMoveEvent(self, event):
        # If polygons have been added enable the geometry and mesh options in the tabview, else keep disabled
        if self.toggle_tab_enabled is not None:
            if not self.poly_list:
                self.toggle_tab_enabled(2, False)
                self.toggle_tab_enabled(3, False)
            else:
                self.toggle_tab_enabled(2, True)
                self.toggle_tab_enabled(3, True)

        # When moving the mouse in the graphicsScene display coords in label
        x = event.scenePos().x()
        y = event.scenePos().y()

        if self.update_labels is not None:
            self.update_labels(x, y)
        if self.grid_snap:
            x = round(x / self.grid_spacing) * self.grid_spacing
            y = round(y / self.grid_spacing) * self.grid_spacing

        if self.mode == "Split Line":

            x_closest = None
            y_closest = None
            edge_closest = None

            print([event.scenePos().x(),event.scenePos().y()] in self.point_coord_list.tolist())

            # Check a area around the mouse point to search for any edges to snap to
            edge_point_list = []
            # Check a square area with width 10 if there is any edge that contains the point, store all edges
            # that contains a point
            # Most inefficient part of the code, noticeable lag when creating many edges in the canvas
            for i, j in itertools.product(range(-10, 10), range(-10, 10)):
                for edge in self.edge_list:
                    p = QPointF(0, 0)
                    p.setX(x + i - edge.scenePos().x())
                    p.setY(y + j - edge.scenePos().y())
                    if[x + i - edge.scenePos().x(), y + j - edge.scenePos().y()] in self.point_coord_list.tolist():
                        pass
                    elif edge.contains(p):
                        edge_point_list.append([x + i, y + j, edge])

            smallest = np.inf
            edge_point_list = np.array(edge_point_list)
            # Loop through all potential points, if they exist, and choose the one closest to the mouse pointer
            # as the point to snap to
            for row in edge_point_list:
                coords = np.array([row[0], row[1]])
                dist = np.linalg.norm(coords - np.array([event.scenePos().x(), event.scenePos().y()]))
                if dist < smallest:
                    smallest = dist
                    x_closest = coords[0]
                    y_closest = coords[1]
                    edge_closest = row[2]

            if x_closest:
                self.node_splitter.setPos(x_closest, y_closest)
                self.node_splitter.setVisible(True)
                self.split_edge = edge_closest
            else:
                self.node_splitter.setVisible(False)
                self.split_edge = None

        # If the mode is set to Arrow, mouse behaves like a normal pointer
        if self.mode == "Arrow":
            if self.selectedItems():
                # If a polygon is selected update the polygons position with the corresponding mouse movement
                if isinstance(self.selectedItems()[0], PyQt5.QtWidgets.QGraphicsPolygonItem):
                    if self.grid_snap:
                        self.selectedItems()[0].moveBy(x - self.grid_snap_last_x, y - self.grid_snap_last_y)
                        self.grid_snap_last_x = x
                        self.grid_snap_last_y = y
                    else:
                        self.selectedItems()[0].moveBy(x - event.lastScenePos().x(), y - event.lastScenePos().y())
                # If a circle is selected update the circles position with the corresponding mouse movement and
                # update the parent polygon with the changed corner
                if isinstance(self.selectedItems()[0], PyQt5.QtWidgets.QGraphicsEllipseItem):

                    circ = self.selectedItems()[0]
                    poly = circ.parentItem()

                    # Check a area around the mouse point to search for any edges to snap to
                    edge_point_list = []
                    temp_edge_list = []
                    for edge in self.edge_list:
                        if edge in poly.childItems():
                            pass  # If the edge is in the parent polygon pass to avoid snapping to self
                        else:
                            temp_edge_list.append(edge)
                    # Check a square area with width 10 if there is any edge that contains the point, store all edges
                    # that contains a point
                    # Most inefficient part of the code, noticeable lag when creating many edges in the canvas
                    for i, j in itertools.product(range(-10, 10), range(-10, 10)):
                        for edge in temp_edge_list:
                            p = QPointF(0, 0)
                            p.setX(x+i - edge.scenePos().x())
                            p.setY(y+j - edge.scenePos().y())
                            if edge.contains(p):
                                edge_point_list.append([x + i, y + j])

                    smallest = np.inf
                    edge_point_list = np.array(edge_point_list)
                    # Loop through all potential points, if they exist, and choose the one closest to the mouse pointer
                    # as the point to snap to
                    for coords in edge_point_list:
                        dist = np.linalg.norm(coords - np.array([event.scenePos().x(), event.scenePos().y()]))
                        if dist < smallest:
                            smallest = dist
                            x = coords[0]
                            y = coords[1]
                            # All points that are at some point in time snapped are added to the potentialEdgeSplitters
                            # to avoid having to loop through all points in later stages
                            if circ not in self.potential_edge_splitters:
                                self.potential_edge_splitters.append(circ)

                    # After check if there are any points to snap to, priority to snap to points over edges
                    # Add a templist and remove own points to avoid snapping with self
                    templist = self.point_coord_list
                    for point in self.poly_to_list(poly, "Global"):
                        if point == circ.scenePos():
                            pass  # This point has already been removed, catch to avoid error in deletion
                        else:
                            templist = np.delete(templist, np.where(np.all(templist == [[point.x(), point.y()]], axis=1))[0][0], axis=0)
                    # Check if any point in the global point list is within snapping threshold, if so snap to that point
                    if (np.linalg.norm(templist - [x, y], axis=1) < 10).any():
                        coords = templist[np.where((np.linalg.norm(templist - [x, y], axis=1) < 10))]
                        x = coords[0][0]
                        y = coords[0][1]

                    # Move corner of the polygon to the new x and y, if no snapping has occurred it is the mouse coords
                    self.move_corner(circ, poly, x, y)

                    #  För att uppdatera nummer när noder flyttas, fungerar inte helt hundra än
                    #if self.tabWidget.currentIndex() == 1:
                    #    if self.numbers_shown:
                    #        self.show_numbers()
                    #        self.show_numbers()

        if self.mode == "Draw Poly":
            # This is to display the line in the polygon from the previous point to see where the new line would occur
            if self.first_draw:
                pass  # Don't draw if the first point has not been initiated
            else:
                # else if there is an existing line update that one with new x,y, if no line create a new one with
                # end point at x,y
                if self.connecting_line:
                    self.connecting_line.setLine(QLineF(self.prev_point, QPointF(x, y)))
                else:
                    self.connecting_line = self.addLine(QLineF(self.prev_point, QPointF(x, y)))

        if self.mode == "Draw Rect":
            # This is to display the rectangle from the previous point to see where the new rectangle would occur
            if self.first_draw:
                pass  # Don't draw if the first point has not been initiated
            else:
                # else if there is an existing rectangle update that one with new x,y, if no rectangle create a new one
                if self.connecting_rect:
                    if self.prev_point.x() > x and self.prev_point.y() > y:
                        self.connecting_rect.setRect(QRectF(QPointF(x, y), self.prev_point))
                    elif self.prev_point.x() > x:
                        self.connecting_rect.setRect(QRectF(QPointF(x, self.prev_point.y()), QPointF(self.prev_point.x(),y)))
                    elif self.prev_point.y() > y:
                        self.connecting_rect.setRect(QRectF(QPointF(self.prev_point.x(),y), QPointF(x, self.prev_point.y())))
                    else:
                        self.connecting_rect.setRect(QRectF(self.prev_point, QPointF(x, y)))
                else:
                    self.connecting_rect = self.addRect(QRectF(self.prev_point, QPointF(x, y)))

    def delete_polygon(self, poly : QGraphicsPolygonItem, delete_from_coord_list=False):
        self.poly_list.remove(poly)
        if poly in self.hole_list:
            self.hole_list.remove(poly)
        for item in poly.childItems():
            if isinstance(item, PyQt5.QtWidgets.QGraphicsLineItem):
                self.edge_list.remove(item)
        if delete_from_coord_list:
            for point in self.poly_to_list(poly, "Global"):
                self.point_coord_list = np.delete(self.point_coord_list, np.where(
                    np.all(self.point_coord_list == [[point.x(), point.y()]], axis=1))[0][0], axis=0)
        poly.hide()

    def reset_scroll_area(self):
        # Clear the scroll area
        scroll_area_widget_contents = QWidget()
        if self.scroll_area:
            self.scroll_area.setWidget(scroll_area_widget_contents)

    def mouseDoubleClickEvent(self, event):
        if self.mode == "Arrow":
            super(EditorScene, self).mouseDoubleClickEvent(event)

            # If in the surface view highlight the polygon to allow updating exact values of the corner points
            if self.view == "Surface View":
                if self.selectedItems():
                    if isinstance(self.selectedItems()[0], PyQt5.QtWidgets.QGraphicsPolygonItem):
                        index = 0
                        poly = self.selectedItems()[0]
                        scroll_area_widget_contents = QWidget()
                        grid = QGridLayout(scroll_area_widget_contents)
                        if self.scroll_area:
                            self.scroll_area.setWidget(scroll_area_widget_contents)

                        # Add a x- and y- editor for each point of the polygon
                        for point in self.poly_to_list(poly, "Global"):
                            label_x = QLineEdit(str(point.x()))
                            label_x.setValidator(QDoubleValidator())
                            label_y = QLineEdit(str(-point.y()))
                            label_y.setValidator(QDoubleValidator())
                            grid.addWidget(label_x, index, 0)
                            grid.addWidget(label_y, index, 1)
                            index += 1

                        def update():
                            # Update the polygon with the new edited values
                            i = 0
                            for child_item in poly.childItems():
                                if isinstance(child_item, PyQt5.QtWidgets.QGraphicsEllipseItem):
                                    if child_item.localIndex == i:
                                        x = float(grid.itemAtPosition(i, 0).widget().text())
                                        y = -float(grid.itemAtPosition(i, 1).widget().text())
                                        circ = child_item
                                        self.move_corner(circ, poly, x, y)
                                        point = circ.scenePos()
                                        self.point_coord_list = np.append(self.point_coord_list, [[point.x(), point.y()]], axis=0)
                                        i += 1

                        update_button = (QPushButton("Update"))
                        grid.addWidget(update_button, index + 1, 1)
                        update_button.clicked.connect(update)

    def add_marker(self, item, marker_text):
        item.__setattr__("marker", self.line_marker_index)
        self.line_marker_dict[self.line_marker_index] = marker_text
        self.line_marker_index += 1

        if isinstance(item, PyQt5.QtWidgets.QGraphicsEllipseItem):
            self.point_marker_list.append(item)
            item.setBrush(QColor("Red"))

            # If there is an old marker remove it and remove as child item
            if item.childItems():
                item.childItems()[0].setVisible(False)
                item.childItems()[0].setParentItem(None)

            text = self.addText(str(marker_text), QFont("Helvetica", 8, QFont.Bold))
            text.setParentItem(item)
            text.setPos(-5, 5)

        if isinstance(item, PyQt5.QtWidgets.QGraphicsLineItem):
            self.line_marker_list.append(item)
            item.setPen(QColor("Red"))

            # If there is an old marker remove it and remove as child item
            if item.childItems():
                item.childItems()[0].setVisible(False)
                item.childItems()[0].setParentItem(None)

            text = self.addText(str(marker_text), QFont("Helvetica", 8, QFont.Bold))
            text.setParentItem(item)
            text.setPos((item.line().x1() + item.line().x2()) / 2 - 15, (item.line().y1() + item.line().y2()) / 2 - 15)

    def mousePressEvent(self, event):
        x = event.scenePos().x()
        y = event.scenePos().y()

        if self.mode == "Split Line":
            if self.split_edge:
                edge = self.split_edge
                poly = edge.parentItem()
                poly_list = self.poly_to_list(poly, "Global")
                line = edge.line()
                poly_is_hole = False

                if poly in self.hole_list:
                    poly_is_hole = True

                line.translate(edge.scenePos())
                p1_index = poly_list.index(line.p1())
                p2_index = poly_list.index(line.p2())
                if abs(p1_index - p2_index) > 1:  # If difference is larger than one means it is the first and last point
                    if p1_index > p2_index:
                        insert_index = p1_index
                    else:
                        insert_index = p2_index
                else:
                    if p1_index < p2_index:
                        insert_index = p1_index
                    else:
                        insert_index = p2_index
                insert_index += 1
                poly_list.insert(insert_index, self.node_splitter.scenePos())
                new_poly = QPolygonF()
                for p in poly_list:
                    new_poly << p

                self.add_poly_to_scene(new_poly, hole_mode=poly_is_hole)
                self.delete_polygon(poly, True)

                self.toggle_border_mode()  # To toggle the newly created polygon to border mode

        if self.mode == "Arrow":
            if event.button() != 1:
                # Return if button clicked is any is any other than left mouse
                return
            super(EditorScene, self).mousePressEvent(event)

            # If a polygon or circle is clicked remove the coordinates from the global coordinate lists to avoid
            # snapping to self etc. point is re-added on mouse release and hence also updating the coord list
            if self.selectedItems():
                if isinstance(self.selectedItems()[0], PyQt5.QtWidgets.QGraphicsPolygonItem):
                    for point in self.poly_to_list(self.selectedItems()[0], "Global"):
                        self.point_coord_list = np.delete(self.point_coord_list, np.where(np.all(self.point_coord_list == [[point.x(), point.y()]], axis=1))[0][0], axis=0)
                        if self.grid_snap:
                            x = round(x / self.grid_spacing) * self.grid_spacing
                            y = round(y / self.grid_spacing) * self.grid_spacing
                            self.grid_snap_last_x = x
                            self.grid_snap_last_y = y
                if isinstance(self.selectedItems()[0], PyQt5.QtWidgets.QGraphicsEllipseItem):
                    self.prev_selected_point = self.selectedItems()[0]
                    point = self.selectedItems()[0].scenePos()
                    self.point_coord_list = np.delete(self.point_coord_list, np.where(np.all(self.point_coord_list == [[point.x(), point.y()]], axis=1))[0][0], axis=0)

        if self.mode == "Draw Poly":
            # The drawing of the polygon is done by adding a point wherever the user clicks in the canvas, support
            # points and lines are being drawn to show the users where the points and lines will be added. The
            # process is finished by right-clicking to close the polygon
            if self.grid_snap:
                x = round(x / self.grid_spacing) * self.grid_spacing
                y = round(y / self.grid_spacing) * self.grid_spacing

            if event.button() == 2:
                # If a polygon is being drawn, finish the polygon by clicking right mouse button. This will close the
                # polygon and remove the lines drawn as support to show the polygon and replace them with the actual
                # edges and points of the polygon
                if self.first_draw:
                    pass
                else:
                    self.add_poly_to_scene(self.drawing_poly, hole_mode=self.hole_mode)
                    self.drawing_poly = QPolygonF()

                    for p in self.drawing_points:
                        p.setVisible(False)

                    for line in self.connecting_line_list:
                        line.setVisible(False)

                    self.connecting_line.setVisible(False)
                    self.connecting_line = None
                    self.first_draw = True
                    if self.set_tooltip:
                        self.set_tooltip("")
            elif event.button() == 1:
                self.append_text_browser("Polygon point added at: (" + str(x) + " , " + str(y) + " )")
                if self.first_draw:
                    self.drawing_poly << QPointF(x, y)
                    point = self.addEllipse(x - 3, y - 3, 6, 6)
                    self.prev_point = QPointF(x, y)
                    self.drawing_points.append(point)
                    self.first_draw = False
                    if self.set_tooltip:
                        self.set_tooltip("Left-click to add new point or right-click to close polygon")
                else:
                    self.drawing_poly << QPointF(x, y)
                    point = self.addEllipse(x-3, y - 3, 6, 6)
                    line = self.addLine(QLineF(self.prev_point, QPointF(x, y)))
                    self.connecting_line_list.append(line)
                    self.prev_point = QPointF(x, y)
                    self.drawing_points.append(point)

        if self.mode == "Draw Rect":
            if self.grid_snap:
                x = round(x / self.grid_spacing) * self.grid_spacing
                y = round(y / self.grid_spacing) * self.grid_spacing

            if self.first_draw:
                self.prev_point = QPointF(x, y)
                self.first_draw = False
            elif self.prev_point == QPointF(x, y):
                pass  # Catch to avoid crash when pressing startpoint twice
            else:
                r = self.connecting_rect.rect()
                self.drawing_rect << QPointF(r.x(), r.y())
                self.drawing_rect << QPointF(r.x() + r.width(), r.y())
                self.drawing_rect << QPointF(r.x() + r.width(), r.y() + r.height())
                self.drawing_rect << QPointF(r.x(), r.y() + r.height())

                self.add_poly_to_scene(self.drawing_rect, hole_mode=self.hole_mode)
                self.drawing_rect = QPolygonF()

                self.connecting_rect.setVisible(False)
                self.connecting_rect = None
                self.first_draw = True

        # If in the edge view allow selection of corner points and edges to add markers
        if self.mode == "Set Marker":
            super(EditorScene, self).mousePressEvent(event)
            if self.selectedItems():
                if isinstance(self.selectedItems()[0], PyQt5.QtWidgets.QGraphicsEllipseItem) or isinstance(
                        self.selectedItems()[0], PyQt5.QtWidgets.QGraphicsLineItem):
                    item = self.selectedItems()[0]

                    e1 = QLineEdit()
                    e1.setFont(QFont("Arial", 10))

                    # If update is pressed and there is a text add the text as a marker to the double clicked item
                    def pressed():
                        if e1.text() != "":
                            self.add_marker(item, e1.text())

                    update_button = QToolButton()
                    update_button.setText("Update")
                    update_button.clicked.connect(pressed)

                    scroll_area_widget_contents = QWidget()
                    flo = QFormLayout(scroll_area_widget_contents)
                    if self.scroll_area:
                        self.scroll_area.setWidget(scroll_area_widget_contents)
                    flo.addRow("Set Marker", e1)
                    flo.addRow("", update_button)

    def mouseReleaseEvent(self, event):
        super(EditorScene, self).mouseReleaseEvent(event)
        # If a point or polygon is selected releasing the mouse will de-select the object and add the
        # current coordinates back to the global coordinate list to update to the new position
        if self.mode == "Arrow":
            if self.selectedItems():
                if isinstance(self.selectedItems()[0], PyQt5.QtWidgets.QGraphicsPolygonItem):
                    for point in self.poly_to_list(self.selectedItems()[0], "Global"):
                        self.point_coord_list = np.append(self.point_coord_list, [[point.x(), point.y()]], axis=0)
                if isinstance(self.selectedItems()[0], PyQt5.QtWidgets.QGraphicsEllipseItem):
                    point = self.selectedItems()[0].scenePos()
                    self.point_coord_list = np.append(self.point_coord_list, [[point.x(), point.y()]], axis=0)
            self.clearSelection()

    def toggle_border_mode(self):
        # Make all polygon surfaces hidden
        for poly in self.poly_list:
            poly.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable, False)
            poly.setBrush(QColor(0, 0, 0, 0))

        # Enable selection of the edges of the polygon, if the edge has a marker display it
        for edge in self.edge_list:
            edge.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable, True)
            if edge.childItems():
                text = edge.childItems()[0]
                text.setVisible(True)

    def toggle_surface_mode(self):
        # Make all polygon surfaces visible
        for poly in self.poly_list:
            poly.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable, True)
            if poly in self.hole_list:
                poly.setBrush(QBrush(QColor(255, 255, 255)))
            else:
                poly.setBrush(QBrush(QColor(0, 0, 0, 50)))

        # Disable the selection of edges and hide the marker if there is one
        for edge in self.edge_list:
            edge.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable, False)
            if edge.childItems():
                text = edge.childItems()[0]
                text.setVisible(False)

    def add_poly_to_scene(self, polygon, point_marker_dict=None, curve_marker_dict=None, hole_mode=False):
        if hole_mode:
            poly = self.addPolygon(polygon, QPen(QColor(0, 0, 0, 0)), QBrush(QColor(255, 255, 255)))
            poly.setZValue(1)
            self.poly_list.append(poly)
            self.hole_list.append(poly)
        else:
            poly = self.addPolygon(polygon, QPen(QColor(0, 0, 0, 0)), QBrush(QColor(0, 0, 0, 50)))
            self.poly_list.append(poly)
        self.add_poly_corners(poly, point_marker_dict)
        self.add_polyedges(poly, curve_marker_dict)
        poly.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        poly.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        return poly

    def add_polyedges(self, poly_item, marker_dict=None):
        """
        Add the edges of a polygon,
        """

        poly = poly_item.polygon()

        for i in range(1, poly.size()+1):
            if i == poly.size():
                p1 = poly.at(i-1)
                p2 = poly.at(0)
                index = -poly.size()

            else:
                p1 = poly.at(i - 1)
                p2 = poly.at(i)
                index = i

            line = self.addLine(QLineF(p1, p2), QPen(self.LUBronze, 2))
            line.__setattr__("localIndex", index)
            line.setParentItem(poly_item)
            self.edge_list.append(line)

            # Used to pass markers when loading a g
            if marker_dict:
                if marker_dict[i] != 0:
                    self.add_marker(line, marker_dict[i])

    def add_poly_corners(self, poly_item, marker_dict=None):
        """
        Add circle-shaped ellipse items at every defined point in the polygon, circles are added as children of the
        polygon and defined in local coordinates of the polygon (not changed when moving in the scene)
        """
        poly = poly_item.polygon()

        for i in range(poly.size()):
            point = poly.at(i)
            p = self.addEllipse(-3, -3,  6, 6, self.LUBronze, self.LUBronze)
            p.setZValue(2)  # Make sure corners always in front of polygon surfaces
            p.setParentItem(poly_item)
            p.__setattr__("localIndex", int(i))
            p.setPos(point.x(), point.y())
            p.setFlag(QGraphicsItem.ItemIsSelectable)
            p.setFlag(QGraphicsItem.ItemIsMovable)
            self.point_coord_list = np.append(self.point_coord_list, [[p.x(), p.y()]], axis=0)

            self.potential_edge_splitters.append(p)

            # Used to pass markers when loading a g
            if marker_dict:
                if marker_dict[i] != 0:
                    self.add_marker(p, marker_dict[i])

    def poly_to_list(self, poly, scope : str):
        """Extract the points from a QGraphicsPolygonItem or a QPolygonF and return a list of all the containing QPointF
        , scope to be chosen as Global return scene coordinates otherwise returns local coordinates """
        inner_list = []
        x = 0
        y = 0

        # To be able to handle input as both QGraphicsPolygonItem and QPolygonF
        if isinstance(poly, PyQt5.QtWidgets.QGraphicsPolygonItem):
            if scope == "Global":
                x = poly.x()
                y = poly.y()
            poly = poly.polygon()

        for i in range(poly.size()):
            inner_list.append(QPointF(poly.at(i).x() + x, poly.at(i).y() + y))
        return inner_list

    def load_from_g(self, g):
        # If no g is passed open the file dialog to let the user select a .cfg file to load
        if not g:
            dlg = QFileDialog()
            dlg.setNameFilter("(*.cfg)")
            dlg.setFileMode(QFileDialog.ExistingFile)

            if dlg.exec_():
                filenames = dlg.selectedFiles()
                name = filenames[0]
                g = cfu.load_geometry(name)

            if not g:
                return  # Abort if no file is loaded

        # For every surface in the g object loop through the points on the surface and make use of the fact that
        # g.stuffOnSurfaces returns points the same order as they should be added to the polygon
        for surfaceID in g.surfaces:
            poly = QPolygonF()
            i = 0
            point_marker_dict = {}
            for pointID in g.stuffOnSurfaces(surfaceID)[0]:
                coords = g.getPointCoords(pointID)
                poly << QPointF(coords[0], -coords[1])  # canvas has flipped y-axis
                point_marker_dict[i] = g.points[pointID][2]
                i += 1

            curve_marker_dict = {}
            i = 1
            for curveID in g.stuffOnSurfaces(surfaceID)[1]:
                curve_marker_dict[i] = g.curves[curveID][2]
                i += 1

            self.add_poly_to_scene(poly, point_marker_dict, curve_marker_dict)

    def toggle_display_numbers(self):
        """
        Toggles the display of numbers of edges and nodes while in the border view
        """
        if self.numbers_shown:
            for text in self.line_number_list:
                text.setVisible(False)
            self.numbers_shown = False
        else:
            point_index = 0
            line_index = 0
            added_points = []
            added_lines = []
            for poly in self.poly_list:
                first = True
                first_index = point_index
                prev_index = point_index
                poly_list = self.on_action_split(poly)

                for point in poly_list:  # loop over all points in the polygon
                    if [point.x(), -point.y()] in added_points:  # if the point already exists don't add again
                        if first:  # If it is the first point of the poly set existing point as first_index
                            first = False
                            current_index = added_points.index([point.x(), -point.y()])
                            first_index = current_index
                        else:  # Check if the line from previous point already exists, in that case don't att a new one
                            current_index = added_points.index([point.x(), -point.y()])
                            if added_lines.__contains__([prev_index, current_index]):
                                pass
                            elif added_lines.__contains__([current_index, prev_index]):
                                pass
                            else:  # else add a line from the previous point to the current pre-existing one
                                p1 = added_points[prev_index]
                                p2 = added_points[current_index]
                                text = self.addText(str(line_index), QFont("Helvetica",8, QFont.Bold))
                                text.setPos((p1[0] + p2[0]) / 2, (-p1[1] - p2[1]) / 2)
                                line_index += 1
                                added_lines.append([prev_index, current_index])

                    else:
                        current_index = point_index
                        text = self.addText(str(current_index), QFont("Helvetica"))
                        text.setPos(point.x(), point.y())
                        self.line_number_list.append(text)
                        point_index += 1
                        added_points.append([point.x(), -point.y()])
                        if first:
                            first = False
                        else:
                            p1 = added_points[prev_index]
                            p2 = added_points[current_index]
                            text = self.addText(str(line_index), QFont("Helvetica",8, QFont.Bold))
                            text.setPos((p1[0] + p2[0]) / 2, (-p1[1] - p2[1]) / 2)
                            circle = self.addEllipse(1,3, 15,15)
                            circle.setParentItem(text)
                            self.line_number_list.append(text)
                            line_index += 1
                            added_lines.append([prev_index, current_index])
                    prev_index = current_index

                if added_lines.__contains__([current_index, first_index]):
                    pass
                elif added_lines.__contains__([first_index, current_index]):
                    pass
                else:
                    p1 = added_points[first_index]
                    p2 = added_points[current_index]
                    text = self.addText(str(line_index), QFont("Helvetica",8, QFont.Bold))
                    text.setPos((p1[0] + p2[0]) / 2, (-p1[1] - p2[1]) / 2)
                    self.line_number_list.append(text)
                    line_index += 1
                    added_lines.append([current_index, first_index])
                self.numbers_shown = True

    def on_action_split(self, polygon):
        """
            Creates a polygon list but includes checking of overlapping points on edges, in case of overlap adds a extra
            point in the list where the overlap occurs.
        """
        added = 0
        polygon_item = polygon.polygon()
        polygon_item.translate(polygon.x(), polygon.y())

        # Comparator to determine which x value of two points is the highest
        def compare_x(item1, item2):
            if item1.x() < item2.x():
                return -1
            elif item1.x() > item2.x():
                return 1
            else:
                return 0

        # Comparator to determine which y value of two points is the highest
        def compare_y(item1, item2):
            if item1.y() < item2.y():
                return -1
            elif item1.y() > item2.y():
                return 1
            else:
                return 0

        # Create two lists, one sorted by ascending x-values, one by ascending y-values
        x_list = sorted(self.potential_edge_splitters, key=cmp_to_key(compare_x))
        y_list = sorted(self.potential_edge_splitters, key=cmp_to_key(compare_y))

        # Loop over all children to the polygon
        for item in polygon.childItems():
            # Look only at edges (overlapping of points is handled elsewhere)
            if isinstance(item, PyQt5.QtWidgets.QGraphicsLineItem):
                edge = item

                p1 = edge.line().p1()
                p2 = edge.line().p2()
                added_this = 0

                # Choose the direction with the largest disparity (to avoid scenario of straight lines)
                # then use the sorted list for that direction
                if abs(p1.x() - p2.x()) > abs(p1.y() - p2.y()):
                    mode = "X"
                    circ_list = x_list
                else:
                    mode = "Y"
                    circ_list = y_list

                for circ in circ_list:
                    poly = circ.parentItem()
                    p = circ.scenePos()

                    # temp_p needed since edge.contains does not account for the edge being moved in the canvas
                    temp_p = circ.scenePos()
                    temp_p.setX(temp_p.x() - edge.scenePos().x())
                    temp_p.setY(temp_p.y() - edge.scenePos().y())

                    # Find the edges to split which contain temp_p, if the edge contains decide the orientation (in x-
                    # or y-direction decided earlier) of p1 and p2, based on this insert the new point in the polygon
                    # in the correct position
                    if edge.contains(temp_p):
                        if edge in poly.childItems():
                            pass  # Ignore if the edge is in the same polygon as the point
                        else:
                            if temp_p == p1 or temp_p == p2:
                                pass  # Don't compare if it contains an edge point, instead handled later by the overlapping points
                            elif mode == "Y":
                                if p1.y() < p2.y():  # Left to right
                                    index = abs(edge.localIndex)
                                    polygon_item.insert(index + added, p)
                                    added += 1
                                elif p1.y() > p2.y():  # Right to left
                                    index = abs(edge.localIndex)
                                    polygon_item.insert(index + added - added_this, p)
                                    added_this += 1
                                    added += 1
                            else:
                                if p1.x() < p2.x():  # Left to right
                                    index = abs(edge.localIndex)
                                    polygon_item.insert(index + added, p)
                                    added += 1
                                elif p1.x() > p2.x():  # Right to left
                                    index = abs(edge.localIndex)
                                    polygon_item.insert(index + added - added_this, p)
                                    added_this += 1
                                    added += 1

        return self.poly_to_list(polygon_item, "Global")

    def polygon_contains(self, poly_outer, poly_inner):
        """
        Checks if a inner polygon is fully contained by a outer polygon, returns a list with boolean values of all
        points in the inner triangle which holds value true if contained and false else. Not that values on the border
        do not count as contained
        """
        inner_list = self.poly_to_list(poly_inner, "Global")
        contain_list = []

        # Loop over all points in the inner polygon to see if they are contained by the outer polygon
        for point in inner_list:
            # Points are defined in local coordinates, move them so they are both in the local coordinates
            # of the outer polygon
            p_x = point.x() - poly_outer.x()
            p_y = point.y() - poly_outer.y()
            point.setX(p_x)
            point.setY(p_y)

            # Check if the outer polygon contains none, some, or all of the points
            if poly_outer.contains(point):
                true_contain = []
                # check a square area around the point to see if the whole square is contained, else the point
                # is on a edge and should not be included
                for i, j in itertools.product(range(-1, 2), range(-1, 2)):
                    point.setX(p_x + i)
                    point.setY(p_y + j)
                    if poly_outer.contains(point):
                        true_contain.append(True)
                    else:
                        true_contain.append(False)
                # Add to contain_list if the whole square area is inside the outer polygon
                if all(true_contain):
                    contain_list.append(True)
                else:
                    contain_list.append(False)
            else:
                contain_list.append(False)
        return contain_list

    def on_action_merge(self):
        """
        Merge together overlapping polygons in the graphics view
        """
        # Loop over all polygons and compare to all other, if two polygons are merged they are removed from the list
        for poly_outer in self.poly_list:
            for poly_inner in self.poly_list:
                if poly_outer == poly_inner:
                    continue  # Ignore comparison to self

                contain_list = self.polygon_contains(poly_outer, poly_inner)

                if all(contain_list):
                    # If all points are inside the outer polygon do not merge (this would remove the inner one)
                    pass
                elif any(contain_list):
                    # If some but not all points are inside the outer polygon the two polygons overlap and should be
                    # merged

                    # Ignore holes
                    if poly_inner in self.hole_list or poly_outer in self.hole_list:
                        pass
                    else:
                        # Move the QPolygonF items to the global coordinates and unite them (merge)
                        p1 = poly_outer.polygon().translated(poly_outer.x(), poly_outer.y())
                        p2 = poly_inner.polygon().translated(poly_inner.x(), poly_inner.y())
                        uni = p1.united(p2)

                        # Unite adds the starting point again as endpoint so we have to remove this duplicate point to avoid
                        # future problems
                        uni = self.poly_to_list(uni, "Global")
                        uni = uni[:-1]

                        # Add the new merged polygon, remove the old polygons from the view and lists
                        self.add_poly_to_scene(QPolygonF(uni))
                        self.delete_polygon(poly_inner, True)
                        self.delete_polygon(poly_outer, True)
                        break
                else:
                    print("No overlap") # for debugging, remove later

    def toggle_gridsnap(self):
        self.grid_snap = not self.grid_snap

    def create_grid(self):

        # If called when a grid already exists create a new grid
        if self.grid:
            self.grid = []

        grid_pen = QPen(QColor(215, 215, 215), 1)
        w = 10000
        h = 10000
        self.addLine(-10000, 0, 10000, 0, QPen(QColor(0, 0, 0), 2))
        self.addLine(0, -10000, 0, 10000, QPen(QColor(0, 0, 0), 2))
        w = int(w / self.grid_spacing) * self.grid_spacing
        h = int(h / self.grid_spacing) * self.grid_spacing
        for i in range(-w, w, self.grid_spacing):
            if i == 0:
                pass
            else:
                self.grid.append(self.addLine(-w, i, w, i, grid_pen))
        for i in range(-h, h, self.grid_spacing):
            if i == 0:
                pass
            else:
                self.grid.append(self.addLine(i, -h, i, h, grid_pen))

        self.grid_built = True

    def toggle_grid(self, state=None):
        """
            Toggles the view of the grid
        """
        if not self.grid_built:
            self.create_grid()

        # If no state is provided toggle depending on current state
        if not state:
            if self.grid_on:
                state = "off"
            else:
                state = "on"

        if state == "on":
            for item in self.grid:
                item.setVisible(True)
                item.setZValue(-1)  # Make sure grid is always furthest back (only items with negative z)
            self.grid_on = True
        elif state == "off":
            for item in self.grid:
                item.setVisible(False)
            self.grid_on = False

    def build_gmsh(self):
        g = cfg.Geometry()
        point_index = 0
        line_index = 0
        added_points = []
        added_lines = []
        point_marker_list = []
        edge_marker_list = []
        surface_index = 0
        ignore_warning = False

        # Create a list of all coordinates where a point has a marker
        for point in self.point_marker_list:
            point_marker_list.append(point.scenePos())

        # Create a list with coordinates of the points of a edge which has a marker
        for edge in self.line_marker_list:
            l = edge.line()
            x = edge.scenePos().x()
            y = edge.scenePos().y()
            edge_marker_list.append([[l.p1().x() + x, -l.p1().y() - y], [l.p2().x() + x, - l.p2().y() - y]])

        # Loop first one time to sort the list in order of polygons containing each other to ensure that polygons which
        # contain other polygons are added first
        for poly in self.poly_list:
            if self.polygon_contains_other_polygon(poly):
                pass
            else:
                self.poly_list.append(self.poly_list.pop(self.poly_list.index(poly)))

        # Loop through all polygons again to add them to the gmsh object
        for poly in self.poly_list:
            if poly in self.hole_list:
                continue

            if self.polygon_overlaps_other_polygon(poly) and not ignore_warning:
                if self.overlap_warning is not None:
                    user_choice = self.overlap_warning()
                    if user_choice == "Cancel":
                        return None
                    elif user_choice == "Ignore":
                        ignore_warning = True

            surface_tot = []
            point_index, line_index, added_points, added_lines, surface = self.add_polygon_to_gmsh(poly, g, point_index, point_marker_list, line_index, added_points, added_lines, edge_marker_list)
            surface_tot.append(surface)
            surface_holes = []

            # Check if there are any holes inside the polygon
            for hole in self.polygon_contains_holes(poly):
                point_index, line_index, added_points, added_lines, surface_hole = self.add_polygon_to_gmsh(hole, g, point_index, point_marker_list, line_index, added_points, added_lines, edge_marker_list)
                surface_holes.append(surface_hole)

            for inner_poly in self.polygon_contains_other_polygon(poly):
                point_index, line_index, added_points, added_lines, surface_inner = self.add_polygon_to_gmsh(inner_poly, g,
                                                                                                            point_index,
                                                                                                            point_marker_list,
                                                                                                            line_index,
                                                                                                            added_points,
                                                                                                            added_lines,
                                                                                                            edge_marker_list)
                surface_holes.append(surface_inner)
                surface_tot.append(surface_inner)
            # Create a surface from the added lines, give the surface a number index
            g.surface(surface, holes=surface_holes, marker=surface_index)
            surface_index += 1
        self.g = g
        return g

    def add_polygon_to_gmsh(self, poly, g, point_index, point_marker_list, line_index, added_points, added_lines, edge_marker_list):
        first = True
        first_index = point_index
        prev_index = point_index
        poly_list = self.on_action_split(poly)
        surface = []

        for point in poly_list:  # loop over all points in the polygon
            point_marker = None

            # If the point is in the marker list store the marker text in point_marker
            if point in point_marker_list:
                i = point_marker_list.index(point)
                point_marker = self.point_marker_list[i].marker

            if [point.x(), -point.y()] in added_points:  # if the point already exists don't add again
                if first:  # If it is the first point of the poly set existing point as first_index
                    first = False
                    current_index = added_points.index([point.x(), -point.y()])
                    first_index = current_index
                elif prev_index == added_points.index([point.x(), -point.y()]):
                    # To avoid occurence of adding zero length line from point to self
                    pass
                else:  # Check if the line from previous point already exists, in that case don't att a new one
                    current_index = added_points.index([point.x(), -point.y()])
                    if [prev_index, current_index] in added_lines:
                        surface.append(added_lines.index([prev_index, current_index]))
                    elif [current_index, prev_index] in added_lines:
                        surface.append(added_lines.index([current_index, prev_index]))
                    else:  # else add a line from the previous point to the current pre-existing one
                        # If the coordinates is in the edge_marker_list also add the marker for that edge
                        if [added_points[prev_index], added_points[current_index]] in edge_marker_list:
                            i = edge_marker_list.index([added_points[prev_index], added_points[current_index]])
                            line_marker = self.line_marker_list[i].marker
                            g.addSpline([prev_index, current_index], marker=line_marker)
                        elif [added_points[current_index], added_points[prev_index]] in edge_marker_list:
                            i = edge_marker_list.index([added_points[current_index], added_points[prev_index]])
                            line_marker = self.line_marker_list[i].marker
                            g.addSpline([prev_index, current_index], marker=line_marker)
                        else:
                            g.addSpline([prev_index, current_index])
                        surface.append(line_index)
                        line_index += 1
                        added_lines.append([prev_index, current_index])
            else:
                # Else add a new point, if the point has a marker also add the marker
                if point_marker:
                    g.addPoint([point.x(), -point.y()], point_index, marker=point_marker)
                else:
                    g.addPoint([point.x(), -point.y()],
                               point_index)  # Negative y since the Graphicsscene has inverted y-axis
                current_index = point_index
                point_index += 1
                added_points.append([point.x(), -point.y()])

                if first:
                    # If it is the first added point do nothing as there is no other point to connect a line to
                    first = False
                else:
                    # Else add a new line between the previous point and the current, if those coordinates
                    # correspond with a marker also add the marker to the line
                    if [added_points[prev_index], added_points[current_index]] in edge_marker_list:
                        i = edge_marker_list.index([added_points[prev_index], added_points[current_index]])
                        line_marker = self.line_marker_list[i].marker
                        g.addSpline([prev_index, current_index], marker=line_marker)
                    elif [added_points[current_index], added_points[prev_index]] in edge_marker_list:
                        i = edge_marker_list.index([added_points[current_index], added_points[prev_index]])
                        line_marker = self.line_marker_list[i].marker
                        g.addSpline([prev_index, current_index], marker=line_marker)
                    else:
                        g.addSpline([prev_index, current_index])
                    surface.append(line_index)
                    line_index += 1
                    added_lines.append([prev_index, current_index])
            prev_index = current_index

        # Finally do the procedure again but add the line from the last index to the first index to close the
        # surface, also here add a line marker if there is one
        if [current_index, first_index] in added_lines:
            surface.append(added_lines.index([current_index, first_index]))
        elif [first_index, current_index] in added_lines:
            surface.append(added_lines.index([first_index, current_index]))
        else:
            if [added_points[first_index], added_points[current_index]] in edge_marker_list:
                i = edge_marker_list.index([added_points[first_index], added_points[current_index]])
                line_marker = self.line_marker_list[i].marker
                g.addSpline([current_index, first_index], marker=line_marker)
            elif [added_points[current_index], added_points[first_index]] in edge_marker_list:
                i = edge_marker_list.index([added_points[current_index], added_points[first_index]])
                line_marker = self.line_marker_list[i].marker
                g.addSpline([current_index, first_index], marker=line_marker)
            else:
                g.addSpline([current_index, first_index])
            surface.append(line_index)
            line_index += 1
            added_lines.append([current_index, first_index])
        return point_index, line_index, added_points, added_lines, surface

    def polygon_contains_holes(self, outer_poly):
        contain_list = []
        for hole_polygon in self.hole_list:
            if all(self.polygon_contains(outer_poly, hole_polygon)):
                contain_list.append(hole_polygon)
        return contain_list

    def polygon_contains_other_polygon(self, outer_poly):
        contain_list = []
        for inner_poly in self.poly_list:
            if outer_poly == inner_poly:
                pass
            elif all(self.polygon_contains(outer_poly, inner_poly)):
                contain_list.append(inner_poly)
        return contain_list

    def polygon_overlaps_other_polygon(self, outer_poly):
        contain_list = []
        for inner_poly in self.poly_list:
            if outer_poly == inner_poly:
                pass
            elif all(self.polygon_contains(outer_poly, inner_poly)):
                pass
            elif any(self.polygon_contains(outer_poly, inner_poly)):
                contain_list.append(inner_poly)
        return contain_list

    def on_save_geometry_action(self):
        if self.g:
            name = QFileDialog.getSaveFileName(self.parent(), filter="*.cfg")
            print(name)
            if name != ('', ''):
                cfu.save_geometry(self.g, name[0])

    def on_save_mesh_action(self):
        if self.g:
            name = QFileDialog.getSaveFileName(self.parent(), filter="*.cfm")
            if name != ('', ''):
                cfu.save_mesh(self.g, name[0])

    def on_save_arrays_action(self):
        if self.mesh:
            name = QFileDialog.getSaveFileName(self.parent(), filter="*.cfma")
            if name != ('', ''):
                mesh = self.mesh
                mesh.return_boundary_elements = True
                coords, edof, dofs, bdofs, elementmarkers, boundary_elements = mesh.create()
                cfu.save_arrays(coords, edof, dofs, bdofs, elementmarkers, boundary_elements, self.line_marker_dict, name[0])

    def on_save_arrays_matlab_action(self):
        if self.mesh:
            name = QFileDialog.getSaveFileName(self.parent(), filter="*.mat")
            if name != ('', ''):
                mesh = self.mesh
                mesh.return_boundary_elements = True
                coords, edof, dofs, bdofs, elementmarkers, boundary_elements = mesh.create()
                cfu.save_matlab_arrays(coords, edof, dofs, bdofs, elementmarkers, boundary_elements, self.line_marker_dict, name[0])

    def show_geom(self):
        g = self.build_gmsh()
        cfv.clf()
        if g:
            cfv.draw_geometry(g, draw_points=self.draw_points, label_points=self.label_points, label_curves=self.label_curves)
            if self.figure_canvas is not None:
                self.figure_canvas.draw()
            else:
                cfv.show_and_wait()
            return None
        else:
            return "Canceled"


    def show_mesh(self):
        g = self.build_gmsh()
        mesh = cfm.GmshMesh(g)
        mesh.el_type = self.el_type

        mesh.dofs_per_node = self.dofs_per_node
        mesh.el_size_factor = self.el_size_factor
        self.mesh = mesh

        coords, edof, dofs, bdofs, elementmarkers = mesh.create()
        print(edof)
        cfv.clf()

        cfv.draw_mesh(
            coords=coords,
            edof=edof,
            dofs_per_node=mesh.dofs_per_node,
            el_type=mesh.el_type,
            filled=True
        )
        if self.figure_canvas is not None:
            self.figure_canvas.draw()
        else:
            cfv.show_and_wait()

    def mesh_scroll_area(self):
        """
        Update the ScrollArea to hold three input lines to allow the user to change element size, element type and
        dofs per node
        """
        scroll_area_widget_contents = QWidget()
        grid = QGridLayout(scroll_area_widget_contents)
        if self.scroll_area:
            self.scroll_area.setWidget(scroll_area_widget_contents)
        label_el_size = QLabel("Element Size")
        edit_el_size = QLineEdit(str(self.el_size_factor))
        #edit_el_size.setValidator(QDoubleValidator())
        locale = QLocale(QLocale.English, QLocale.UnitedStates)

        validator = QDoubleValidator()
        validator.setLocale(locale)
        validator.setNotation(QDoubleValidator.StandardNotation)
        edit_el_size.setValidator(validator)

        grid.addWidget(label_el_size, 0, 0)
        grid.addWidget(edit_el_size, 0, 1)

        label_el_type = QLabel("Element Type")
        edit_el_type = QComboBox()
        edit_el_type.addItem("Triangle")
        edit_el_type.addItem("Quadrangle")
        if self.el_type == 3:
            edit_el_type.setCurrentIndex(1)
        grid.addWidget(label_el_type, 1, 0)
        grid.addWidget(edit_el_type, 1, 1)

        label_dofs_per_node = QLabel("DOFs per node")
        edit_dofs_per_node = QLineEdit(str(self.dofs_per_node))
        edit_dofs_per_node.setValidator(QIntValidator())
        grid.addWidget(label_dofs_per_node, 2, 0)
        grid.addWidget(edit_dofs_per_node, 2, 1)

        def update():
            if edit_dofs_per_node.text() != "":
                self.edit_dofs_per_node = int(edit_dofs_per_node.text())
            if edit_el_size.text() != "":
                string = edit_el_size.text()
                string = string.replace(",", ".")
                self.el_size_factor = float(string)
            if edit_el_type.currentText() == "Triangle":
                self.el_type = 2
            elif edit_el_type.currentText() == "Quadrangle":
                self.el_type = 3
            self.show_mesh()
        update_button = (QPushButton("Update"))
        grid.addWidget(update_button, 3, 1)
        update_button.clicked.connect(update)


def init_app():
    app = QApplication.instance()

    if app is None:
        print("No QApplication instance found. Creating one.")
        # if it does not exist then a QApplication is created
        app = QApplication(sys.argv)
    else:
        print("QApplication instance found.")

    return app


def main_loop():
    app = QApplication.instance()

    if app is None:
        print("No QApplication instance found. Creating one.")
        # if it does not exist then a QApplication is created
        app = QApplication(sys.argv)
    else:
        print("QApplication instance found.")

    # För matplotlib kompatibilitet

    # plt.show(block=False)

    app.exec_()

def run_editor():

    app = init_app()

    widget = EditorWindow()
    widget.show()

    app.exec_()


def run_and_load(g : cfg.Geometry):
    # --- Skapa applikationsinstans

    app = init_app()

    # --- Skapa och visa huvudfönster

    widget = EditorWindow()
    widget.load_scene_from_g(g)
    widget.show()

    app.exec_()


def edit_geometry(g : cfg.Geometry = None):
    app = init_app()

    widget = EditorWindow()
    if g is not None:
        widget.load_scene_from_g(g)
    widget.return_g = True
    widget.show()
    app.exec_()
    new_g = widget.scene.g
    line_marker_dict = widget.scene.line_marker_dict
    if new_g:
        return new_g, line_marker_dict
    else:
        return g


if __name__ == '__main__':
    app = init_app()

    widget = EditorWindow()
    widget.show()

    sys.exit(app.exec_())


