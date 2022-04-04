# -*- coding: iso-8859-15 -*-
"""
CALFEM Editor Module

Contains functions for implementing a interactive geometry editor.

Written by Karl Eriksson
"""

import os, sys

from PyQt5.QtCore import Qt

from PyQt5.QtGui import QPainter, QCloseEvent

from PyQt5.QtWidgets import QApplication, QMainWindow, QGraphicsView, QButtonGroup, QMessageBox

from PyQt5.uic import loadUi

import calfem.geometry as cfg
import calfem.vis_mpl as cfv
import calfem.editor_scene as editor_scene

app = None


class EditorWindow(QMainWindow):
    """MainWindow-klass som hanterar vårt huvudfönster"""

    def __init__(self):
        super(QMainWindow, self).__init__()
        self.app = app

        root = os.path.dirname(os.path.realpath(__file__))
        loadUi(os.path.join(root, 'editor.ui'), self)

        # loadUi('editor.ui', self)   loadUI kan ladd ui-fil och lägga till objekt direkt i klassen.
        self.setWindowTitle("CALFEM Geometry Editor")

        scene = editor_scene.EditorScene(self)

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

        # Create buttongroup for buttons in surface view
        self.arrowButtonSurface.setDown(True)
        self.selectGroupSurface = QButtonGroup(self)
        self.selectGroupSurface.addButton(self.arrowButtonSurface)
        self.selectGroupSurface.addButton(self.rectButton)
        self.selectGroupSurface.addButton(self.polyButton)
        self.selectGroupSurface.addButton(self.addPolyHoleButton)
        self.selectGroupSurface.addButton(self.addRectHoleButton)
        self.selectGroupSurface.addButton(self.deleteButton)
        self.selectGroupSurface.addButton(self.panningButtonSurface)

        # Create buttongroup for buttons in border view
        self.selectGroupBorder = QButtonGroup(self)
        self.selectGroupBorder.addButton(self.arrowButtonBorder)
        self.selectGroupBorder.addButton(self.panningButtonBorder)
        self.selectGroupBorder.addButton(self.setMarkerButton)
        self.selectGroupBorder.addButton(self.splitEdgeButton)
        self.selectGroupBorder.buttonClicked.connect(self.clear_selection)

        # Connect all buttons to corresponding action to perform when pressed
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
        self.deleteButton.clicked.connect(self.on_action_delete)
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

        # Connect methods that affect the window based on actions in the graphics scene
        self.gridSpacingSpinBox.setValue(self.scene.grid_spacing)
        self.scene.update_labels = self.update_labels
        self.scene.toggle_tab_enabled = self.toggle_tab_enabled
        self.scene.get_graphics_view_size = self.get_graphics_view_size
        self.scene.overlap_warning = self.overlap_warning
        self.scene.marker_removal_warning = self.marker_removal_warning
        self.scene.intersection_error = self.intersection_error
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

        # Set the minimum and maximum allowed zooming distance, factors by 2
        self.max_view_scale = 8
        self.min_view_scale = 0.25

        self.graphicsView.wheelEvent = self.override_scroll

        self.text_browser_active = True
        self.showTextBrowserButton.setVisible(False)

        self.show()
        self.raise_()

    def toggle_text_browser(self):
        """ Toggle the text browser to be expanded or collapsed """
        if self.text_browser_active:
            self.clearTextBrowserButton.setVisible(False)
            self.hideTextBrowserButton.setVisible(False)
            self.showTextBrowserButton.setVisible(True)
            self.textBrowser.setMaximumHeight(4)

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
        # Control the amount of zoom from one scroll call
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

    def clear_selection(self):
        self.scene.clearSelection()
        self.scene.reset_scroll_area()

    def set_view_scale(self, factor):

        new_scale = self.graphicsView.transform().m11() * factor

        if self.min_view_scale <= new_scale <= self.max_view_scale:
            self.graphicsView.scale(factor, factor)

    def set_tooltip(self, string):
        self.labelTooltip.setText(string)

    def toggle_grid_snap(self):
        self.grid_snap_on = not self.grid_snap_on
        self.scene.toggle_gridsnap()

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Shift:
            if self.scene.mode == "Arrow" and not self.scene.selectedItems():
                self.scene.mode = "None"
                self.graphicsView.setDragMode(QGraphicsView.ScrollHandDrag)

    def keyReleaseEvent(self, event):
        if self.scene.mode == "None" and not self.scene.selectedItems():
            if event.key() == Qt.Key_Shift:
                self.scene.mode = "Arrow"
                self.graphicsView.setDragMode(QGraphicsView.NoDrag)

    def overlap_warning(self):
        msg = QMessageBox()
        msg.setWindowTitle("Warning!")
        msg.setText("There exists overlapping geometries")
        msg.setIcon(QMessageBox.Warning)
        msg.setStandardButtons(QMessageBox.Cancel | QMessageBox.Ignore)

        msg.buttonClicked.connect(self.popup_button)
        msg.exec_()
        return self.overlap_warning_choice

    def marker_removal_warning(self):
        msg = QMessageBox()
        msg.setWindowTitle("Warning!")
        msg.setText("Performing merging will remove markers")
        msg.setIcon(QMessageBox.Warning)
        msg.setStandardButtons(QMessageBox.Cancel | QMessageBox.Ignore)

        msg.buttonClicked.connect(self.popup_button)
        msg.exec_()
        return self.overlap_warning_choice

    def intersection_error(self):
        msg = QMessageBox()
        msg.setWindowTitle("Error!")
        msg.setText("There exists a polygon overlapping itself")
        msg.setIcon(QMessageBox.Critical)
        msg.setStandardButtons(QMessageBox.Cancel)

        msg.exec_()

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
        """Action performed when clicking one of the tabs in the tab bar, index of tab clicked determines action"""
        self.set_tooltip("")
        index = self.tabWidget.currentIndex()
        self.scene.set_view(index)
        self.on_action_arrow()  # Set arrow pointer as the standard tool

        if self.scene.canceled:  # Set to surface mode in case of cancelling due to errors
            index = 0
            self.scene.canceled = False
            self.tabWidget.setCurrentIndex(0)

        if index == 0:  # Set to surface mode
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

        elif index == 1:  # Set to border mode
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
        # If closed and no CALFEM Geometry has been generated make sure one is with the currently drawn geometry
        if self.return_g:
            self.scene.build_gmsh()

    def load_scene_from_g(self, g):
        self.scene.load_from_g(g)

    def on_action_panning(self):
        if self.scene.mode == "Draw Poly":
            self.scene.remove_drawing_poly()
        elif self.scene.mode == "Draw Rect":
            self.scene.remove_drawing_rect()
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
        self.graphicsView.viewport().setCursor(Qt.CrossCursor)
        self.graphicsView.setCursor(Qt.CrossCursor)

        self.set_tooltip("Press left mouse button to select starting point of polygon")
        for button in self.selectGroupSurface.buttons():
            button.setDown(False)
        self.polyButton.setDown(True)
        if self.scene.mode == "Draw Poly" and self.scene.hole_mode:
            self.scene.remove_drawing_poly()
        elif self.scene.mode == "Draw Rect":
            self.scene.remove_drawing_rect()
        self.scene.mode = "Draw Poly"
        self.scene.hole_mode = False

        if self.panning_on:
            self.toggle_panning()

    def on_action_poly_hole(self):
        """
        De-selects all other toolButtons and sets the current mode to drawing polygon hole
        """
        self.graphicsView.viewport().setCursor(Qt.CrossCursor)
        self.graphicsView.setCursor(Qt.CrossCursor)
        self.set_tooltip("Press left mouse button to select starting point of polygon hole")
        for button in self.selectGroupSurface.buttons():
            button.setDown(False)
        self.addPolyHoleButton.setDown(True)
        if self.scene.mode == "Draw Poly" and not self.scene.hole_mode:
            self.scene.remove_drawing_poly()
        elif self.scene.mode == "Draw Rect":
            self.scene.remove_drawing_rect()
        self.scene.mode = "Draw Poly"
        self.scene.hole_mode = True

        if self.panning_on:
            self.toggle_panning()

    def on_action_rect(self):
        """
        De-selects all other toolButtons and sets the current mode to drawing rectangle
        """
        self.graphicsView.viewport().setCursor(Qt.CrossCursor)
        self.graphicsView.setCursor(Qt.CrossCursor)
        self.set_tooltip("Press left mouse button to select starting point of rectangle")
        for button in self.selectGroupSurface.buttons():
            button.setDown(False)
        self.rectButton.setDown(True)
        if self.scene.mode == "Draw Poly":
            self.scene.remove_drawing_poly()
        elif self.scene.mode == "Draw Rect" and self.scene.hole_mode:
            self.scene.remove_drawing_rect()
        self.scene.mode = "Draw Rect"
        self.scene.hole_mode = False

        if self.panning_on:
            self.toggle_panning()

    def on_action_rect_hole(self):
        """
        De-selects all other toolButtons and sets the current mode to drawing rectangle hole
        """
        self.graphicsView.viewport().setCursor(Qt.CrossCursor)
        self.graphicsView.setCursor(Qt.CrossCursor)
        self.set_tooltip("Press left mouse button to select starting point of rectangular hole")
        for button in self.selectGroupSurface.buttons():
            button.setDown(False)
        self.addRectHoleButton.setDown(True)
        if self.scene.mode == "Draw Poly":
            self.scene.remove_drawing_poly()
        elif self.scene.mode == "Draw Rect" and not self.scene.hole_mode:
            self.scene.remove_drawing_rect()
        self.scene.mode = "Draw Rect"

        self.scene.hole_mode = True

        if self.panning_on:
            self.toggle_panning()

    def on_action_delete(self):
        """
        De-selects all other toolButtons and sets the current mode to drawing polygon
        """
        self.graphicsView.viewport().setCursor(Qt.PointingHandCursor)
        self.graphicsView.setCursor(Qt.PointingHandCursor)
        self.set_tooltip("Press item to delete")
        for button in self.selectGroupSurface.buttons():
            button.setDown(False)
        self.deleteButton.setDown(True)
        if self.scene.mode == "Draw Poly":
            self.scene.remove_drawing_poly()
        elif self.scene.mode == "Draw Rect":
            self.scene.remove_drawing_rect()
        self.scene.mode = "Delete"

        if self.panning_on:
            self.toggle_panning()

    def on_action_split(self):
        """
        De-selects all other toolButtons and sets the current mode to split line
        """
        self.set_tooltip("Click on existing edge to split by adding an extra node")
        self.scene.mode = "Split Line"
        self.graphicsView.viewport().setCursor(Qt.ArrowCursor)
        self.graphicsView.setCursor(Qt.ArrowCursor)

        for button in self.selectGroupBorder.buttons():
            button.setDown(False)
        self.splitEdgeButton.setDown(True)

        if self.panning_on:
            self.toggle_panning()

    def on_action_set_marker(self):
        """
        De-selects all other toolButtons and sets the current mode to set marker
        """
        self.set_tooltip("Select edge to add marker")
        self.scene.mode = "Set Marker"
        self.graphicsView.viewport().setCursor(Qt.PointingHandCursor)
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
        self.graphicsView.viewport().setCursor(Qt.ArrowCursor)
        self.graphicsView.setCursor(Qt.ArrowCursor)

        for button in self.selectGroupSurface.buttons():
            button.setDown(False)
        for button in self.selectGroupBorder.buttons():
            button.setDown(False)
        self.arrowButtonSurface.setDown(True)
        self.arrowButtonBorder.setDown(True)

        if self.scene.mode == "Draw Poly":
            self.scene.remove_drawing_poly()
        elif self.scene.mode == "Draw Rect":
            self.scene.remove_drawing_rect()

        self.scene.mode = "Arrow"

        if self.panning_on:
            self.toggle_panning()


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

    # --- Skapa och visa huvudf�nster

    widget = EditorWindow()
    widget.load_scene_from_g(g)
    widget.show()

    app.exec_()


def edit_geometry(g : cfg.Geometry = None):
    """
    Start interactive geometry editor
    
    :param g cfg.Geometry: Geometry to modify interactively ()
    :return cfg.Geometry new_g, line_marker_dict: Modfied geometry and marker dictionary.
    """
    app = init_app()

    widget = EditorWindow()
    if g is not None:
        widget.load_scene_from_g(g)
    widget.return_g = True
    widget.show()
    app.exec_()
    new_g = widget.scene.g
    line_marker_dict = widget.scene.marker_dict
    if new_g:
        return new_g, line_marker_dict
    else:
        return g


if __name__ == '__main__':
    app = init_app()

    widget = EditorWindow()
    widget.show()

    sys.exit(app.exec_())


