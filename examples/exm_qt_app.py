# -*- coding: utf-8 -*-

'''Example 11

Simple example illustrating how to integrate calfem for Python in 
a complex Qt user interface.
'''

import sys

from calfem.qt5 import *

import calfem.core as cfc
import calfem.vis as cfv
import calfem.utils as cfu
import calfem.shapes as cfs
import calfem.solver as cfslv

class PlaneStress2DProblem(object):
    def __init__(self, width=5.0, height=1.0, t=1.0, v=0.2, E=2e9, maxArea=0.08):
        """Init model parameters"""
        self.w = width
        self.h = height
        self.t = t
        self.v = v
        self.E = E
        self.maxArea = 0.08
        
    def update_geometry(self):
        """Update geometry from set parameters"""        
        rect = cfs.Rectangle(self.w, self.h, element_type=3, dofs_per_node=2, max_area=self.maxArea)
        rect.t = self.t
        rect.v = self.v
        rect.E = self.E

        rect.ptype = 1
        rect.ep = [rect.ptype, rect.t]
        rect.D = cfc.hooke(rect.ptype, rect.E, rect.v)
        
        self.rect = rect
        
    def update_mesh(self):
        """Generate mesh"""
        self.update_geometry()
        self.mesh = cfs.ShapeMesh(self.rect)    
        
    def solve(self):
        """Solve problem"""
        
        cfu.enableLogging()
        
        solver = cfslv.Plan2DSolver(self.mesh)
        
        solver.addBC(self.rect.left_id, 0.0)
        solver.addForceTotal(self.rect.top_id, -10e5, dimension=2)
        
        self.results = solver.execute()       
        
        cfu.disableLogging()

        
    def draw_geometry(self, figGeometry):
        """Draw geometry in provided figure"""
        cfv.figure(figGeometry.nr) 
        cfv.clf()
        cfv.draw_geometry(self.rect.geometry(), title="Geometry")
        
    def draw_mesh(self, figMesh):
        """Draw mesh in provided figure"""
        cfv.figure(figMesh.nr)
        cfv.clf()
        cfv.draw_mesh(self.mesh.coords, self.mesh.edof, self.rect.dofs_per_node, self.rect.element_type, 
                     filled=True, title="Mesh") #Draws the mesh.
        
    def draw_displacements(self, figDisplacements):
        """Draw displacements in provided figure"""
        cfv.figure(figDisplacements.nr)
        cfv.clf()
        cfv.draw_displacements(self.results.a, self.mesh.coords, self.mesh.edof, self.rect.dofs_per_node, self.rect.element_type, draw_undisplaced_mesh=False, title="Displacements", magnfac=1)
        
    def draw_elementValues(self, figElementValues):
        """Draw element values in provided figure"""
        cfv.figure(figElementValues.nr)
        cfv.clf()
        cfv.draw_element_values(self.results.el_forces, self.mesh.coords, self.mesh.edof, self.rect.dofs_per_node, self.rect.element_type, self.results.a, draw_elements=True, draw_undisplaced_mesh=False, title="Effective Stress", magnfac=1)
                      

class MainWindow(QMainWindow):
    """Main window class of our UI"""
    def __init__(self):
        """Constructor"""
        super(MainWindow, self).__init__()

        # Load user interface from UI-file

        loadUi('exm_qt_app.ui', self)

        # Query for figure class name

        figure = cfv.figureClass()

        # Create figure widgets to insert in UI

        self.fig_geometry = figure(self)
        self.fig_mesh = figure(self)
        self.fig_element_values = figure(self)
        self.fig_displacements = figure(self)

        # Insert widgets in gridLayout

        self.middleLayout.addWidget(self.fig_geometry._widget, 20)
        self.middleLayout.addWidget(self.fig_mesh._widget, 20)
        self.resultLayout.addWidget(self.fig_element_values._widget, 20)
        self.resultLayout.addWidget(self.fig_displacements._widget, 20)

        # Create our problem instance

        self.problem = PlaneStress2DProblem()
        self.update_view()
        
    def update_view(self):
        """Update controls with values from problem"""
        self.lengthEdit.setText(str(self.problem.w))
        self.heightEdit.setText(str(self.problem.h))
        self.thicknessEdit.setText(str(self.problem.t))
        self.elasticModulusEdit.setText(str(self.problem.E))
        self.youngEdit.setText(str(self.problem.v))
        self.maxAreaEdit.setText(str(self.problem.maxArea))
        
    def update_problem(self):
        """Update problem with values from controls"""
        self.problem.w = float(self.lengthEdit.text())
        self.problem.h = float(self.heightEdit.text())
        self.problem.t = float(self.thicknessEdit.text())
        self.problem.E = float(self.elasticModulusEdit.text())
        self.problem.v = float(self.youngEdit.text())
        self.problem.maxArea = float(self.maxAreaEdit.text())
        
    @Slot()
    def on_updateButton_clicked(self):
        """Update geometry"""

        self.update_problem()        
        self.problem.update_mesh()

        # Draw geometry

        self.problem.draw_geometry(self.fig_geometry)
        
        # Draw mesh
        
        self.problem.draw_mesh(self.fig_mesh)
        
    @Slot(int)
    def on_tabWidget_currentChanged(self, tabIndex):
        """Handle tab change"""

        if tabIndex == 1:        
            
            # Only calculate when on tab 1 (results)
       
            self.update_problem()        
            self.problem.update_mesh()
            self.problem.solve()
            
            # Draw geometry
    
            self.problem.draw_geometry(self.fig_geometry)
            
            # Draw mesh
            
            self.problem.draw_mesh(self.fig_mesh)
            
            # Draw results
            
            self.problem.draw_elementValues(self.fig_element_values)
            self.problem.draw_displacements(self.fig_displacements)

if __name__ == "__main__":

    app = QApplication(sys.argv)
    widget = MainWindow()
    widget.show()
    sys.exit(app.exec_())
