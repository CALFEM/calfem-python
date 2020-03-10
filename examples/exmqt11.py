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
        
    def updateGeometry(self):
        """Update geometry from set parameters"""        
        rect = cfs.Rectangle(self.w, self.h, elementType=3, dofsPerNode=2, maxArea=self.maxArea)
        rect.t = self.t
        rect.v = self.v
        rect.E = self.E

        rect.ptype = 1
        rect.ep = [rect.ptype, rect.t]
        rect.D = cfc.hooke(rect.ptype, rect.E, rect.v)
        
        self.rect = rect
        
    def updateMesh(self):
        """Generate mesh"""
        self.updateGeometry()
        self.mesh = cfs.ShapeMesh(self.rect)    
        
    def solve(self):
        """Solve problem"""
        
        cfu.enableLogging()
        
        solver = cfslv.Plan2DSolver(self.mesh)
        
        solver.addBC(self.rect.leftId, 0.0)
        solver.addForceTotal(self.rect.topId, -10e5, dimension=2)
        
        self.results = solver.execute()       
        
        cfu.disableLogging()

        
    def drawGeometry(self, figGeometry):
        """Draw geometry in provided figure"""
        cfv.figure(figGeometry.nr) 
        cfv.clf()
        cfv.draw_geometry(self.rect.geometry(), title="Geometry")
        
    def drawMesh(self, figMesh):
        """Draw mesh in provided figure"""
        cfv.figure(figMesh.nr)
        cfv.clf()
        cfv.draw_mesh(self.mesh.coords, self.mesh.edof, self.rect.dofsPerNode, self.rect.elementType, 
                     filled=True, title="Mesh") #Draws the mesh.
        
    def drawDisplacements(self, figDisplacements):
        """Draw displacements in provided figure"""
        cfv.figure(figDisplacements.nr)
        cfv.clf()
        cfv.draw_displacements(self.results.a, self.mesh.coords, self.mesh.edof, self.rect.dofsPerNode, self.rect.elementType, draw_undisplaced_mesh=False, title="Displacements", magnfac=1)
        
    def drawElementValues(self, figElementValues):
        """Draw element values in provided figure"""
        cfv.figure(figElementValues.nr)
        cfv.clf()
        cfv.draw_element_values(self.results.elForces, self.mesh.coords, self.mesh.edof, self.rect.dofsPerNode, self.rect.elementType, self.results.a, draw_elements=True, draw_undisplaced_mesh=False, title="Effective Stress", magnfac=1)
                      

class MainWindow(QMainWindow):
    """Main window class of our UI"""
    def __init__(self):
        """Constructor"""
        super(MainWindow, self).__init__()

        # Load user interface from UI-file

        loadUi('exmqt11.ui', self)

        # Query for figure class name

        Figure = cfv.figureClass()

        # Create figure widgets to insert in UI

        self.figGeometry = Figure(self)
        self.figMesh = Figure(self)
        self.figElementValues = Figure(self)
        self.figDisplacements = Figure(self)

        # Insert widgets in gridLayout

        self.middleLayout.addWidget(self.figGeometry._widget, 20)
        self.middleLayout.addWidget(self.figMesh._widget, 20)
        self.resultLayout.addWidget(self.figElementValues._widget, 20)
        self.resultLayout.addWidget(self.figDisplacements._widget, 20)
        #self.gridLayout.addWidget(self.figElementValues._widget, 0, 1)
        #self.gridLayout.addWidget(self.figDisplacements._widget, 1, 0)

        # Create our problem instance

        self.problem = PlaneStress2DProblem()
        self.updateView()
        
    def updateView(self):
        """Update controls with values from problem"""
        self.lengthEdit.setText(str(self.problem.w))
        self.heightEdit.setText(str(self.problem.h))
        self.thicknessEdit.setText(str(self.problem.t))
        self.elasticModulusEdit.setText(str(self.problem.E))
        self.youngEdit.setText(str(self.problem.v))
        self.maxAreaEdit.setText(str(self.problem.maxArea))
        
    def updateProblem(self):
        """Update problem with values from controls"""
        self.problem.w = float(self.lengthEdit.text())
        self.problem.h = float(self.heightEdit.text())
        self.problem.t = float(self.thicknessEdit.text())
        self.problem.E = float(self.elasticModulusEdit.text())
        self.problem.v = float(self.youngEdit.text())
        self.problem.maxArea = float(self.maxAreaEdit.text())
        
    @pyqtSlot()
    def on_updateButton_clicked(self):
        """Update geometry"""

        self.updateProblem()        
        self.problem.updateMesh()

        # Draw geometry

        self.problem.drawGeometry(self.figGeometry)
        
        # Draw mesh
        
        self.problem.drawMesh(self.figMesh)
        
    @pyqtSlot(int)
    def on_tabWidget_currentChanged(self, tabIndex):
        """Handle tab change"""

        if tabIndex == 1:        
            
            # Only calculate when on tab 1 (results)
       
            self.updateProblem()        
            self.problem.updateMesh()
            self.problem.solve()
            
            # Draw geometry
    
            self.problem.drawGeometry(self.figGeometry)
            
            # Draw mesh
            
            self.problem.drawMesh(self.figMesh)
            
            # Draw results
            
            self.problem.drawElementValues(self.figElementValues)
            self.problem.drawDisplacements(self.figDisplacements)

if __name__ == "__main__":

    app = QApplication(sys.argv)
    widget = MainWindow()
    widget.show()
    sys.exit(app.exec_())
