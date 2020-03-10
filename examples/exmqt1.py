# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 16:23:37 2017

@author: Jonas Lindemann
"""

import sys
import time

from calfem.qt5 import *

import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis as cfv

class MainWindow(QMainWindow):
    def __init__(self):
        """Constructor"""
        super(MainWindow, self).__init__()

        # Load user interface from UI-file

        loadUi('exmqt1.ui', self)
        
        Figure = cfv.figureClass()
        
        self.fig1 = Figure(self)
        self.fig2 = Figure(self)
        
        self.gridLayout.addWidget(self.fig1._widget, 0, 0)
        self.gridLayout.addWidget(self.fig2._widget, 0, 1)
        
        print(self.fig1._widget)
        
      
    @pyqtSlot()
    def on_executeButton_clicked(self):
        self.solveProblem()        
        self.drawGeometry()
        self.drawMesh()
        
        
    def solveProblem(self):
                
        g = cfg.Geometry() #Create a GeoData object that holds the geometry.

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

        g.ellipse([7,8,9,10], marker=50)      # 0 - An ellipse arc. Read the function 
        g.spline([0, 1], marker=80)           # 1 - A spline. Splines pass through the 
        g.spline([2, 1])                      # 2
        g.spline([3, 2])                      # 3
        g.spline([0, 3])                      # 4
        g.spline([7, 9], marker=50)           # 5
        g.spline([10, 9])                     # 6
        g.spline([4, 5, 6, 4])                # 7 - This is a closed spline. 

        g.surface([4,3,2,1], [[7], [5,6,0]])
        
        meshGen = cfm.GmshMeshGenerator(g)
              
        meshGen.elType = 3 # Degrees of freedom per node.
        meshGen.dofsPerNode = 1 # Factor that changes element sizes.
        meshGen.elSizeFactor = 0.05     
                
        self.coords, self.edof, self.dofs, self.bdofs, self.elementmarkers = meshGen.create()
        self.meshGen = meshGen       
        self.g = g
               
    def drawGeometry(self):
        cfv.figure(self.fig1.nr) 
        cfv.clf()
        cfv.drawGeometry(self.g)
    
    def drawMesh(self):
        cfv.figure(self.fig2.nr) 
        cfv.clf()
        cfv.draw_mesh(
            coords=self.coords, 
            edof=self.edof, 
            dofs_per_node=self.meshGen.dofsPerNode, 
            el_type=self.meshGen.elType, 
            filled=True, 
            title="Example 01"
            ) 
        
        cfv.addText("This is a Text", pos=(1, -0.3), angle=45)  #Adds a text in world space
        
        ourLabel = cfv.label("This is a Label", pos=(100,200), angle=-45) #Adds a label in the screen space
        ourLabel.text = "Label, changed." #We can change the attributes of labels and texts, such as color, text, and position.
        ourLabel.textColor = 'r'  #Make it red. (1,0,0) would also have worked.
        ourLabel.position = (20,30)
        
        
if __name__ == "__main__":
    
    app = QApplication(sys.argv)
    widget = MainWindow()
    widget.show()
    sys.exit(app.exec_())  
