# -*- coding: utf-8 -*-

'''Example 06

Solves a plane stress 2D problem using a structured mesh.
Shows how to draw von Mises effective stress as an element value with
drawElementValues(). Shows use of GmshMesher attribute 'nodesOnCurve'
(dictionary that says which nodes are on a given geometry curve)
'''

import sys
from math import sqrt

from calfem.qt5 import *

import numpy as np

import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis as cfv
import calfem.utils as cfu
import calfem.core as cfc

#cfu.enableLogging()

class PlaneStress2DProblem(object):
    """Encapsulates our plane stress problem in a single class"""
    def __init__(self, t=0.2, v=0.35, E=2.1e9, ptype=1):
        """Setup intial problem values"""
        self.t = t
        self.v = v
        self.E = E
        self.ptype = ptype

    def updateGeometry(self):
        """Update geometry"""

        cfu.info("Creating geometry...")

        self.g = cfg.geometry()
        g = self.g

        # Just a shorthand. We use this to make the circle arcs.

        s2 = 1/sqrt(2)

        points = [[0, 3], [2.5, 3], [3, 3], [4-s2, 3-s2], [4, 2],     #0-4
                  [4+s2, 3-s2], [5, 3], [5.5, 3], [8, 3], [0, 1.5],    #5-9
                  [2.5, 1.5], [4, 1.5], [5.5, 1.5], [8, 1.5], [0, 0], #10-14
                  [2.5, 0], [3, 0], [4-s2, s2], [4, 1], [4+s2, s2],   #15-19
                  [5, 0], [5.5, 0], [8, 0], [4, 3], [4, 0]]              #20-24

        for xp, yp in points:
            g.point([xp*0.1, yp*0.1])

        splines = [[0, 1], [1, 2], [6, 7], [7, 8], [8, 13],          #0-4
                   [13, 22], [22, 21], [21, 20], [16, 15], [15, 14], #5-9
                   [14, 9], [9, 0], [9, 10], [10, 1], [10, 15],     #10-14
                   [10, 11], [11, 4], [11, 18], [11, 12], [12, 7],   #15-19
                   [12, 21], [12, 13], [3, 10], [5, 12], [10, 17],   #20-24
                   [12, 19]]                                     #25

        for s in splines:
            g.spline(s, el_on_curve=5)

        g.curveMarker(ID=4, marker=7) #Assign marker 7 to the splines on the right.
        g.curveMarker(ID=5, marker=7) # We will apply a force on nodes with marker 7.
        g.curveMarker(ID=10, marker=5) #Assign marker 5 to the splines on the left.
        g.curveMarker(ID=11, marker=5) # The nodes with marker 5 will be locked in place.

        # Points in circle arcs are [start, center, end]

        circlearcs = [[2, 23, 3], [3, 23, 4], [4, 23, 5], [5, 23, 6],           #26-29
                      [16, 24, 17], [17, 24, 18], [18, 24, 19], [19, 24, 20]]   #30-33

        for c in circlearcs:
            g.circle(c, el_on_curve=5)

        g.structuredSurface([11, 12, 13, 0]) #0
        g.structuredSurface([14, 12, 10, 9])
        g.structuredSurface([8, 30, 24, 14])
        g.structuredSurface([24, 31, 17, 15])
        g.structuredSurface([15, 16, 27, 22]) #4
        g.structuredSurface([22, 26, 1, 13])
        g.structuredSurface([16, 18, 23, 28])
        g.structuredSurface([19, 2, 29, 23])
        g.structuredSurface([19, 21, 4, 3]) #8
        g.structuredSurface([20, 6, 5, 21])
        g.structuredSurface([25, 20, 7, 33])
        g.structuredSurface([32, 17, 18, 25]) #11

    def updateMesh(self):
        """Update mesh"""

        cfu.info("Meshing geometry...")

        # 3 Quads

        self.el_type = 3
        self.dofs_per_node = 2

        # Create mesh

        self.mesh = cfm.GmshMeshGenerator(geometry=self.g)
        self.mesh.el_type = self.el_type
        self.mesh.dofs_per_node = self.dofs_per_node

        self.coords, self.edof, self.dofs, self.bdofs, self.elementmarkers = self.mesh.create()

    def solve(self):
        """Solve problem"""

        self.updateGeometry()
        self.updateMesh()

        self.ep = [self.ptype, self.t]
        self.D = cfc.hooke(self.ptype, self.E, self.v)

        cfu.info("Assembling system matrix...")

        nDofs = np.size(self.dofs)
        ex, ey = cfc.coordxtr(self.edof, self.coords, self.dofs)
        K = np.zeros([nDofs, nDofs])

        for eltopo, elx, ely in zip(self.edof, ex, ey):
            Ke = cfc.planqe(elx, ely, self.ep, self.D)
            cfc.assem(eltopo, K, Ke)

        cfu.info("Solving equation system...")

        f = np.zeros([nDofs, 1])

        bc = np.array([], 'i')
        bcVal = np.array([], 'i')

        bc, bcVal = cfu.applybc(self.bdofs, bc, bcVal, 5, 0.0, 0)

        cfu.applyforce(self.bdofs, f, 7, 10e5, 1)

        self.a, self.r = cfc.solveq(K, f, bc, bcVal)

        cfu.info("Computing element forces...")

        ed = cfc.extractEldisp(self.edof, self.a)
        self.vonMises = []

        # For each element:

        for i in range(self.edof.shape[0]):

            # Determine element stresses and strains in the element.

            es, et = cfc.planqs(ex[i, :], ey[i, :], self.ep, self.D, ed[i, :])

            # Calc and append effective stress to list.

            self.vonMises.append(sqrt(pow(es[0], 2) - es[0]*es[1] + pow(es[1], 2) + 3*es[2]))

            ## es: [sigx sigy tauxy]

    def drawGeometry(self, figGeometry):
        """Draws geometry in specified figure"""
        cfv.figure(figGeometry.nr)
        cfv.clf()
        cfv.draw_geometry(self.g, draw_points=False, label_curves=True)

    def drawElementValues(self, figElementValues):
        """Draws element values in specified figure"""
        cfv.figure(figElementValues.nr)
        cfv.clf()
        cfv.draw_element_values(self.vonMises, self.coords, self.edof,
                              self.mesh.dofs_per_node, self.mesh.el_type,
                              self.a, draw_elements=True,
                              draw_undisplaced_mesh=False,
                              title="Example 06 effective stress")

    def drawDisplacements(self, figDisplacements):
        """Draws displacements in specified figure"""
        cfv.figure(figDisplacements.nr)
        cfv.clf()
        cfv.draw_displacements(self.a, self.coords, self.edof,
                              self.mesh.dofs_per_node, self.mesh.el_type,
                              draw_undisplaced_mesh=True, title="Example 06")

        # Make use of attribute 'nodesOnCurve' in GmshMesher to draw some arrows on
        # the right hand side of the mesh:

        rightSideNodes = set()

        # 4 and 5 are the IDs of the curves where we applied the forces.

        for curveID in [4, 5]:

            # Get the nodes, without duplicates.

            rightSideNodes = rightSideNodes.union(set(self.mesh.nodesOnCurve[curveID]))

        for i in rightSideNodes:

            # Position of the node with displacements.

            x = self.coords[i, 0] + self.a[i*2, 0]
            y = self.coords[i, 1] + self.a[i*2+1, 0]

            # A poor man's force indicator. Could also use vv.plot()

            cfv.addText("\rightarrow", (x, y), font_size=20, color='g')


class MainWindow(QMainWindow):
    """Main window class of our UI"""
    def __init__(self):
        """Constructor"""
        super(MainWindow, self).__init__()

        # Load user interface from UI-file

        loadUi('exmqt6.ui', self)

        # Query for figure class name

        Figure = cfv.figureClass()

        # Create figure widgets to insert in UI

        self.figGeometry = Figure(self)
        self.figElementValues = Figure(self)
        self.figDisplacements = Figure(self)

        # Insert widgets in gridLayout

        self.gridLayout.addWidget(self.figGeometry._widget, 0, 0)
        self.gridLayout.addWidget(self.figElementValues._widget, 0, 1)
        self.gridLayout.addWidget(self.figDisplacements._widget, 1, 0)

        # Create our problem instance

        self.problem = PlaneStress2DProblem()

    @pyqtSlot()
    def on_executeButton_clicked(self):
        """Execute calculation"""

        # Solve problem

        self.problem.solve()

        # Draw results

        self.problem.drawGeometry(self.figGeometry)
        self.problem.drawElementValues(self.figElementValues)
        self.problem.drawDisplacements(self.figDisplacements)


if __name__ == "__main__":

    app = QApplication(sys.argv)
    widget = MainWindow()
    widget.show()
    sys.exit(app.exec_())
