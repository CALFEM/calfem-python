import numpy as np
import calfem.core as cfc
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis as cfv

import logging as cflog

def error(msg):
    cflog.error("calfem.shapes: "+msg)

def info(msg):
    cflog.info("calfem.shapes: "+msg)


class Shape:
    """Base class for shapes"""

    def __init__(self, elementType = 3, dofsPerNode = 1, maxArea = 1.0):
        self.maxArea = maxArea
        self.elementType = elementType
        self.dofsPerNode = dofsPerNode

    def geometry(self):
        """Return geometry of shape"""
        return None



class Rectangle(Shape):
    """Rectangle geometry"""
    
    def __init__(self, width=1.0, height=1.0, elementType = 3, dofsPerNode = 1, maxArea = -1):
        Shape.__init__(self, elementType, dofsPerNode, maxArea)
        self.width = width
        self.height = height
    
        if (maxArea<0):
            self.maxArea = self.width * self.height * 0.05
        else:
            self.maxArea = maxArea
        
        self.leftId = 101
        self.rightId = 102
        self.topId = 103
        self.bottomId = 104

    def geometry(self):
        """Return geometry of shape"""
        self.g = cfg.Geometry()

        w = self.width
        h = self.height

        self.g.point([0, 0])
        self.g.point([w, 0])
        self.g.point([w, h])
        self.g.point([0, h])

        self.g.spline([0, 1], marker=self.bottomId)
        self.g.spline([1, 2], marker=self.rightId)
        self.g.spline([2, 3], marker=self.topId)
        self.g.spline([3, 0], marker=self.leftId)

        self.g.surface([0,1,2,3])

        return self.g

class ShapeMesh:
    """Mesh generator for shapes"""
    def __init__(self, shape):
        """Initialise mesh generator"""
        self.shape = shape

        self.create()

    def create(self):
        meshGen = cfm.GmshMeshGenerator(self.shape.geometry())
        meshGen.elType = self.shape.elementType
        meshGen.elSizeFactor = self.shape.maxArea
        meshGen.dofsPerNode = self.shape.dofsPerNode

        self.coords, self.edof, self.dofs, self.bdofs, self.markers = meshGen.create()

        # --- Ber├ñkna element koordinater

        self.ex, self.ey = cfc.coordxtr(self.edof, self.coords, self.dofs)
        self.pointDofs = self.dofs[list(self.shape.g.points.keys()),:]

if __name__ == "__main__":


    print("Creating rectangle")
    rect = Rectangle(5.0, 1.0, elementType=3, dofsPerNode=2, maxArea=0.05)

    print("Creating mesh...")

    mesh = ShapeMesh(rect)

    print(mesh.edof)
    print(mesh.dofs)
    print(mesh.bdofs[rect.leftId])
    print(mesh.bdofs[rect.rightId])
    print(mesh.bdofs[rect.topId])
    print(mesh.bdofs[rect.bottomId])
    print(mesh.pointDofs)
