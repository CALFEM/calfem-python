# -*- coding: utf-8 -*-
"""
Created on Tue May 31 00:18:15 2016

@author: Jonas Lindemann
"""

import calfem.core as cfc
import calfem.geometry as cfg
import calfem.mesh as cfm

class Shape:
    """Base class for shapes"""
    
    def __init__(self):
        self.maxArea = 1.0
        self.elementType = 3
        self.dofsPerNode = 1
    
    def geometry(self):
        """Return geometry of shape"""
        return None
    
    
class Rectangle(Shape):
    """Rectangle geometry"""
    def __init__(self, width, height):
        Shape.__init__(self)
        self.width = width
        self.height = height
        self.maxArea = self.width * self.height * 0.05
        
    def geometry(self):
        """Return geometry of shape"""
        g = cfg.Geometry()
        
        w = self.width
        h = self.height

        g.point([0, 0])
        g.point([w, 0])
        g.point([w, h])
        g.point([0, h])
        
        g.spline([0, 1])
        g.spline([1, 2])
        g.spline([2, 3])
        g.spline([3, 0])
        
        g.surface([0,1,2,3])
        
        return g
        
class ShapeMesher:
    """Mesh generator for shapes"""
    def __init__(self, shape):
        """Initialise mesh generator"""
        self.shape = shape
        
    def create(self):
        meshGen = cfm.GmshMeshGenerator(self.shape.geometry())
        meshGen.elSizeFactor = self.shape.maxArea
        meshGen.elType = self.shape.elementType
        meshGen.dofsPerNode = self.shape.dofsPerNode
        
        self.coords, self.edof, self.dofs, self.bdofs, self.markers = meshGen.create()
        #self.topo = meshGen.topo
    
        # --- Beräkna element koordinater
        
        self.ex, self.ey = cfc.coordxtr(self.edof, self.coords, self.dofs)        

if __name__ == "__main__":
    
    print("Creating rectangle")
    rect = Rectangle(2.0, 1.0)
    
    print("Creating mesh...")

    mesher = ShapeMesher(rect)    
    mesher.create()
    
    print(mesher.edof)

    
    
