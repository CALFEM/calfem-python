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

    def __init__(self, element_type = 3, dofs_per_node = 1, max_area = 1.0):
        self.max_area = max_area
        self.element_type = element_type
        self.dofs_per_node = dofs_per_node

    def geometry(self):
        """Return geometry of shape"""
        return None



class Rectangle(Shape):
    """Rectangle geometry"""
    
    def __init__(self, width=1.0, height=1.0, element_type = 3, dofs_per_node = 1, max_area = -1):
        Shape.__init__(self, element_type, dofs_per_node, max_area)
        self.width = width
        self.height = height
    
        if (max_area<0):
            self.max_area = self.width * self.height * 0.05
        else:
            self.max_area = max_area
        
        self.left_id = 101
        self.right_id = 102
        self.top_id = 103
        self.bottom_id = 104

    def geometry(self):
        """Return geometry of shape"""
        self.g = cfg.Geometry()

        w = self.width
        h = self.height

        self.g.point([0, 0])
        self.g.point([w, 0])
        self.g.point([w, h])
        self.g.point([0, h])

        self.g.spline([0, 1], marker=self.bottom_id)
        self.g.spline([1, 2], marker=self.right_id)
        self.g.spline([2, 3], marker=self.top_id)
        self.g.spline([3, 0], marker=self.left_id)

        self.g.surface([0,1,2,3])

        return self.g

class RectangleWithHole(Rectangle):
    def __init__(self, width=1.0, height=1.0, element_type = 3, dofs_per_node = 1, max_area = -1):
        super().__init__(self, width, height, element_type, dofs_per_node, max_area)

class ShapeMesh:
    """Mesh generator for shapes"""
    def __init__(self, shape):
        """Initialise mesh generator"""
        self.shape = shape

        self.create()

    def create(self):
        meshGen = cfm.GmshMeshGenerator(self.shape.geometry())
        meshGen.el_type = self.shape.element_type
        meshGen.el_size_factor = self.shape.max_area
        meshGen.dofs_per_node = self.shape.dofs_per_node

        self.coords, self.edof, self.dofs, self.bdofs, self.markers = meshGen.create()

        # --- Ber├ñkna element koordinater

        self.ex, self.ey = cfc.coordxtr(self.edof, self.coords, self.dofs)
        self.pointDofs = self.dofs[list(self.shape.g.points.keys()),:]

if __name__ == "__main__":


    print("Creating rectangle")
    rect = Rectangle(5.0, 1.0, element_type=3, dofs_per_node=2, max_area=0.05)

    print("Creating mesh...")

    mesh = ShapeMesh(rect)

    print(mesh.edof)
    print(mesh.dofs)
    print(mesh.bdofs[rect.left_id])
    print(mesh.bdofs[rect.right_id])
    print(mesh.bdofs[rect.top_id])
    print(mesh.bdofs[rect.bottom_id])
    print(mesh.pointDofs)
