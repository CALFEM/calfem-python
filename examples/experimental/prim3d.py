import calfem.geometry as cfg

import math

class GeomBase:
    def __init__(self):
        self.id = -1
        self.marker = -1

    def enumerate(self, next_id):
        self.id = next_id
        return next_id + 1


class Node(GeomBase):
    def __init__(self, x=0.0, y=0.0, z=0.0):
        super().__init__()

        self.x = x
        self.y = y
        self.z = z

    def is_close_to(self, n):
        return math.isclose(self.x, n.x) and \
               math.isclose(self.y, n.y) and \
               math.isclose(self.z, n.z)

    def __str__(self):
        return "Node(%10.4f, %10.4f, %10.4f, %d, %d)" % (self.x, self.y, self.z, self.id, self.marker)

class Curve(GeomBase):
    def __init__(self, p0=None, p1=None, p2=None, elements=10):
        super().__init__()
        self.element = elements
        self.__nodes = []

        if p0!=None:
            self.add_node(p0)
        if p1!=None:
            self.add_node(p1)
        if p2!=None:
            self.add_node(p2)

    def add_node(self, p):
        self.__nodes.append(p)

    def clear(self):
        self.__nodes.clear()

    def is_same_as(self, c):
        if len(self.__nodes) == len(c.nodes):
            
            count = 0
            for n in self.__nodes:
                for cn in c.nodes:
                    if n == cn:
                        count +=1
            return count == 2
        else:
            return False

    def __str__(self):
        return "Curve(" + ",".join(str(n.id) for n in self.__nodes) + ", id="+str(self.id) + ")"

    @property
    def nodes(self):
        return self.__nodes

    @property
    def node_count(self):
        return len(self.nodes)

class Surface(GeomBase):
    def __init__(self, curves=[]):
        super().__init__()
        self.__curves = curves

    def add_curve(self, c):
        self.__curves.append(c)

    def add_curves(self, lc):
        for c in lc:
            self.add_curve(c)

    def clear(self):
        self.__curves.clear()

    def __str__(self):
        return "Surface(" + ",".join(str(c.id) for c in self.__curves) + ", id="+str(self.id) + ")"

    @property
    def curves(self):
        return self.__curves

    @property
    def curve_count(self):
        return len(self.nodes)

class Volume(GeomBase):
    def __init__(self, surfaces = []):
        super().__init__()
        self.__surfaces = surfaces

    def add_surface(self, s):
        self.__surfaces.append(s)

    def clear(self):
        self.__surfaces.clear()

    def __str__(self):
        return "Volume(" + ",".join(str(s.id) for s in self.__surfaces) + ", id="+str(self.id) + ")"

    @property
    def surfaces(self):
        return self.__surfaces

    @property
    def surface_count(self):
        return len(self.__surfaces)

class GeometryModel:
    def __init__(self):
        self.__nodes = []  
        self.__curves = []
        self.__surfaces = []
        self.__volumes = []

    def find_shared_nodes(self):

        n1_dict = {}
        n2_dict = {}
        
        for n1 in self.__nodes:
            for n2 in self.__nodes:
                if n1!=n2:
                    if n1.is_close_to(n2):
                        if (n1 not in n1_dict) and (n1 not in n2_dict):
                            n1_dict[n1] = n2
                            n2_dict[n2] = n1
 
        #shared_node_list = list(set(shared_node_list))

        print("Found", len(n1_dict), "shared nodes.")

        print("----")
        for n in n1_dict.keys():
            print(n.id, "-", n1_dict[n].id)
        print("----")
        for n in n2_dict.keys():
            print(n.id, "-", n2_dict[n].id)

        return n1_dict

    def find_shared_curves(self):
        c1_dict = {}
        c2_dict = {}
        
        for c1 in self.__curves:
            for c2 in self.__curves:
                if c1!=c2:
                    if c1.is_same_as(c2):
                        if (c1 not in c1_dict) and (c1 not in c2_dict):
                            c1_dict[c1] = c2
                            c2_dict[c2] = c1
 
        print("Found", len(c1_dict), "shared curves.")

        print("----")
        for c in c1_dict.keys():
            print(c.id, "-", c1_dict[c].id)
        print("----")
        for c in c2_dict.keys():
            print(c.id, "-", c2_dict[c].id)

        return c1_dict

    def remove_duplicates(self):

        print("Before:")
        print("points = ", len(self.__nodes))
        print("curves = ", len(self.__curves))
        print("surfaces = ", len(self.__surfaces))

        self.__nodes = list(set(self.__nodes))
        self.__curves = list(set(self.__curves))
        self.__surfaces = list(set(self.__surfaces))

        print("After:")
        print("points = ", len(self.__nodes))
        print("curves = ", len(self.__curves))
        print("surfaces = ", len(self.__surfaces)) 

    
    def enumerate_nodes(self, start_count=0):
        node_count = start_count

        for node in self.__nodes:
            node_count = node.enumerate(node_count)

    def enumerate_curves(self, start_count=0):
        curve_count = start_count

        for curve in self.__curves:
            curve_count = curve.enumerate(curve_count)
    
    def enumerate_surfaces(self, start_count=0):
        surface_count = start_count

        for surface in self.__surfaces:
            surface_count = surface.enumerate(surface_count)

    def add_volume(self, volume):

        self.__volumes.append(volume)

        for surface in volume.surfaces:
            self.__surfaces.append(surface)
            for curve in surface.curves:
                self.__curves.append(curve)
                for node in curve.nodes:
                    self.__nodes.append(node)

        print("Removing duplicates...")        

        self.remove_duplicates()

        print("Enumerating nodes")

        self.enumerate_nodes()

        print("Looking for shared nodes...")

        shared_node_dict = self.find_shared_nodes()

        print("Replacing shared nodes in node list and curves...")

        for n in shared_node_dict.keys():
            replace_node = shared_node_dict[n]

            for i, nn in enumerate(self.__nodes):
                if nn == replace_node:
                    self.__nodes[i] = replace_node

            for i, c in enumerate(self.__curves):
                for j, nn in enumerate(c.nodes):
                    if nn == replace_node:
                        c.nodes[j] = replace_node

        print("Looking for shared curves...")

        shared_curve_dict = self.find_shared_curves()

        self.enumerate_nodes()
        self.enumerate_curves()
        self.enumerate_surfaces()

    @property
    def geometry(self):

        g = cfg.Geometry()

        for n in self.nodes:
            g.point([n.x, n.y, n.z])

        for i, c in enumerate(self.curves):
            print(i, c)
            g.line([n.id for n in c.nodes], ID=i, el_on_curve=10)

        for i, s in enumerate(self.surfaces):
            print(i, s)
            g.structured_surface([c.id for c in s.curves], ID=i)

        for i, v in enumerate(self.volumes):
            print(i, v)
            g.structured_volume([s.id for s in v.surfaces], ID=i)
        
        return g

    @property
    def nodes(self):
        return self.__nodes

    @nodes.setter
    def nodes(self, nodes):
        self.__nodes = nodes

    @property
    def curves(self):
        return self.__curves

    @curves.setter
    def curves(self, curves):
        self.__curves = curves

    @property
    def surfaces(self):
        return self.__surfaces

    @surfaces.setter
    def surfaces(self, surfaces):
        self.__surfaces = surfaces

    @property
    def volumes(self):
        return self.__volumes

    @volumes.setter
    def volumes(self, volumes):
        self.__volumes = volumes


class Box:
    def __init__(self, x=0.0, y=0.0, z=0.0, width=1.0, height=1.0, depth=1.0, r=10, c=10, s=10):
        self.x = x
        self.y = y
        self.z = z
        self.width = width
        self.height = height
        self.depth = depth
        self.rows = r
        self.cols = c
        self.stacks = s

        self.top_surface = 100
        self.bottom_surface = 101
        self.front_surface = 102
        self.back_surface = 103
        self.left_surface = 104
        self.right_surface = 105

        self.create()
        self.create_geom()

    def create_geom(self):

        w = self.width
        h = self.height
        d = self.depth

        r = self.rows
        c = self.cols
        s = self.stacks

        x = self.x
        y = self.y
        z = self.z

        self.p0 = Node(x, y, z) # 0
        self.p1 = Node(x+w,   y, z) # 1
        self.p2 = Node(x+w,   y+h,   z) # 2
        self.p3 = Node(x, y+h,   z) # 3

        self.p4 = Node(x, y, z+d  ) # 4
        self.p5 = Node(x+w,   y, z+d  ) # 5
        self.p6 = Node(x+w,   y+h,   z+d  ) # 6
        self.p7 = Node(x, y+h,   z+d  ) # 7

        self.c0 = Curve(self.p0, self.p1)
        self.c1 = Curve(self.p1, self.p2)
        self.c2 = Curve(self.p2, self.p3)
        self.c3 = Curve(self.p3, self.p0)

        self.c4 = Curve(self.p4, self.p5)
        self.c5 = Curve(self.p5, self.p6)
        self.c6 = Curve(self.p6, self.p7)
        self.c7 = Curve(self.p7, self.p4)

        self.c8 = Curve(self.p0, self.p4)
        self.c9 = Curve(self.p1, self.p5)
        self.c10 = Curve(self.p2, self.p6)
        self.c11 = Curve(self.p3, self.p7)

        self.s0 = Surface([self.c0, self.c1, self.c2, self.c3])
        self.s1 = Surface([self.c4, self.c5, self.c6, self.c7])
        self.s2 = Surface([self.c0, self.c9, self.c4, self.c8])
        self.s3 = Surface([self.c1, self.c10, self.c5, self.c9])
        self.s4 = Surface([self.c2, self.c10, self.c6, self.c11])
        self.s5 = Surface([self.c3, self.c11, self.c7, self.c8])

        self.volume = Volume([self.s0, self.s1, self.s2, self.s3, self.s4, self.s5])


    def create(self):

        self.geometry = cfg.geometry()
        g = self.geometry

        # Add Points

        w = self.width
        h = self.height
        d = self.depth

        r = self.rows
        c = self.cols
        s = self.stacks

        #
        #       7   (6)     6
        #        o---------o
        #    (7)/|        /|
        #      / |    (5)/ |
        #   4 /  |(4)   /5 |
        #    o---------o   | 
        #    | 3 o-----|---o 2
        #    |  /   (2)|  /
        #    | /(3)    | /(1)
        #    |/        |/
        #    o---------o
        #   0    (0)    1
        #
        #       7           8
        #        o---------o
        #       /|        /|
        #      / | (11)  / |(10)
        #     /  |      /  |
        #  4 o---------o 5 | 
        #    | 3 o-----|---o 2
        # (8)|  /   (9)|  /
        #    | /       | /
        #    |/        |/
        #    o---------o
        #   0            1
        #
        #

        g.point([0.0, 0.0, 0.0]) # 0
        g.point([w,   0.0, 0.0]) # 1
        g.point([w,   h,   0.0]) # 2
        g.point([0.0, h,   0.0]) # 3

        g.point([0.0, 0.0, d  ]) # 4
        g.point([w,   0.0, d  ]) # 5
        g.point([w,   h,   d  ]) # 6
        g.point([0.0, h,   d  ]) # 7

        # Add splines

        g.line([0, 1], 0, el_on_curve=self.rows)
        g.line([1, 2], 1, el_on_curve=self.cols)
        g.line([2, 3], 2, el_on_curve=self.rows)
        g.line([3, 0], 3, el_on_curve=self.cols)

        g.line([4, 5], 4, el_on_curve=self.rows, marker=90)
        g.line([5, 6], 5, el_on_curve=self.cols, marker=91)
        g.line([6, 7], 6, el_on_curve=self.rows, marker=92)
        g.line([7, 4], 7, el_on_curve=self.cols, marker=93)

        g.line([0, 4], 8, el_on_curve=s)
        g.line([1, 5], 9, el_on_curve=s)
        g.line([2, 6], 10, el_on_curve=s)
        g.line([3, 7], 11, el_on_curve=s)

        # Add surfaces

        g.structured_surface([0, 1,  2, 3 ], 0, marker = self.bottom_surface)
        g.structured_surface([4, 5,  6, 7  ], 1, marker = self.top_surface)
        g.structured_surface([0, 9,  4, 8 ], 2, marker = self.front_surface)
        g.structured_surface([1, 10, 5, 9 ], 3, marker = self.right_surface)
        g.structured_surface([2, 10, 6, 11], 4, marker = self.back_surface)
        g.structured_surface([3, 11, 7, 8 ], 5, marker = self.left_surface)

        # Add Volume:
        #  addStructuredVolume() takes three args. The first is a list of surface IDs 
        #  (structured surfaces). The surfaces should make a hexahedron 
        #  (i.e. 6 surfaces). Other kinds of structured volumes than hexahedra will
        #  not work for hexahedral elements, which is the only type of 3D element that 
        #  CALFEM handles. The two optional parameters are the volume ID and 
        #  volume marker.

        g.structured_volume([0,1,2,3,4,5], 0, marker=90)


if __name__ == "__main__":

    
    box1 = Box(0.0, 0.0, 0.0)
    box2 = Box(1.0, 0.0, 0.0)

    model = GeometryModel()
    model.add_volume(box1.volume)
    model.add_volume(box2.volume)

    g = model.geometry