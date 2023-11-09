import gmsh
import sys
import numpy as np
import calfem.mesh as cfm
import calfem.vis_mpl as cfv

if __name__ == "__main__":

    gmsh.initialize(sys.argv)

    gmsh.model.add("t1")
    gmsh.model.geo.add_point(0.0, 0.0, 0.0)
    gmsh.model.geo.add_point(1.0, 0.0, 0.0)
    gmsh.model.geo.add_point(1.0, 1.0, 0.0)
    gmsh.model.geo.add_point(0.0, 1.0, 0.0)

    gmsh.model.geo.add_line(1, 2)
    gmsh.model.geo.add_line(2, 3)
    gmsh.model.geo.add_line(3, 4)
    gmsh.model.geo.add_line(4, 1)

    gmsh.model.geo.add_curve_loop([1, 2, 3, 4], 1)
    gmsh.model.geo.add_plane_surface([1], 1)

    #gmsh.model.geo.add_surface_loop([1, 2, 3, 4])

    gmsh.model.geo.synchronize()

    #gmsh.option.setNumber("Mesh.ElementOrder", 5)
    # gmsh.option.setNumber("Mesh.HighOrderOptimize", 2)
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    #gmsh.option.setNumber('Mesh.MeshSizeMin', 0.025)
    #gmsh.option.setNumber('Mesh.MeshSizeMax', 0.025)


    gmsh.model.mesh.generate(2)

    print_entities()

    #gmsh.write("t1.msh")

    gmsh.fltk.run()

    gmsh.finalize()

