import gmsh
import sys

if __name__ == "__main__":

    gmsh.initialize(sys.argv)
    print(sys.argv)
    gmsh.finalize()

