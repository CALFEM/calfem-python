'''
Handle reading and writing of geometry and generated mesh from the program

Kept for backwards compatibility, use calfem.io instead.
'''

from calfem.io import (
    loadGeometry,
    saveGeometry,
    loadMesh,
    saveMesh,
    saveArrays,
    loadArrays,
    saveMatlabArrays,
)
