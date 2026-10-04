__version__ = '3.6.16'
VERSION = __version__

import sys as _sys

if _sys.platform == "win32":
    # Workaround for a crash (0xC06D007F) on Windows with MKL based NumPy
    # (e.g. conda-forge). MKL loads its OpenMP runtime at the first BLAS call.
    # If Gmsh or VTK have been loaded before that and an incompatible
    # libwinpthread-1.dll (e.g. from Git for Windows) is found on PATH, the
    # first matrix multiplication crashes. Doing a small matrix multiplication
    # here, before calfem.mesh or calfem.vis_vedo are imported, avoids this.
    try:
        import numpy as _np

        _a = _np.ones((8, 8))
        _a @ _a
        del _np, _a
    except Exception:
        pass

del _sys
