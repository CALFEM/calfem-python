# -*- coding: iso-8859-15 -*-
"""
CALFEM Core module

Contains all the functions implementing CALFEM standard functionality.

Copyright (c) Division of Structural Mechanics and
Division of Solid Mechanics, Lund University.
"""

from __future__ import annotations
from typing import Optional, Tuple, Union

from scipy.sparse.linalg import dsolve
from scipy.sparse import csc_matrix, csr_matrix, linalg, lil_matrix
from scipy.linalg import eig, lu
import numpy as np
from numpy.typing import ArrayLike, NDArray

from calfem.matrix_compat import MatrixCompat

import logging as cflog
import sys

def disable_friendly_errors():
    """
    Disable friendly error logging and restore the previous exception hook.
    """
    sys.excepthook = __prev_exception_hook

def enable_friendly_errors():
    """
    Enable friendly error logging by setting a custom exception hook.
    """
    def __friendly_exception_hook(exc_type, exc_value, exc_traceback):
        if issubclass(exc_type, KeyboardInterrupt):
            sys.__excepthook__(exc_type, exc_value, exc_traceback)
            return
        cflog.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))
        sys.exit(1)

    global __prev_exception_hook
    __prev_exception_hook = sys.excepthook
    sys.excepthook = __friendly_exception_hook


easy_on = enable_friendly_errors
easy_off = disable_friendly_errors


def check_list_array(v, error_string):
    """
    Check if the input is a list or numpy array, raise TypeError if not.

    Parameters
    ----------
    v : object
        Object to check.
    error_string : str
        Error message to display if check fails.
    """
    fname = sys._getframe(1).f_code.co_name
    if (type(v) != list) and (type(v) != np.ndarray):
        raise TypeError("%s (%s)" % (error_string, fname))


def check_length(v, length, error_string):
    """
    Check if the input has the specified length, raise ValueError if not.

    Parameters
    ----------
    v : object
        Object to check (must support len()).
    """

    fname = sys._getframe(1).f_code.co_name

    if len(v) != length:
        raise ValueError("%s (%s)" % (error_string, fname))

def user_warning(msg: str) -> None:
    """
    Print a user warning message.

    Parameters
    ----------

    msg : str
        Warning message to display.
    """
    fname = sys._getframe(1).f_code.co_name
    print("Warning: %s (%s)" % (msg, fname))

def error(msg: str) -> None:
    """
    Log an error message.

    Parameters
    ----------

    msg : str
        Error message to log.
    """
    cflog.error(" calfem.core: "+msg)

def info(msg: str) -> None:
    """
    Log an informational message.

    Parameters
    ----------
    msg : str
        Informational message to log.
    """
    cflog.info(" calfem.core: "+msg)        
    
def spring1e(ep: ArrayLike) -> NDArray[np.floating]:
    """
    Compute the element stiffness matrix for a spring element.

    Parameters
    ----------
    ep : array_like
        Spring stiffness or analog quantity [k].

    Returns
    -------
    Ke : ndarray
        Spring stiffness matrix, shape (2, 2).

    Examples
    --------
    >>> spring1e(100)
    array([[ 100, -100],
           [-100,  100]])

    History
    -------
    LAST MODIFIED: P-E Austrell 1994-11-02
                   O Dahlblom   2022-11-15 (Python version)
    """
    k = ep  

    Ke = k * np.array([
        [1, -1],
        [-1, 1]
    ])

    return Ke


def spring1s(ep: ArrayLike, ed: ArrayLike) -> float:
    """
    Compute the element force in a spring element.

    Parameters
    ----------
    ep : array_like
        Spring stiffness or analog quantity [k].
    ed : array_like
        Element displacement vector [u1, u2].

    Returns
    -------
    es : float
        Element force.

    Examples
    --------
    >>> spring1s(100, [0.1, 0.2])
    10.0

    History
    -------
    LAST MODIFIED: P-E Austrell 1994-11-02
                   O Dahlblom   2022-11-14 (Python version)
    """
    k = ep

    N = k*(ed[1]-ed[0])
    es = N

    return es


def bar1e(ex: ArrayLike, ep: ArrayLike, eq: Optional[ArrayLike] = None) -> Union[NDArray[np.floating], Tuple[NDArray[np.floating], NDArray[np.floating]]]:
    """
    Compute the stiffness matrix (and optionally load vector) for a 1D bar element.

    Parameters
    ----------
    ex : array_like
        Element node coordinates [x1, x2].
    ep : array_like
        Element properties [E, A], where E is Young's modulus and A is cross-sectional area.
    eq : array_like, optional
        Distributed load [qX].

    Returns
    -------
    Ke : ndarray
        Bar stiffness matrix, shape (2, 2).
    fe : ndarray, optional
        Element load vector, shape (2, 1), if eq is not None.

    Examples
    --------
    >>> bar1e([0, 2], [210e9, 0.01])
    array([[ 1.05e+09, -1.05e+09],
           [-1.05e+09,  1.05e+09]])

    History
    -------
    LAST MODIFIED: O Dahlblom   2015-10-22
                   O Dahlblom   2022-11-14 (Python version)
    """
    E, A = ep
    DEA=E*A

    qX=0.
    if not eq is None:
        qX=eq[0]

    x1, x2 = ex
    dx = x2-x1
    L = abs(dx)
    
    Ke = DEA/L*np.array([
        [1, -1],
        [-1, 1]
    ])
    
    fe = qX*L*np.array([1/2, 1/2]).reshape(2,1)
  
    if eq is None:
        return Ke
    else:
        return Ke, fe


def bar1s(ex, ep, ed, eq=None, nep=None):
    """
    Compute section forces in a 1D bar element.

    Parameters
    ----------
    ex : array_like
        Element node coordinates [x1, x2].
    ep : array_like
        Element properties [E, A], where E is Young's modulus and A is cross-sectional area.
    ed : array_like
        Element displacement vector [u1, u2].
    eq : array_like, optional
        Distributed load [qX].
    nep : int, optional
        Number of evaluation points (default is 2).

    Returns
    -------
    es : ndarray
        Section forces at evaluation points, shape (nep, 1).
    edi : ndarray, optional
        Element displacements at evaluation points, shape (nep, 1), if nep is given.
    eci : ndarray, optional
        Evaluation points on the local x-axis, shape (nep, 1), if nep is given.

    Examples
    --------
    >>> bar1s([0, 2], [210e9, 0.01], [0.0, 0.001])
    array([[1.05e+06],
           [1.05e+06]])

    History
    -------
    LAST MODIFIED: O Dahlblom  2021-02-25
                   O Dahlblom  2022-11-14 (Python version)
    """
    E, A = ep
    DEA=E*A
  
    qX=0.
    if not eq is None:  
       qX=eq[0] 
     
    ne=2
    if nep != None: 
       ne=nep
     
    x1, x2 = ex
    dx = x2-x1
    L = abs(dx)
       
    a1 = ed.reshape(2,1)

    C1 = np.array([
        [1.,      0.],
        [-1/L,   1/L]
    ]) 
    
    C1a = C1 @ a1
    
    X = np.arange(0., L+L/(ne-1), L/(ne-1)).reshape(ne,1) 
    zero = np.zeros(ne).reshape(ne,1)    
    one = np.ones(ne).reshape(ne,1)
  
    u = np.concatenate((one,  X), 1) @ C1a
    du = np.concatenate((zero,  one), 1) @ C1a
  
    if DEA != 0:
       u = u -(X**2-L*X)*qX/(2*DEA)
       du = du -(2*X-L)*qX/(2*DEA)
 
    N = DEA*du
    es = N
    edi=u
    eci=X

    if nep is None:
        return es
    else:
        return es, edi, eci


def bar1we(ex, ep, eq=None):
    """
    Compute the stiffness matrix (and optionally load vector) for a 1D bar element with axial springs.

    Parameters
    ----------
    ex : array_like
        Element node coordinates [x1, x2].
    ep : array_like
        Element properties [E, A, kX], where E is Young's modulus, A is cross-sectional area, and kX is axial spring stiffness.
    eq : array_like, optional
        Distributed load [qX].

    Returns
    -------
    Ke : ndarray
        Bar stiffness matrix, shape (2, 2).
    fe : ndarray, optional
        Element load vector, shape (2, 1), if eq is not None.

    History
    -------
    LAST MODIFIED: O Dahlblom   2015-12-17
                   O Dahlblom   2022-10-19 (Python version)
    """
    E, A, kX = ep
    DEA = E*A;

    qX = 0.
    if not eq is None:
        qX = eq[0]

    x1, x2 = ex
    dx = x2-x1
    L = abs(dx)

    K1 = DEA/L*np.array([
        [1, -1],
        [-1, 1]
    ])

    K2 = kX*L/6*np.array([
        [2,  1],
        [1,  2]
    ])
   
    Ke = K1+K2

    fe = qX*L*np.array([1/2, 1/2]).reshape(2,1)
  
    if eq is None:
        return Ke
    else:
        return Ke, fe


def bar1ws(ex, ep, ed, eq=None, nep=None):
    """
    Compute section forces in a 1D bar element with axial springs.

    Parameters
    ----------
    ex : array_like
        Element node coordinates [x1, x2].
    ep : array_like
        Element properties [E, A, kX], where E is Young's modulus, A is cross-sectional area, and kX is axial spring stiffness.
    ed : array_like
        Element displacement vector [u1, u2].
    eq : array_like, optional
        Distributed load [qX].
    nep : int, optional
        Number of evaluation points (default is 2).

    Returns
    -------
    es : ndarray
        Section forces at evaluation points, shape (nep, 1).
    edi : ndarray, optional
        Element displacements at evaluation points, shape (nep, 1), if nep is given.
    eci : ndarray, optional
        Evaluation points on the local x-axis, shape (nep, 1), if nep is given.

    History
    -------
    LAST MODIFIED: O Dahlblom  2021-02-25
                   O Dahlblom  2022-11-14 (Python version)
    """
    E, A, kX = ep
    DEA = E*A
  
    qX = 0.
    if not eq is None:  
       qX = eq[0] 
     
    ne = 2
    if nep != None: 
       ne = nep
     
    x1, x2 = ex
    dx = x2-x1
    L = abs(dx)
       
    a1 = ed.reshape(2,1)

    C1 = np.array([
        [1.,      0.],
        [-1/L,   1/L]
    ]) 
    
    C1a = C1 @ a1
    
    X = np.arange(0., L+L/(ne-1), L/(ne-1)).reshape(ne,1) 
    zero = np.zeros(ne).reshape(ne,1)    
    one = np.ones(ne).reshape(ne,1)
  
    u = np.concatenate((one,  X), 1) @ C1a
    du = np.concatenate((zero,  one), 1) @ C1a
  
    if DEA != 0:
       u = u +kX/DEA*np.concatenate(((X**2-L*X)/2, (X**3-L**2*X)/6),1) @ C1a-(X**2-L*X)*qX/(2*DEA)
       du = du +kX/DEA*np.concatenate(((2*X-L)/2, (3*X**2-L**2)/6),1) @ C1a-(2*X-L)*qX/(2*DEA)
 
    N = DEA*du
    es = N
    edi = u
    eci = X

    if nep is None:
        return es
    else:
        return es, edi, eci


def bar2e(ex, ey, ep, eq=None):
    """
    Compute the element stiffness matrix (and optionally load vector) for a 2D bar element.

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2].
    ey : array_like
        Element node y-coordinates [y1, y2].
    ep : array_like
        Element properties [E, A], where E is Young's modulus and A is cross-sectional area.
    eq : array_like, optional
        Distributed load [qX].

    Returns
    -------
    Ke : ndarray
        Bar stiffness matrix, shape (4, 4).
    fe : ndarray, optional
        Element load vector, shape (4, 1), if eq is not None.

    History
    -------
    LAST MODIFIED: O Dahlblom   2015-10-20
                   O Dahlblom   2022-11-16 (Python version)
    """
    E, A = ep
    DEA = E*A

    qX = 0.
    if not eq is None:
        qX=eq[0]

    x1, x2 = ex
    y1, y2 = ey
    dx = x2-x1
    dy = y2-y1
    L = np.sqrt(dx*dx+dy*dy)
   
    Kle = DEA/L*np.array([
        [1, -1],
        [-1, 1]
    ])

    fle = qX*L*np.array([1/2, 1/2]).reshape(2,1)
    
    nxX=dx/L
    nyX=dy/L
    G = np.array([
        [nxX, nyX,   0,   0],  
        [  0,   0, nxX, nyX]
    ])
        
    Ke = G.T @ Kle @ G   
    fe = G.T @ fle

    if eq is None:
        return Ke
    else:
        return Ke, fe


def bar2s(ex, ey, ep, ed, eq=None, nep=None):
    """
    Compute normal force in a 2D bar element.

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2].
    ey : array_like
        Element node y-coordinates [y1, y2].
    ep : array_like
        Element properties [E, A], where E is Young's modulus and A is cross-sectional area.
    ed : array_like
        Element displacement vector [u1, ..., u4].
    eq : array_like, optional
        Distributed load [qX].
    nep : int, optional
        Number of evaluation points (default is 2).

    Returns
    -------
    es : ndarray
        Section forces at evaluation points, shape (nep, 1).
    edi : ndarray, optional
        Element displacements at evaluation points, shape (nep, 1), if nep is given.
    eci : ndarray, optional
        Evaluation points on the local x-axis, shape (nep, 1), if nep is given.

    History
    -------
    LAST MODIFIED: O Dahlblom  2015-12-04
                   O Dahlblom  2022-11-16 (Python version)
    """
    E, A = ep
    DEA = E*A
  
    qX = 0.
    if not eq is None:  
       qX = eq[0] 
     
    ne = 2
    if nep != None: 
       ne=nep

    x1, x2 = ex
    y1, y2 = ey
    dx = x2-x1
    dy = y2-y1
    L = np.sqrt(dx*dx+dy*dy)

    nxX = dx/L
    nyX = dy/L

    G = np.array([
        [nxX, nyX,   0,   0],  
        [  0,   0, nxX, nyX]
    ])
   
    a1 = G @ ed.reshape(4,1)

    C1 = np.array([
        [1.,      0.],
        [-1/L,   1/L]
    ]) 
    
    C1a = C1 @ a1

    X = np.arange(0., L+L/(ne-1), L/(ne-1)).reshape(ne,1) 
    zero = np.zeros(ne).reshape(ne,1)    
    one = np.ones(ne).reshape(ne,1)
  
    u = np.concatenate((one,  X), 1) @ C1a
    du = np.concatenate((zero,  one), 1) @ C1a
  
    if DEA != 0:
       u = u -(X**2-L*X)*qX/(2*DEA)
       du = du -(2*X-L)*qX/(2*DEA)
 
    N = DEA*du
    es = N
    edi = u
    eci = X

    if nep is None:
        return es
    else:
        return es, edi, eci


def bar2ge(ex, ey, ep, QX):
    """
    Compute the element stiffness matrix for a 2D bar element with additional axial force QX.

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2].
    ey : array_like
        Element node y-coordinates [y1, y2].
    ep : array_like
        Element properties [E, A], where E is Young's modulus and A is cross-sectional area.
    QX : float
        Additional axial force.

    Returns
    -------
    Ke : ndarray
        Bar stiffness matrix, shape (4, 4).

    History
    -------
    LAST MODIFIED: O Dahlblom   2015-12-17
                   O Dahlblom   2022-11-16 (Python version)
    """
    E, A = ep
    DEA = E*A

    x1, x2 = ex
    y1, y2 = ey
    dx = x2-x1
    dy = y2-y1
    L = np.sqrt(dx*dx+dy*dy)
   
    K0le = DEA/L*np.array([
        [ 1,  0, -1,  0],
        [ 0,  0,  0,  0],
        [-1,  0,  1,  0],
        [ 0,  0,  0,  0]
    ])          

    Ksle = QX/L*np.array([
        [ 0,  0,  0,  0],
        [ 0,  1,  0, -1],
        [ 0,  0,  0,  0],
        [ 0, -1,  0,  1]
    ])

    Kle = K0le + Ksle

    nxX = dx/L
    nyX = dy/L
    nxY = -dy/L
    nyY = dx/L

    G = np.array([
        [nxX, nyX,   0,   0],
        [nxY, nyY,   0,   0],  
        [  0,   0, nxX, nyX],
        [  0,   0, nxY, nyY]
    ])

    Ke = G.T @ Kle @ G   

    return Ke


def bar2gs(ex, ey, ep, ed, nep=None):
    """
    Compute normal force and axial force in a 2D bar element (bar2ge).

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2].
    ey : array_like
        Element node y-coordinates [y1, y2].
    ep : array_like
        Element properties [E, A], where E is Young's modulus and A is cross-sectional area.
    ed : array_like
        Element displacement vector [u1, ..., u4].
    nep : int, optional
        Number of evaluation points (default is 2).

    Returns
    -------
    es : ndarray
        Section forces at evaluation points, shape (nep, 1).
    QX : float
        Axial force.
    edi : ndarray, optional
        Element displacements at evaluation points, shape (nep, 1), if nep is given.
    eci : ndarray, optional
        Evaluation points on the local x-axis, shape (nep, 1), if nep is given.

    History
    -------
    LAST MODIFIED: O Dahlblom  2015-10-20
                   O Dahlblom  2022-11-16 (Python version)
    """
    E, A = ep
    DEA = E*A
  
    ne=2
    if nep != None: 
       ne=nep

    x1, x2 = ex
    y1, y2 = ey
    dx = x2-x1
    dy = y2-y1
    L = np.sqrt(dx*dx+dy*dy)

    nxX = dx/L
    nyX = dy/L
    nxY = -dy/L
    nyY = dx/L

    G = np.array([
        [nxX, nyX,   0,   0],
        [nxY, nyY,   0,   0],  
        [  0,   0, nxX, nyX],
        [  0,   0, nxY, nyY]
    ])
   
    edl = G @ ed.reshape(4,1)
    a1 = np.array([
        edl[0],
        edl[2]
    ])
      
    C1 = np.array([
        [1.,      0.],
        [-1/L,   1/L]
    ]) 
    C1a = C1 @ a1

    X = np.arange(0., L+L/(ne-1), L/(ne-1)).reshape(ne,1) 
    zero = np.zeros(ne).reshape(ne,1)    
    one = np.ones(ne).reshape(ne,1)
  
    u = np.concatenate((one,  X), 1) @ C1a
    du = np.concatenate((zero,  one), 1) @ C1a
  
    N = DEA*du
    QX = N[0]
    es = N
    edi=u
    eci=X

    if nep is None:
        return es, QX
    else:
        return es, QX, edi, eci


def bar3e(ex, ey, ez, ep, eq=None):
    """
    Compute the element stiffness matrix (and optionally load vector) for a 3D bar element.

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2].
    ey : array_like
        Element node y-coordinates [y1, y2].
    ez : array_like
        Element node z-coordinates [z1, z2].
    ep : array_like
        Element properties [E, A], where E is Young's modulus and A is cross-sectional area.
    eq : array_like, optional
        Distributed load [qX].

    Returns
    -------
    Ke : ndarray
        Bar stiffness matrix, shape (6, 6).
    fe : ndarray, optional
        Element load vector, shape (6, 1), if eq is not None.

    History
    -------
    LAST MODIFIED: O Dahlblom   2015-10-19
                   O Dahlblom   2022-11-18 (Python version)
    """
    E, A = ep
    DEA=E*A

    qX=0.
    if not eq is None:
        qX=eq[0]

    x1, x2 = ex
    y1, y2 = ey
    z1, z2 = ez
    dx = x2-x1
    dy = y2-y1
    dz = z2-z1
    L = np.sqrt(dx*dx+dy*dy+dz*dz)
   
    Kle = DEA/L*np.array([
        [1, -1],
        [-1, 1]
    ])

    fle = qX*L*np.array([1/2, 1/2]).reshape(2,1)
    
    nxX=dx/L
    nyX=dy/L
    nzX=dz/L
    G = np.array([
        [nxX, nyX, nzX,  0,   0,   0],  
        [  0,   0,   0, nxX, nyX, nzX]
    ])
        
    Ke = G.T @ Kle @ G   
    fe = G.T @ fle

    if eq is None:
        return Ke
    else:
        return Ke, fe


def bar3s(ex, ey, ez, ep, ed, eq=None, nep=None):
    """
    Compute normal force in a 3D bar element.

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2].
    ey : array_like
        Element node y-coordinates [y1, y2].
    ez : array_like
        Element node z-coordinates [z1, z2].
    ep : array_like
        Element properties [E, A], where E is Young's modulus and A is cross-sectional area.
    ed : array_like
        Element displacement vector [u1, ..., u6].
    eq : array_like, optional
        Distributed load [qX].
    nep : int, optional
        Number of evaluation points (default is 2).

    Returns
    -------
    es : ndarray
        Section forces at evaluation points, shape (nep, 1).
    edi : ndarray, optional
        Element displacements at evaluation points, shape (nep, 1), if nep is given.
    eci : ndarray, optional
        Evaluation points on the local x-axis, shape (nep, 1), if nep is given.

    History
    -------
    LAST MODIFIED: O Dahlblom  2021-09-01
                   O Dahlblom  2022-11-18 (Python version)
    """
    E, A = ep
    DEA = E*A
  
    qX = 0.
    if not eq is None:  
       qX = eq[0] 
     
    ne = 2
    if nep != None: 
       ne = nep

    x1, x2 = ex
    y1, y2 = ey
    z1, z2 = ez
    dx = x2-x1
    dy = y2-y1
    dz = z2-z1
    L = np.sqrt(dx*dx+dy*dy+dz*dz)

    nxX = dx/L
    nyX = dy/L
    nzX  =dz/L

    G = np.array([
        [nxX, nyX, nzX,   0,   0,   0],  
        [  0,   0,   0, nxX, nyX, nzX]
    ])
   
    a1 = G @ ed.reshape(6,1)

    C1 = np.array([
        [1.,      0.],
        [-1/L,   1/L]
    ]) 
    
    C1a = C1 @ a1

    X = np.linspace(0., L+L/(ne-1), ne).reshape(ne,1) 
    #X = np.arange(0., L+L/(ne-1), L/(ne-1)).reshape(ne,1) 
    zero = np.zeros(ne).reshape(ne,1)    
    one = np.ones(ne).reshape(ne,1)
  
    u = np.concatenate((one,  X), 1) @ C1a
    du = np.concatenate((zero,  one), 1) @ C1a
  
    if DEA != 0:
       u = u -(X**2-L*X)*qX/(2*DEA)
       du = du -(2*X-L)*qX/(2*DEA)
 
    N = DEA*du
    es = N
    edi=u
    eci=X

    if nep is None:
        return es
    else:
        return es, edi, eci
 

def beam1e(ex, ep, eq=None):
    """
    Compute the stiffness matrix (and optionally load vector) for a 1D beam element.

    Parameters
    ----------
    ex : array_like
        Element node coordinates [x1, x2].
    ep : array_like
        Element properties [E, I], where E is Young's modulus and I is moment of inertia.
    eq : array_like, optional
        Distributed load [qY].

    Returns
    -------
    Ke : ndarray
        Beam stiffness matrix, shape (4, 4).
    fe : ndarray, optional
        Element load vector, shape (4, 1), if eq is not None.

    History
    -------
    LAST MODIFIED: O Dahlblom   2019-01-09
                   O Dahlblom   2022-10-25 (Python version)
    """
    E, I = ep
    DEI = E*I

    qY = 0.
    if not eq is None:
        qY = eq[0]

    x1, x2 = ex
    dx = x2-x1
    L = abs(dx)

    Ke = DEI/L**3*np.array([
        [12, 6*L, -12, 6*L],
        [6*L,  4*L**2, -6*L,  2*L**2],
        [-12, -6*L,  12, -6*L],
        [6*L,  2*L**2,  -6*L, 4*L**2]
    ])
 
    fe = qY*np.array([L/2, L**2/12, L/2, -L**2/12]).reshape(4,1)
  
    if eq is None:
        return Ke
    else:
        return Ke, fe


def beam1s(ex, ep, ed, eq=None, nep=None):
    """
    Compute section forces in a 1D beam element (beam1e).

    Parameters
    ----------
    ex : array_like
        Element node coordinates [x1, x2].
    ep : array_like
        Element properties [E, I], where E is Young's modulus and I is moment of inertia.
    ed : array_like
        Element displacement vector [u1, ..., u4].
    eq : array_like, optional
        Distributed loads [qy], local directions.
    nep : int, optional
        Number of evaluation points (default is 2).

    Returns
    -------
    es : ndarray
        Section forces at evaluation points, shape (nep, 2).
    edi : ndarray, optional
        Element displacements at evaluation points, shape (nep, 1), if nep is given.
    eci : ndarray, optional
        Evaluation points on the local x-axis, shape (nep, 1), if nep is given.

    History
    -------
    LAST MODIFIED: O Dahlblom  2021-09-01
                   O Dahlblom  2022-10-25 (Python version)
    """
 
    E, I = ep
    DEI = E*I
  
    qY=0.
  
    if not eq is None:  
       qY = eq[0] 
     
    ne = 2
    if nep != None: 
       ne = nep
     
    x1, x2 = ex
    dx = x2-x1
    L = abs(dx)
       
    a2 = ed.reshape(4,1)
    
    C2 = np.array([
        [1.,      0.,    0.,     0.],
        [0.,      1.,    0.,     0.],
        [-3/L**2, -2/L,  3/L**2, -1/L],
        [2/L**3,  1/L**2, -2/L**3, 1/L**2]
    ]) 
    
    C2a = C2 @ a2
    
    X = np.arange(0., L+L/(ne-1), L/(ne-1)).reshape(ne,1)   
    zero = np.zeros(ne).reshape(ne,1)    
    one = np.ones(ne).reshape(ne,1)
  
    v = np.concatenate((one,  X, X**2, X**3), 1) @ C2a
#   dv = np.concatenate((zero,  one, 2*X, 3*X**2), 1) @ C2a
    d2v = np.concatenate((zero, zero, 2*one, 6*X), 1) @ C2a
    d3v = np.concatenate((zero, zero, zero, 6*one), 1) @ C2a
   
    if DEI != 0:
       v = v+(X**4 - 2*L*X**3 + L**2*X**2)*qY/(24*DEI)
#      dv = dv+(2*X**3 - 3*L*X**2 + L**2*X)*qY/(12*DEI)
       d2v = d2v+(6*X**2 - 6*L*X + L**2*one)*qY/(12*DEI)
       d3v = d3v+(2*X - L*one)*qY/(2*DEI)
 
    M = DEI*d2v
    V = -DEI*d3v 
    es = np.concatenate((V, M), 1)
    edi = v
    eci = X

    if nep is None:
        return es
    else:
        return es, edi, eci
 

def beam1we(ex, ep, eq=None):
    """
    Compute the stiffness matrix (and optionally load vector) for a 1D beam element on elastic foundation.

    Parameters
    ----------
    ex : array_like
        Element node coordinates [x1, x2].
    ep : array_like
        Element properties [E, I, kY], where E is Young's modulus, I is moment of inertia, and kY is transversal foundation stiffness.
    eq : array_like, optional
        Distributed load [qY].

    Returns
    -------
    Ke : ndarray
        Beam stiffness matrix, shape (4, 4).
    fe : ndarray, optional
        Element load vector, shape (4, 1), if eq is not None.

    History
    -------
    LAST MODIFIED: O Dahlblom   2016-02-17
                   O Dahlblom   2022-10-18 (Python version)
    """
    E, I, kY =ep
    DEI = E*I;

    qY = 0
    if not eq is None:
        qY = eq[0]

    x1, x2 = ex
    dx = x2-x1
    L = abs(dx)
    
    K0 = DEI/L**3*np.array([
        [12,  6*L,    -12,  6*L],
        [6*L, 4*L**2, -6*L, 2*L**2],
        [-12, -6*L,   12,   -6*L],
        [6*L, 2*L**2, -6*L, 4*L**2]
    ])    
         
    Ks = kY*L/420*np.array([
        [156,   22*L,    54,   -13*L],
        [22*L,  4*L**2,  13*L, -3*L**2],
        [54,    13*L,    156,  -22*L],
        [-13*L, -3*L**2, -22*L, 4*L**2]
    ])
   
    Ke = K0+Ks
 
    fe = qY*np.array([L/2, L**2/12, L/2, -L**2/12]).reshape(4,1)
    if eq is None:
        return Ke
    else:
        return Ke, fe


def beam1ws(ex, ep, ed, eq=None, nep=None):
    """
    Compute section forces in a 1D beam element on elastic foundation (beam1we).

    Parameters
    ----------
    ex : array_like
        Element node coordinates [x1, x2].
    ep : array_like
        Element properties [E, I, kY], where E is Young's modulus, I is moment of inertia, and kY is transversal foundation stiffness.
    ed : array_like
        Element displacement vector [u1, ..., u4].
    eq : array_like, optional
        Distributed loads [qy], local directions.
    nep : int, optional
        Number of evaluation points (default is 2).

    Returns
    -------
    es : ndarray
        Section forces at evaluation points, shape (nep, 2).
    edi : ndarray, optional
        Element displacements at evaluation points, shape (nep, 1), if nep is given.
    eci : ndarray, optional
        Evaluation points on the local x-axis, shape (nep, 1), if nep is given.

    History
    -------
    LAST MODIFIED: O Dahlblom  2021-09-01
                   O Dahlblom  2022-10-18 (Python version)
    """
    E, I, kY = ep
    DEI = E*I
  
    qY = 0.
    
    if not eq is None:  
       qY = eq[0] 
           
    ne = 2
    if nep != None: 
        ne = nep
  
    x1, x2 = ex
    dx = x2-x1
    L = abs(dx)

    a2 = ed.reshape(4,1)  
    
    C2 = np.array([
        [1.,      0.,    0.,     0.],
        [0.,      1.,    0.,     0.],
        [-3/L**2, -2/L,  3/L**2, -1/L],
        [2/L**3,  1/L**2, -2/L**3, 1/L**2]
    ]) 

    C2a = C2 @ a2
  
    X = np.arange(0., L+L/(ne-1), L/(ne-1)).reshape(ne,1) 
    zero = np.zeros(ne).reshape(ne,1)    
    one = np.ones(ne).reshape(ne,1)

    v = np.concatenate((one,  X, X**2, X**3), 1) @ C2a
    d2v = np.concatenate((zero, zero, 2*one, 6*X), 1) @ C2a
    d3v = np.concatenate((zero, zero, zero, 6*one), 1) @ C2a

    if DEI != 0:
       v = v - kY/DEI*np.concatenate((
       (X**4 - 2*L*X**3 + L**2*X**2)/24,
       (X**5 - 3*L**2*X**3 + 2*L**3*X**2)/120,
       (X**6 - 4*L**3*X**3 + 3*L**4*X**2)/360, 
       (X**7 - 5*L**4*X**3 + 4*L**5*X**2)/840), 1) @ C2a + \
       (X**4 - 2*L*X**3 + L**2*X**2)*qY/(24*DEI)
       d2v = d2v - kY/DEI*np.concatenate((
       (6*X**2 - 6*L*X + L**2*one)/12, 
       (10*X**3 - 9*L**2*X + 2*L**3*one)/60, 
       (5*X**4 - 4*L**3*X + L**4*one)/60, 
       (21*X**5 - 15*L**4*X + 4*L**5*one)/420),1) @ C2a + \
       (6*X**2 - 6*L*X + L**2*one)*qY/(12*DEI)
       d3v = d3v - kY/DEI*np.concatenate((
       (2*X - L*one)/2,
       (10*X**2 - 3*L**2)/20, 
       (5*X**3 - L**3*one)/15,
       (7*X**4 - L**4*one)/28), 1) @ C2a + \
       (2*X - L*one)*qY/(2*DEI)
    M = DEI*d2v
    V = -DEI*d3v 
    es = np.concatenate((V, M), 1)
    edi = v
    eci = X

    if nep is None:
        return es
    else:
        return es, edi, eci
 

def beam2e(ex: ArrayLike, ey: ArrayLike, ep: ArrayLike, eq: Optional[ArrayLike] = None) -> Union[NDArray[np.floating], Tuple[NDArray[np.floating], NDArray[np.floating]]]:
    """
    Compute the stiffness matrix (and optionally load vector) for a 2D beam element.

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2].
    ey : array_like
        Element node y-coordinates [y1, y2].
    ep : array_like
        Element properties [E, A, I], where E is Young's modulus, A is cross-sectional area, and I is moment of inertia.
    eq : array_like, optional
        Distributed loads [qX, qY], local directions.

    Returns
    -------
    Ke : ndarray
        Element stiffness matrix, shape (6, 6).
    fe : ndarray, optional
        Element load vector, shape (6, 1), if eq is not None.

    History
    -------
    LAST MODIFIED: O Dahlblom   2015-08-17
                   O Dahlblom   2022-11-21 (Python version)
    """
    E, A, I = ep
    DEA = E*A
    DEI = E*I
  
    qX = 0.
    qY = 0.
    if not eq is None:
        qX, qY = eq

    x1, x2 = ex
    y1, y2 = ey
    dx = x2-x1
    dy = y2-y1
    L = np.sqrt(dx*dx+dy*dy)
   
    Kle = np.array([
        [ DEA/L,           0.,          0., -DEA/L,           0.,          0.],
        [    0.,  12*DEI/L**3,  6*DEI/L**2,     0., -12*DEI/L**3,  6*DEI/L**2],
        [    0.,   6*DEI/L**2,     4*DEI/L,     0.,  -6*DEI/L**2,     2*DEI/L],
        [-DEA/L,           0.,          0.,  DEA/L,           0.,          0.],
        [    0., -12*DEI/L**3, -6*DEI/L**2,     0.,  12*DEI/L**3, -6*DEI/L**2],
        [    0.,   6*DEI/L**2,     2*DEI/L,     0.,  -6*DEI/L**2,     4*DEI/L]
    ])

    fle = L*np.array([qX/2, qY/2, qY*L/12, qX/2, qY/2, -qY*L/12]).reshape(6,1)

    nxX = dx/L
    nyX = dy/L
    nxY = -dy/L
    nyY = dx/L
    G = np.array([
        [nxX, nyX,   0,   0,   0,   0],
        [nxY, nyY,   0,   0,   0,   0],
        [  0,   0,   1,   0,   0,   0],
        [  0,   0,   0, nxX, nyX,   0],
        [  0,   0,   0, nxY, nyY,   0],
        [  0,   0,   0,   0,   0,   1]
    ])

    Ke = G.T @ Kle @ G
    fe = G.T @ fle

    if eq is None:
        return Ke
    else:
        return Ke, fe


def beam2s(ex, ey, ep, ed, eq=None, nep=None):
    """
    Compute section forces in a 2D beam element (beam2e).

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2].
    ey : array_like
        Element node y-coordinates [y1, y2].
    ep : array_like
        Element properties [E, A, I], where E is Young's modulus, A is cross-sectional area, and I is moment of inertia.
    ed : array_like
        Element displacement vector [u1, ..., u6].
    eq : array_like, optional
        Distributed loads [qX, qY], local directions.
    nep : int, optional
        Number of evaluation points (default is 2).

    Returns
    -------
    es : ndarray
        Section forces at evaluation points, shape (nep, 3).
    edi : ndarray, optional
        Element displacements at evaluation points, shape (nep, 2), if nep is given.
    eci : ndarray, optional
        Evaluation points on the local x-axis, shape (nep, 1), if nep is given.

    History
    -------
    LAST MODIFIED: O Dahlblom   2015-08-17
                   O Dahlblom   2022-11-21 (Python version)
    """
    E, A, I = ep
    DEA = E*A
    DEI = E*I
  
    qX = 0.
    qY = 0.
    if not eq is None:
        qX, qY = eq

    ne=2
    if nep != None: 
       ne=nep
    
    x1, x2 = ex
    y1, y2 = ey
    dx = x2-x1
    dy = y2-y1
    L = np.sqrt(dx*dx+dy*dy)

    nxX = dx/L
    nyX = dy/L
    nxY = -dy/L
    nyY = dx/L
    G = np.array([
        [nxX, nyX,   0,   0,   0,   0],
        [nxY, nyY,   0,   0,   0,   0],
        [  0,   0,   1,   0,   0,   0],
        [  0,   0,   0, nxX, nyX,   0],
        [  0,   0,   0, nxY, nyY,   0],
        [  0,   0,   0,   0,   0,   1]
    ])

    edl = G @ ed.reshape(6,1)

    a1 = np.array([
        edl[0],
        edl[3]
    ])
    C1 = np.array([
        [1.,      0.],
        [-1/L,   1/L]
    ]) 
    C1a = C1 @ a1

    a2 = np.array([
        edl[1],
        edl[2],
        edl[4],
        edl[5]
    ])
    C2 = np.array([
        [1.,      0.,    0.,     0.],
        [0.,      1.,    0.,     0.],
        [-3/L**2, -2/L,  3/L**2, -1/L],
        [2/L**3,  1/L**2, -2/L**3, 1/L**2]
    ]) 
    C2a = C2 @ a2

    X = np.arange(0., L+L/(ne-1), L/(ne-1)).reshape(ne,1) 
    zero = np.zeros(ne).reshape(ne,1)    
    one = np.ones(ne).reshape(ne,1)
  
    u = np.concatenate((one,  X), 1) @ C1a
    du = np.concatenate((zero,  one), 1) @ C1a
    if DEA != 0:
       u = u -(X**2-L*X)*qX/(2*DEA)
       du = du -(2*X-L)*qX/(2*DEA)

    v = np.concatenate((one,  X, X**2, X**3), 1) @ C2a
#   dv = np.concatenate((zero,  one, 2*X, 3*X**2), 1) @ C2a
    d2v=np.concatenate((zero, zero, 2*one, 6*X), 1) @ C2a
    d3v = np.concatenate((zero, zero, zero, 6*one), 1) @ C2a
    if DEI != 0:
       v = v+(X**4 - 2*L*X**3 + L**2*X**2)*qY/(24*DEI)
#      dv = dv+(2*X**3 - 3*L*X**2 + L**2*X)*qY/(12*DEI)
       d2v = d2v+(6*X**2 - 6*L*X + L**2*one)*qY/(12*DEI)
       d3v = d3v+(2*X - L*one)*qY/(2*DEI)
 
    N = DEA*du
    M = DEI*d2v
    V = -DEI*d3v 
    es = np.concatenate((N, V, M), 1)
    edi = np.concatenate((u, v), 1)
    eci = X

    if nep is None:
        return es
    else:
        return es, edi, eci
  

def beam2we(ex, ey, ep, eq=None):
    """
    Compute the stiffness matrix (and optionally load vector) for a 2D beam element on elastic foundation.

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2].
    ey : array_like
        Element node y-coordinates [y1, y2].
    ep : array_like
        Element properties [E, A, I, kX, kY], where E is Young's modulus, A is cross-sectional area, I is moment of inertia, kX is axial foundation stiffness, and kY is transversal foundation stiffness.
    eq : array_like, optional
        Distributed loads [qX, qY], local directions.

    Returns
    -------
    Ke : ndarray
        Element stiffness matrix, shape (6, 6).
    fe : ndarray, optional
        Element load vector, shape (6, 1), if eq is not None.

    History
    -------
    LAST MODIFIED: O Dahlblom   2015-08-07
                   O Dahlblom   2022-11-21 (Python version)
    """
    E, A, I, kX, kY = ep
    DEA = E*A
    DEI = E*I

    qX = 0
    qY = 0
    if not eq is None:
        qX, qY = eq

    x1, x2 = ex
    y1, y2 = ey
    dx = x2-x1
    dy = y2-y1
    L = np.sqrt(dx*dx+dy*dy)
   
    K0 = np.array([
        [ DEA/L,           0.,          0., -DEA/L,           0.,          0.],
        [    0.,  12*DEI/L**3,  6*DEI/L**2,     0., -12*DEI/L**3,  6*DEI/L**2],
        [    0.,   6*DEI/L**2,     4*DEI/L,     0.,  -6*DEI/L**2,     2*DEI/L],
        [-DEA/L,           0.,          0.,  DEA/L,           0.,          0.],
        [    0., -12*DEI/L**3, -6*DEI/L**2,     0.,  12*DEI/L**3, -6*DEI/L**2],
        [     0.,  6*DEI/L**2.,     2*DEI/L,     0.,  -6*DEI/L**2,     4*DEI/L]
    ])
  
    Ks = L/420*np.array([
        [140*kX, 0,       0,          70*kX,  0,       0],
        [0,      156*kY,  22*kY*L,    0,      54*kY,  -13*kY*L],
        [0,      22*kY*L, 4*kY*L**2,  0,      13*kY*L, -3*kY*L**2],
        [70*kX,  0,       0,          140*kX, 0,       0],
        [0,      54*kY,   13*kY*L,    0,      156*kY, -22*kY*L],
        [0,     -13*kY*L, -3*kY*L**2, 0,     -22*kY*L, 4*kY*L**2]
    ])
   
    Kle = K0+Ks
    
    fle = L*np.array([qX/2, qY/2, qY*L/12, qX/2, qY/2, -qY*L/12]).reshape(6,1)

    nxX = dx/L
    nyX = dy/L
    nxY = -dy/L
    nyY = dx/L
    G = np.array([
        [nxX, nyX,   0,   0,   0,   0],
        [nxY, nyY,   0,   0,   0,   0],
        [  0,   0,   1,   0,   0,   0],
        [  0,   0,   0, nxX, nyX,   0],
        [  0,   0,   0, nxY, nyY,   0],
        [  0,   0,   0,   0,   0,   1]
    ])

    Ke = G.T @ Kle @ G
    fe = G.T @ fle

    if eq is None:
        return Ke
    else:
        return Ke, fe


def beam2ws(ex, ey, ep, ed, eq=None, nep=None):
    """
    Compute section forces in a 2D beam element on elastic foundation.

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2].
    ey : array_like
        Element node y-coordinates [y1, y2].
    ep : array_like
        Element properties [E, A, I, kX, kY], where E is Young's modulus, A is cross-sectional area, I is moment of inertia, kX is axial foundation stiffness, and kY is transversal foundation stiffness.
    ed : array_like
        Element displacement vector [u1, ..., u6].
    eq : array_like, optional
        Distributed loads [qX, qY], local directions.
    nep : int, optional
        Number of evaluation points (default is 2).

    Returns
    -------
    es : ndarray
        Section forces at evaluation points, shape (nep, 3).
    edi : ndarray, optional
        Element displacements at evaluation points, shape (nep, 2), if nep is given.
    eci : ndarray, optional
        Evaluation points on the local x-axis, shape (nep, 1), if nep is given.

    History
    -------
    LAST MODIFIED: O Dahlblom   2022-09-30
                   O Dahlblom   2022-11-21 (Python version)
    """
    E, A, I, kX, kY = ep
    DEA = E*A
    DEI = E*I

    qX = 0
    qY = 0
    if not eq is None:
        qX, qY = eq

    x1, x2 = ex
    y1, y2 = ey
    dx = x2-x1
    dy = y2-y1
    L = np.sqrt(dx*dx+dy*dy)
   
    ne = 2
    if nep != None: 
       ne = nep

    nxX = dx/L
    nyX = dy/L
    nxY = -dy/L
    nyY = dx/L
    G = np.array([
        [nxX, nyX,   0,   0,   0,   0],
        [nxY, nyY,   0,   0,   0,   0],
        [  0,   0,   1,   0,   0,   0],
        [  0,   0,   0, nxX, nyX,   0],
        [  0,   0,   0, nxY, nyY,   0],
        [  0,   0,   0,   0,   0,   1]
    ])

    edl = G @ ed.reshape(6,1)

    a1 = np.array([
        edl[0],
        edl[3]
    ])
    C1 = np.array([
        [1.,      0.],
        [-1/L,   1/L]
    ]) 
    C1a = C1 @ a1

    a2 = np.array([
        edl[1],
        edl[2],
        edl[4],
        edl[5]
    ])
    C2 = np.array([
        [1.,      0.,    0.,     0.],
        [0.,      1.,    0.,     0.],
        [-3/L**2, -2/L,  3/L**2, -1/L],
        [2/L**3,  1/L**2, -2/L**3, 1/L**2]
    ]) 
    C2a = C2 @ a2

    X = np.arange(0., L+L/(ne-1), L/(ne-1)).reshape(ne,1) 
    zero = np.zeros(ne).reshape(ne,1)    
    one = np.ones(ne).reshape(ne,1)
  
    u = np.concatenate((one,  X), 1) @ C1a
    du = np.concatenate((zero,  one), 1) @ C1a
    if DEA != 0:
       u = u +kX/DEA*np.concatenate(((X**2-L*X)/2, (X**3-L**2*X)/6),1) @ C1a-(X**2-L*X)*qX/(2*DEA)
       du = du +kX/DEA*np.concatenate(((2*X-L)/2, (3*X**2-L**2)/6),1) @ C1a-(2*X-L)*qX/(2*DEA)

    v = np.concatenate((one,  X, X**2, X**3), 1) @ C2a
#   dv = np.concatenate((zero,  one, 2*X, 3*X**2), 1) @ C2a
    d2v = np.concatenate((zero, zero, 2*one, 6*X), 1) @ C2a
    d3v = np.concatenate((zero, zero, zero, 6*one), 1) @ C2a
    if DEI != 0:
       v = v - kY/DEI*np.concatenate((
       (X**4 - 2*L*X**3 + L**2*X**2)/24,
       (X**5 - 3*L**2*X**3 + 2*L**3*X**2)/120,
       (X**6 - 4*L**3*X**3 + 3*L**4*X**2)/360, 
       (X**7 - 5*L**4*X**3 + 4*L**5*X**2)/840), 1) @ C2a + \
       (X**4 - 2*L*X**3 + L**2*X**2)*qY/(24*DEI)
       d2v = d2v - kY/DEI*np.concatenate((
       (6*X**2 - 6*L*X + L**2*one)/12, 
       (10*X**3 - 9*L**2*X + 2*L**3*one)/60, 
       (5*X**4 - 4*L**3*X + L**4*one)/60, 
       (21*X**5 - 15*L**4*X + 4*L**5*one)/420),1) @ C2a + \
       (6*X**2 - 6*L*X + L**2*one)*qY/(12*DEI)
       d3v = d3v - kY/DEI*np.concatenate((
       (2*X - L*one)/2,
       (10*X**2 - 3*L**2)/20, 
       (5*X**3 - L**3*one)/15,
       (7*X**4 - L**4*one)/28), 1) @ C2a + \
       (2*X - L*one)*qY/(2*DEI)

    N = DEA*du
    M = DEI*d2v
    V = -DEI*d3v 
    es = np.concatenate((N, V, M), 1)
    edi = np.concatenate((u, v), 1)
    eci = X

    if nep is None:
        return es
    else:
        return es, edi, eci
  

def beam2ge(ex, ey, ep, QX, eq=None):
    """
    Compute the element stiffness matrix (and optionally load vector) for a 2D beam element with geometric nonlinearity.

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2].
    ey : array_like
        Element node y-coordinates [y1, y2].
    ep : array_like
        Element properties [E, A, I], where E is Young's modulus, A is cross-sectional area, and I is moment of inertia.
    QX : float
        Axial force in the beam.
    eq : array_like, optional
        Distributed transverse load [qY].

    Returns
    -------
    Ke : ndarray
        Element stiffness matrix, shape (6, 6).
    fe : ndarray, optional
        Element load vector, shape (6, 1), if eq is not None.

    History
    -------
    LAST MODIFIED: O Dahlblom   2015-12-17
                   O Dahlblom   2022-12-08 (Python version)
    """
    E, A, I = ep
    DEA = E*A
    DEI = E*I

    if eq != None:
        if np.size(eq) > 1:
            error("eq should be a scalar !!!")
            return
        else:
            qY = eq[0]
    else:
        qY = 0

    x1, x2 = ex
    y1, y2 = ey
    dx = x2-x1
    dy = y2-y1
    L = np.sqrt(dx*dx+dy*dy)

    K0le = np.array([
        [ DEA/L,           0.,          0., -DEA/L,           0.,          0.],
        [    0.,  12*DEI/L**3,  6*DEI/L**2,     0., -12*DEI/L**3,  6*DEI/L**2],
        [    0.,   6*DEI/L**2,     4*DEI/L,     0.,  -6*DEI/L**2,     2*DEI/L],
        [-DEA/L,           0.,          0.,  DEA/L,           0.,          0.],
        [    0., -12*DEI/L**3, -6*DEI/L**2,     0.,  12*DEI/L**3, -6*DEI/L**2],
        [    0.,   6*DEI/L**2,     2*DEI/L,     0.,  -6*DEI/L**2,     4*DEI/L]
    ])

    Ksle = QX/(30*L)*np.array([
        [   0.,   0.,     0.,   0.,   0.,     0.],
        [   0.,  36.,    3*L,   0., -36.,    3*L],
        [   0.,  3*L, 4*L**2,   0., -3*L,  -L**2],
        [   0.,   0.,     0.,   0.,   0.,     0.],
        [   0., -36.,   -3*L,   0.,  36.,   -3*L],
        [   0.,  3*L,  -L**2,   0., -3*L, 4*L**2]
    ])
    
    fle = qY*L*np.array([0, 1/2, L/12, 0, 1/2, -L/12]).reshape(6,1)

    nxX = dx/L
    nyX = dy/L
    nxY = -dy/L
    nyY = dx/L
    G = np.array([
        [nxX, nyX,   0,   0,   0,   0],
        [nxY, nyY,   0,   0,   0,   0],
        [  0,   0,   1,   0,   0,   0],
        [  0,   0,   0, nxX, nyX,   0],
        [  0,   0,   0, nxY, nyY,   0],
        [  0,   0,   0,   0,   0,   1]
    ])

    Kle = K0le+Ksle 
    Ke = G.T @ Kle @ G
    fe = G.T @ fle

    if eq is None:
        return Ke
    else:
        return Ke, fe


def beam2gs(ex, ey, ep, ed, QX, eq=None, nep=None):
    """
    Calculate section forces in a 2D nonlinear beam element (beam2ge).

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2].
    ey : array_like
        Element node y-coordinates [y1, y2].
    ep : array_like
        Element properties [E, A, I], where E is Young's modulus, A is cross-sectional area, and I is moment of inertia.
    ed : array_like
        Element displacement vector [u1, ..., u6].
    QX : float
        Axial force in the beam.
    eq : array_like, optional
        Distributed transverse load [qY].
    nep : int, optional
        Number of evaluation points (default is 2).

    Returns
    -------
    es : ndarray
        Section forces at evaluation points, shape (nep, 3).
    QX : float
        Axial force.
    edi : ndarray, optional
        Element displacements at evaluation points, shape (nep, 2), if nep is given.
    eci : ndarray, optional
        Evaluation points on the local x-axis, shape (nep, 1), if nep is given.

    History
    -------
    LAST MODIFIED: O Dahlblom   2021-09-01
                   O Dahlblom   2022-12-06 (Python version)
    """
    E, A, I = ep
    DEA = E*A
    DEI = E*I

    if eq != None:
        if np.size(eq) > 1:
            error("eq should be a scalar !!!")
            return
        else:
            qY = eq[0]
    else:
        qY = 0

    x1, x2 = ex
    y1, y2 = ey
    dx = x2-x1
    dy = y2-y1
    L = np.sqrt(dx*dx+dy*dy)

    ne = 2
    if nep != None: 
       ne = nep

    nxX = dx/L
    nyX = dy/L
    nxY = -dy/L
    nyY = dx/L
    G = np.array([
        [nxX, nyX,   0,   0,   0,   0],
        [nxY, nyY,   0,   0,   0,   0],
        [  0,   0,   1,   0,   0,   0],
        [  0,   0,   0, nxX, nyX,   0],
        [  0,   0,   0, nxY, nyY,   0],
        [  0,   0,   0,   0,   0,   1]
    ])

    edl = G @ ed.reshape(6,1)

    X = np.arange(0., L+L/(ne-1), L/(ne-1)).reshape(ne,1) 
    zero = np.zeros(ne).reshape(ne,1)    
    one = np.ones(ne).reshape(ne,1)
  
    a1 = np.array([
        edl[0],
        edl[3]
    ])
    C1 = np.array([
        [1.,      0.],
        [-1/L,   1/L]
    ]) 
    C1a = C1 @ a1

    a2 = np.array([
        edl[1],
        edl[2],
        edl[4],
        edl[5]
    ])
    C2 = np.array([
        [1.,      0.,    0.,     0.],
        [0.,      1.,    0.,     0.],
        [-3/L**2, -2/L,  3/L**2, -1/L],
        [2/L**3,  1/L**2, -2/L**3, 1/L**2]
    ]) 
    C2a = C2 @ a2

    X = np.arange(0., L+L/(ne-1), L/(ne-1)).reshape(ne,1) 
    zero = np.zeros(ne).reshape(ne,1)    
    one = np.ones(ne).reshape(ne,1)
  
    u = np.concatenate((one,  X), 1) @ C1a
    du = np.concatenate((zero,  one), 1) @ C1a
    v = np.concatenate((one,  X, X**2, X**3), 1) @ C2a
    dv = np.concatenate((zero,  one, 2*X, 3*X**2), 1) @ C2a
    d2v = np.concatenate((zero, zero, 2*one, 6*X), 1) @ C2a
    d3v = np.concatenate((zero, zero, zero, 6*one), 1) @ C2a
    if DEI != 0:
       v = v + QX/DEI*np.concatenate((zero, zero, (X**4-2*L*X**3+L**2*X**2)/12, (X**5-3*L**2*X**3+2*L**3*X**2)/20), 1) @ C2a + \
       (X**4 - 2*L*X**3 + L**2*X**2)*qY/(24*DEI)
       dv = dv + QX/DEI*np.concatenate((zero, zero, (2*X**3-3*L*X**2+L**2*X)/6, (5*X**4-9*L**2*X**2+4*L**3*X)/20), 1) @ C2a + \
       (2*X**3 - 3*L*X**2 + L**2*X)*qY/(12*DEI)
       d2v = d2v + QX/DEI*np.concatenate((zero, zero, (6*X**2-6*L*X+L**2*one)/6, (10*X**3-9*L**2*X+2*L**3*one)/10), 1) @ C2a + \
       (6*X**2 - 6*L*X + L**2*one)*qY/(12*DEI)
       d3v = d3v + QX/DEI*np.concatenate((zero, zero, (2*X-L*one), (30*X**2-9*L**2*one)/10), 1) @ C2a + \
       (2*X - L*one)*qY/(2*DEI)

    QX = DEA*du.item(0)
    M = DEI*d2v
    V = -DEI*d3v
    N=QX+dv*V
    es = np.concatenate((N, V, M), 1)
    edi = np.concatenate((u, v), 1)
    eci = X
    
    if nep is None:
        return es, QX
    else:
        return es, QX, edi, eci


def beam2gxe(ex, ey, ep, QX, eq=None):
    """
    Compute the element stiffness matrix (and optionally load vector) for a 2D beam element with geometric nonlinearity (exact solution).

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2].
    ey : array_like
        Element node y-coordinates [y1, y2].
    ep : array_like
        Element properties [E, A, I], where E is Young's modulus, A is cross-sectional area, and I is moment of inertia.
    QX : float
        Axial force in the beam.
    eq : array_like, optional
        Distributed transverse load.

    Returns
    -------
    Ke : ndarray
        Element stiffness matrix, shape (6, 6).
    fe : ndarray, optional
        Element load vector, shape (6, 1), if eq is not None.

    History
    -------
    LAST MODIFIED: O Dahlblom   2021-06-21
                   O Dahlblom   2022-12-06 (Python version)
    """
    E, A, I = ep
    DEA = E*A
    DEI = E*I

    if eq != None:
        if np.size(eq) > 1:
            error("eq should be a scalar !!!")
            return
        else:
            qY = eq[0]
    else:
        qY = 0

    x1, x2 = ex
    y1, y2 = ey
    dx = x2-x1
    dy = y2-y1
    L = np.sqrt(dx*dx+dy*dy)

    eps = 1e-12

    if QX < -eps*DEI/L**2: 
        kL = np.sqrt(-QX/DEI)*L
        f1 = (kL/2)/np.tan(kL/2)
        f2 = kL**2/(12*(1-f1))
        f3 = f1/4+3*f2/4
        f4 = -f1/2+3*f2/2
        f5 = f1*f2
        h = 6*(2/kL**2-(1+np.cos(kL))/(kL*np.sin(kL)))
    elif QX > eps*DEI/L**2:
        kL = np.sqrt(QX/DEI)*L
        f1 = (kL/2)/np.tanh(kL/2)
        f2 = -(1/12.)*kL**2/(1-f1)
        f3 = f1/4+3*f2/4
        f4 = -f1/2+3*f2/2
        f5 = f1*f2
        h = -6*(2/kL**2-(1+np.cosh(kL))/(kL*np.sinh(kL)))
    else:
        f1 = f2 = f3 = f4 = f5 = h = 1

    Kle = np.array([
         [DEA/L,              0.,             0., -DEA/L, 0.,                          0.],
        [    0.,  12*DEI*f5/L**3,  6*DEI*f2/L**2,     0., -12*DEI*f5/L**3,  6*DEI*f2/L**2],
        [    0.,   6*DEI*f2/L**2,     4*DEI*f3/L,     0.,  -6*DEI*f2/L**2,     2*DEI*f4/L],
        [-DEA/L,              0.,             0.,  DEA/L,              0.,             0.],
        [    0., -12*DEI*f5/L**3, -6*DEI*f2/L**2,     0.,  12*DEI*f5/L**3, -6*DEI*f2/L**2],
            [0.,   6*DEI*f2/L**2,     2*DEI*f4/L,     0.,  -6*DEI*f2/L**2,     4*DEI*f3/L]
    ])

    fle = qY*L*np.array([0., 1/2., L*h/12, 0., 1/2., -L*h/12]).reshape(6,1)

    nxX = dx/L
    nyX = dy/L
    nxY = -dy/L
    nyY = dx/L
    G = np.array([
        [nxX, nyX,   0,   0,   0,   0],
        [nxY, nyY,   0,   0,   0,   0],
        [  0,   0,   1,   0,   0,   0],
        [  0,   0,   0, nxX, nyX,   0],
        [  0,   0,   0, nxY, nyY,   0],
        [  0,   0,   0,   0,   0,   1]
    ])

    Ke = G.T @ Kle @ G
    fe = G.T @ fle

    if eq is None:
        return Ke
    else:
        return Ke, fe


def beam2gxs(ex, ey, ep, ed, QX, eq=None, nep=None):
    """

    Calculate section forces in a 2D nonlinear beam element (beam2gxe).

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2].
    ey : array_like
        Element node y-coordinates [y1, y2].
    ep : array_like
        Element properties [E, A, I], where E is Young's modulus, A is cross-sectional area, and I is moment of inertia.
    ed : array_like
        Element displacement vector [u1, ..., u6].
    QX : float
        Axial force in the beam.
    eq : array_like, optional
        Distributed transverse load.
    nep : int, optional
        Number of evaluation points (default is 2).

    Returns
    -------
    es : ndarray
        Section forces at evaluation points, shape (nep, 3).
    QX : float
        Axial force.
    edi : ndarray, optional
        Element displacements at evaluation points, shape (nep, 2), if nep is given.
    eci : ndarray, optional
        Evaluation points on the local x-axis, shape (nep, 1), if nep is given.

    History
    -------
    LAST MODIFIED: O Dahlblom   2021-06-21
                   O Dahlblom   2022-12-06 (Python version)
    """
    E, A, I = ep
    DEA = E*A
    DEI = E*I

    if eq != None:
        if np.size(eq) > 1:
            error("eq should be a scalar !!!")
            return
        else:
            qY = eq[0]
    else:
        qY = 0

    x1, x2 = ex
    y1, y2 = ey
    dx = x2-x1
    dy = y2-y1
    L = np.sqrt(dx*dx+dy*dy)

    ne = 2
    if nep != None: 
       ne = nep

    nxX = dx/L
    nyX = dy/L
    nxY = -dy/L
    nyY = dx/L
    G = np.array([
        [nxX, nyX,   0,   0,   0,   0],
        [nxY, nyY,   0,   0,   0,   0],
        [  0,   0,   1,   0,   0,   0],
        [  0,   0,   0, nxX, nyX,   0],
        [  0,   0,   0, nxY, nyY,   0],
        [  0,   0,   0,   0,   0,   1]
    ])

    edl = G @ ed.reshape(6,1)

    X = np.arange(0., L+L/(ne-1), L/(ne-1)).reshape(ne,1) 
    zero = np.zeros(ne).reshape(ne,1)    
    one = np.ones(ne).reshape(ne,1)
  
    a1 = np.array([
        edl[0],
        edl[3]
    ])
    C1 = np.array([
        [1.,      0.],
        [-1/L,   1/L]
    ]) 
    C1a = C1 @ a1
    u = np.concatenate((one,  X), 1) @ C1a
    du = np.concatenate((zero,  one), 1) @ C1a

    a2 = np.array([
        edl[1],
        edl[2],
        edl[4],
        edl[5]
    ])

    eps = 1e-12

    if QX < -eps*DEI/L**2:
        k = np.sqrt(-QX/DEI)
        kL = k*L
        C2 = 1/(k*(-2*(1-np.cos(kL))+kL*np.sin(kL)))*np.array([
            [k*(kL*np.sin(kL)+np.cos(kL)-1), -kL*np.cos(kL)+np.sin(kL), -k*(1-np.cos(kL)), -np.sin(kL)+kL],
            [-(k**2)*np.sin(kL), -k*(1-np.cos(kL)), (k**2)*np.sin(kL), -k*(1-np.cos(kL))],
            [-k*(1-np.cos(kL)), kL*np.cos(kL)-np.sin(kL),  k*(1-np.cos(kL)), np.sin(kL)-kL],
            [k*np.sin(kL), (kL*np.sin(kL)+np.cos(kL)-1), -k*np.sin(kL), (1-np.cos(kL))]
        ]) 

        C2a = C2 @ a2
        v = np.concatenate((one,  X, np.cos(k*X), np.sin(k*X)), 1) @ C2a
        dv = np.concatenate((zero, one, -k*np.sin(k*X), k*np.cos(k*X)), 1) @ C2a
        d2v = np.concatenate((zero, zero, -k**2*np.cos(k*X), -k**2*np.sin(k*X)), 1) @ C2a
        d3v = np.concatenate((zero, zero, k**3*np.sin(k*X), -k**3*np.cos(k*X)), 1) @ C2a
        if DEI != 0:
          v = v+qY*L**4/(2*DEI)*((1+np.cos(kL))/(kL**3*np.sin(kL))*(-1+np.cos(k*X))+np.sin(k*X)/kL**3+X*(-1+X/L)/(kL**2*L))
          dv = dv+qY*L**3/(2*DEI)*((1+np.cos(kL))/(kL**2*np.sin(kL))*(-np.sin(k*X))+np.cos(k*X)/kL**2+(-1+2*X/L)/kL**2)
          d2v = d2v+qY*L**2/(2*DEI)*((1+np.cos(kL))/(kL*np.sin(kL))*(-np.cos(k*X))-np.sin(k*X)/kL+2/(kL**2)) 
          d3v = d3v+qY*L/(2*DEI)*((1+np.cos(kL))/np.sin(kL)*(np.sin(k*X))-np.cos(k*X))
    elif QX > eps*DEI/L**2:
        k = np.sqrt(QX/DEI)
        kL = k*L
        C2 = 1/(k*(-2*(1-np.cosh(kL))-kL*np.sinh(kL)))*np.array([
            [k*(-kL*np.sinh(kL)+np.cosh(kL)-1), -kL*np.cosh(kL)+np.sinh(kL), -k*(1-np.cosh(kL)), -np.sinh(kL)+kL],
            [(k**2)*np.sinh(kL), -k*(1-np.cosh(kL)), -(k**2)*np.sinh(kL), -k*(1-np.cosh(kL))],
            [-k*(1-np.cosh(kL)), kL*np.cosh(kL)-np.sinh(kL),  k*(1-np.cosh(kL)), np.sinh(kL)-kL],
            [-k*np.sinh(kL), (-kL*np.sinh(kL)+np.cosh(kL)-1), k*np.sinh(kL), (1-np.cosh(kL))]
        ]) 
        C2a = C2 @ a2
        v = np.concatenate((one,  X, np.cosh(k*X), np.sinh(k*X)), 1) @ C2a
        dv = np.concatenate((zero, one, k*np.sinh(k*X), k*np.cosh(k*X)), 1) @ C2a
        d2v = np.concatenate((zero, zero, k**2*np.cosh(k*X), k**2*np.sinh(k*X)), 1) @ C2a
        d3v = np.concatenate((zero, zero, k**3*np.sinh(k*X), k**3*np.cosh(k*X)), 1) @ C2a
        if DEI != 0:
          v = v+qY*L**4/(2*DEI)*((1+np.cosh(kL))/(kL**3*np.sinh(kL))*(-1+np.cosh(k*X))-np.sinh(k*X)/kL**3+X*(1-X/L)/(kL**2*L))
          dv = dv+qY*L**3/(2*DEI)*((1+np.cosh(kL))/(kL**2*np.sinh(kL))*(np.sinh(k*X))-np.cosh(k*X)/kL**2+(1-2*X/L)/kL**2)
          d2v = d2v+qY*L**2/(2*DEI)*((1+np.cosh(kL))/(kL*np.sinh(kL))*(np.cosh(k*X))-np.sinh(k*X)/kL-2/(kL**2)) 
          d3v = d3v+qY*L/(2*DEI)*((1+np.cosh(kL))/np.sinh(kL)*(np.sinh(k*X))-np.cosh(k*X))
    else:
        C2 = np.array([
            [1.,      0.,    0.,     0.],
            [0.,      1.,    0.,     0.],
            [-3/L**2, -2/L,  3/L**2, -1/L],
            [2/L**3,  1/L**2, -2/L**3, 1/L**2]
        ]) 
        C2a = C2 @ a2
        v = np.concatenate((one,  X, X**2, X**3), 1) @ C2a
        dv = np.concatenate((zero,  one, 2*X, 3*X**2), 1) @ C2a
        d2v = np.concatenate((zero, zero, 2*one, 6*X), 1) @ C2a
        d3v = np.concatenate((zero, zero, zero, 6*one), 1) @ C2a
        if DEI != 0:
           v = v+(X**4 - 2*L*X**3 + L**2*X**2)*qY/(24*DEI)
           dv = dv+(2*X**3 -3*L*X**2 +L**2*X)*qY/(12*DEI)
           d2v = d2v +(6*X**2 - 6*L*X + L**2*one)*qY/(12*DEI)
           d3v = d3v +(2*X - L*one)*qY/(2*DEI)

    QX = DEA*du.item(0)
    M = DEI*d2v
    V = -DEI*d3v
    N=QX+dv*V
    es = np.concatenate((N, V, M), 1)
    edi = np.concatenate((u, v), 1)
    eci = X
    
    if nep is None:
        return es, QX
    else:
        return es, QX, edi, eci


def beam2te(ex, ey, ep, eq=None):
    """
    Compute the stiffness matrix (and optionally load vector) for a 2D Timoshenko beam element.

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2].
    ey : array_like
        Element node y-coordinates [y1, y2].
    ep : array_like
        Element properties [E, Gm, A, I, ks], where E is Young's modulus, Gm is shear modulus, 
        A is cross-sectional area, I is moment of inertia, and ks is shear correction factor.
    eq : array_like, optional
        Distributed loads in local directions [qX, qY].

    Returns
    -------
    Ke : ndarray
        Element stiffness matrix, shape (6, 6).
    fe : ndarray, optional
        Element load vector, shape (6, 1), if eq is not None.

    History
    -------
    LAST MODIFIED: O Dahlblom   2021-11-05
                   O Dahlblom   2022-12-08 (Python version)
    """
    E, Gm, A, I, ks = ep
    DEA = E*A
    DEI = E*I
    DGAks = Gm*A*ks
  
    qX = 0.
    qY = 0.
    if not eq is None:
        qX, qY = eq

    x1, x2 = ex
    y1, y2 = ey
    dx = x2-x1
    dy = y2-y1
    L = np.sqrt(dx*dx+dy*dy)
    m = (12*DEI)/(L**2*DGAks)
    f1=1/(1+m)
    f2=f1*(1+m/4)
    f3=f1*(1-m/2)


    Kle = np.array([
        [ DEA/L,              0.,             0., -DEA/L,              0.,              0.],
        [    0.,  12*DEI*f1/L**3,  6*DEI*f1/L**2,     0., -12*DEI*f1/L**3,   6*DEI*f1/L**2],
        [    0.,   6*DEI*f1/L**2,     4*DEI*f2/L,     0.,  -6*DEI*f1/L**2,      2*DEI*f3/L],
        [-DEA/L,              0.,             0.,  DEA/L,              0.,              0.],
        [    0., -12*DEI*f1/L**3, -6*DEI*f1/L**2,     0.,   12*DEI*f1/L**3, -6*DEI*f1/L**2],
        [    0.,   6*DEI*f1/L**2,     2*DEI*f3/L,     0.,   -6*DEI*f1/L**2,     4*DEI*f2/L]
    ])

    fle = L*np.array([qX/2, qY/2, qY*L/12, qX/2, qY/2, -qY*L/12]).reshape(6,1)

    nxX = dx/L
    nyX = dy/L
    nxY = -dy/L
    nyY = dx/L
    G = np.array([
        [nxX, nyX,   0,   0,   0,   0],
        [nxY, nyY,   0,   0,   0,   0],
        [  0,   0,   1,   0,   0,   0],
        [  0,   0,   0, nxX, nyX,   0],
        [  0,   0,   0, nxY, nyY,   0],
        [  0,   0,   0,   0,   0,   1]
    ])

    Ke = G.T @ Kle @ G
    fe = G.T @ fle

    if eq is None:
        return Ke
    else:
        return Ke, fe
    

def beam2ts(ex, ey, ep, ed, eq=None, nep=None):
    """
    Compute section forces in a 2D Timoshenko beam element (beam2te).

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2].
    ey : array_like
        Element node y-coordinates [y1, y2].
    ep : array_like
        Element properties [E, G, A, I, ks], where E is Young's modulus, G is shear modulus,
        A is cross-sectional area, I is moment of inertia, and ks is shear correction factor.
    ed : array_like
        Element displacement vector [u1, ..., u6].
    eq : array_like, optional
        Distributed loads in local directions [qx, qy].
    nep : int, optional
        Number of evaluation points (default is 2).

    Returns
    -------
    es : ndarray
        Section forces in local directions at evaluation points, shape (nep, 3).
        Each row contains [N, V, M] (normal force, shear force, moment).
    edi : ndarray, optional
        Element displacements in local directions at evaluation points, shape (nep, 3), if nep is given.
        Each row contains [u, v, theta]. Note: For Timoshenko beam element, the rotation of the 
        cross section is not equal to dv/dx.
    eci : ndarray, optional
        Local x-coordinates of the evaluation points, shape (nep, 1), if nep is given.

    History
    -------
    LAST MODIFIED: O Dahlblom   2021-11-05
                   O Dahlblom   2022-12-08 (Python version)
    """
    E, Gm, A, I, ks = ep
    DEA = E*A
    DEI = E*I
    DGAks = Gm*A*ks
    alpha=DEI/DGAks
  
    qX = 0.
    qY = 0.
    if not eq is None:
        qX, qY = eq

    ne=2
    if nep != None: 
       ne=nep
    
    x1, x2 = ex
    y1, y2 = ey
    dx = x2-x1
    dy = y2-y1
    L = np.sqrt(dx*dx+dy*dy)

    nxX = dx/L
    nyX = dy/L
    nxY = -dy/L
    nyY = dx/L
    G = np.array([
        [nxX, nyX,   0,   0,   0,   0],
        [nxY, nyY,   0,   0,   0,   0],
        [  0,   0,   1,   0,   0,   0],
        [  0,   0,   0, nxX, nyX,   0],
        [  0,   0,   0, nxY, nyY,   0],
        [  0,   0,   0,   0,   0,   1]
    ])

    edl = G @ ed.reshape(6,1)

    a1 = np.array([
        edl[0],
        edl[3]
    ])
    C1 = np.array([
        [1.,      0.],
        [-1/L,   1/L]
    ]) 
    C1a = C1 @ a1

    a2 = np.array([
        edl[1],
        edl[2],
        edl[4],
        edl[5]
    ])
    C2 = 1/(L**2+12*alpha)*np.array([
        [L**2+12*alpha,             0.,         0.,           0.],
        [  -12*alpha/L,   L**2+6*alpha, 12*alpha/L,     -6*alpha],
        [          -3., -2*L-6*alpha/L,         3., -L+6*alpha/L],
        [          2/L,             1.,       -2/L,           1.]
    ]) 
    C2a = C2 @ a2

    X = np.arange(0., L+L/(ne-1), L/(ne-1)).reshape(ne,1) 
    zero = np.zeros(ne).reshape(ne,1)    
    one = np.ones(ne).reshape(ne,1)
  
    u = np.concatenate((one,  X), 1) @ C1a
    du = np.concatenate((zero,  one), 1) @ C1a
    if DEA != 0:
       u = u -(X**2-L*X)*qX/(2*DEA)
       du = du -(2*X-L)*qX/(2*DEA)

    v = np.concatenate((one,  X, X**2, X**3), 1) @ C2a
    dv = np.concatenate((zero,  one, 2*X, 3*X**2), 1) @ C2a
    theta = np.concatenate((zero, one, 2*X, 3*X**2+6*alpha), 1) @ C2a
    dtheta=np.concatenate((zero, zero, 2*one, 6*X), 1) @ C2a
    if DEI != 0:
       v = v+(X**4 - 2*L*X**3 + L**2*X**2)*qY/(24*DEI)+(-X**2/2+L*X/2)*qY/DGAks
       dv = dv+(2*X**3 - 3*L*X**2 + L**2*X)*qY/(12*DEI)+(-X+L/2)*qY/DGAks
       theta = theta+(2*X**3-3*L*X**2+L**2*X)*qY/(12*DEI)
       dtheta = dtheta+(6*X**2-6*L*X+L**2)*qY/(12*DEI)

    N = DEA*du
    M = DEI*dtheta
    V = DGAks*(dv-theta) 
    es = np.concatenate((N, V, M), 1)
    edi = np.concatenate((u, v, theta), 1)
    eci = X

    if nep is None:
        return es
    else:
        return es, edi, eci
  

def beam2de(ex, ey, ep):
    """
    Calculate the stiffness matrix Ke, mass matrix Me, and optionally damping matrix Ce for a 2D elastic Bernoulli beam element.

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2].
    ey : array_like
        Element node y-coordinates [y1, y2].
    ep : array_like
        Element properties [E, A, I, m, (a, b)], where E is Young's modulus, A is cross-sectional area,
        I is moment of inertia, m is mass per unit length, and a, b are optional damping coefficients
        where Ce = a*Me + b*Ke.

    Returns
    -------
    Ke : ndarray
        Element stiffness matrix, shape (6, 6).
    Me : ndarray
        Element mass matrix, shape (6, 6).
    Ce : ndarray, optional
        Element damping matrix, shape (6, 6), if damping coefficients are provided.

    History
    -------
    LAST MODIFIED: K Persson    1995-08-23
                   O Dahlblom   2022-12-08 (Python version)
    """
    b = np.array([
        [ex[1]-ex[0]],
        [ey[1]-ey[0]]
    ])
    L = np.sqrt(b.T @ b).item(0)
    n = np.array(b/L).reshape(2,)

    a = 0
    b = 0
    if np.size(ep) == 4:
        E, A, I, m = ep
    elif np.size(ep) == 6:
        E, A, I, m, a, b = ep

    print ("E, A, I, m, a, b, L")    
    print (E, A, I, m, a, b, L) 


    Kle = np.array([
        [E*A/L, 0,           0,         -E*A/L, 0,           0],
        [0,     12*E*I/L**3, 6*E*I/L**2, 0,    -12*E*I/L**3, 6*E*I/L**2],
        [0,     6*E*I/L**2,  4*E*I/L,    0,    -6*E*I/L**2,  2*E*I/L],
        [-E*A/L, 0,           0,          E*A/L, 0,           0],
        [0,    -12*E*I/L**3, -6*E*I/L**2, 0,     12*E*I/L**3, -6*E*I/L**2],
        [0,     6*E*I/L**2,  2*E*I/L,    0,    -6*E*I/L**2,  4*E*I/L]
    ])

    Mle = m*L/420*np.array([
        [140, 0,    0,      70,  0,    0],
        [0,   156,  22*L,   0,   54,  -13*L],
        [0,   22*L, 4*L**2, 0,   13*L, -3*L**2],
        [70,  0,    0,      140, 0,    0],
        [0,   54,   13*L,   0,   156, -22*L],
        [0,  -13*L, -3*L**2, 0,  -22*L, 4*L**2]
    ])

    Cle = a*Mle+b*Kle

    G = np.array([
        [n[0], n[1], 0, 0,    0,    0],
        [-n[1], n[0], 0, 0,    0,    0],
        [0,    0,    1, 0,    0,    0],
        [0,    0,    0, n[0], n[1], 0],
        [0,    0,    0, -n[1], n[0], 0],
        [0,    0,    0, 0,    0,    1]
    ])

    Ke = G.T @ Kle @ G
    Me = G.T @ Mle @ G
    Ce = G.T @ Cle @ G

    if np.size(ep) == 4:
        return Ke, Me
    elif np.size(ep) == 6:
        return Ke, Me, Ce


def beam2ds(ex, ey, ep, ed, ev, ea):
    """
    Calculate element forces for multiple identical 2D Bernoulli beam elements in dynamic analysis.

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2].
    ey : array_like
        Element node y-coordinates [y1, y2].
    ep : array_like
        Element properties [E, A, I, m, (a, b)], where E is Young's modulus, A is cross-sectional area,
        I is moment of inertia, m is mass per unit length, and a, b are optional damping coefficients
        where Ce = a*Me + b*Ke.
    ed : array_like
        Element displacement matrix.
    ev : array_like
        Element velocity matrix.
    ea : array_like
        Element acceleration matrix.

    Returns
    -------
    es : ndarray
        Element forces in local directions, shape (nie, 6).
        Each row contains [-N1, -V1, -M1, N2, V2, M2].

    History
    -------
    LAST MODIFIED: K Persson    1995-08-23
                   O Dahlblom   2022-12-08 (Python version)
    """
    b = np.array([
        [ex[1]-ex[0]],
        [ey[1]-ey[0]]
    ])
    L = np.sqrt(b.T @ b).item(0)
    n = np.array(b/L).reshape(2,)

    a = 0
    b = 0
    if np.size(ep) == 4:
        E, A, I, m = ep
    elif np.size(ep) == 6:
        E, A, I, m, a, b = ep

    Kle = np.array([
        [E*A/L, 0,           0,         -E*A/L, 0,           0],
        [0,     12*E*I/L**3, 6*E*I/L**2, 0,    -12*E*I/L**3, 6*E*I/L**2],
        [0,     6*E*I/L**2,  4*E*I/L,    0,    -6*E*I/L**2,  2*E*I/L],
        [-E*A/L, 0,           0,          E*A/L, 0,           0],
        [0,    -12*E*I/L**3, -6*E*I/L**2, 0,     12*E*I/L**3, -6*E*I/L**2],
        [0,     6*E*I/L**2,  2*E*I/L,    0,    -6*E*I/L**2,  4*E*I/L]
    ])

    Mle = m*L/420*np.array([
        [140, 0,    0,      70,  0,    0],
        [0,   156,  22*L,   0,   54,  -13*L],
        [0,   22*L, 4*L**2, 0,   13*L, -3*L**2],
        [70,  0,    0,      140, 0,    0],
        [0,   54,   13*L,   0,   156, -22*L],
        [0,  -13*L, -3*L**2, 0,  -22*L, 4*L**2]
    ])

    Cle = a*Mle+b*Kle

    G = np.array([
        [n[0], n[1], 0, 0,    0,    0],
        [-n[1], n[0], 0, 0,    0,    0],
        [0,    0,    1, 0,    0,    0],
        [0,    0,    0, n[0], n[1], 0],
        [0,    0,    0, -n[1], n[0], 0],
        [0,    0,    0, 0,    0,    1]
    ])

    nie, ned = ed.shape
    es = np.array(np.zeros((nie, 6)))
    for i in range(nie):
       d = ed[i,:].reshape(6,1)
       v = ev[i,:].reshape(6,1)
       a = ea[i,:].reshape(6,1)
       es[i,:] = (Kle @ G @ d + Cle @ G @ v + Mle @ G @ a).T

#   [nie,ned]=size(ed);
#   for i=1:nie
#        d=ed(i,:)';
#        v=ev(i,:)';
#        a=ea(i,:)';
#        es(i,:)=(Kle*G*d+Cle*G*v+Mle*G*a)';

    return es


def beam3e(ex, ey, ez, eo, ep, eq=None):
    """
    Calculate the stiffness matrix (and optionally load vector) for a 3D elastic Bernoulli beam element.

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2].
    ey : array_like
        Element node y-coordinates [y1, y2].
    ez : array_like
        Element node z-coordinates [z1, z2].
    eo : array_like
        Orientation of local z-axis [xz, yz, zz].
    ep : array_like
        Element properties [E, G, A, Iy, Iz, Kv], where E is Young's modulus, G is shear modulus,
        A is cross-sectional area, Iy is moment of inertia about local y-axis, 
        Iz is moment of inertia about local z-axis, and Kv is Saint-Venant's torsion constant.
    eq : array_like, optional
        Distributed loads in local directions [qX, qY, qZ, qW].

    Returns
    -------
    Ke : ndarray
        Element stiffness matrix, shape (12, 12).
    fe : ndarray, optional
        Element load vector, shape (12, 1), if eq is not None.

    History
    -------
    LAST MODIFIED: O Dahlblom   2015-10-19
                   O Dahlblom   2022-11-21 (Python version)
    """
    E, Gs, A, Iy, Iz, Kv = ep
    DEA = E*A 
    DEIz = E*Iz
    DEIy = E*Iy
    DGK = Gs*Kv

    qX = 0.
    qY = 0.
    qZ = 0.
    qW = 0.
    if not eq is None:
        qX, qY, qZ, qW = eq

    x1, x2 = ex
    y1, y2 = ey
    z1, z2 = ez
    dx = x2-x1
    dy = y2-y1
    dz = z2-z1
    L = np.sqrt(dx*dx+dy*dy+dz*dz)

    a = DEA/L
    b = 12*DEIz/L**3
    c = 6*DEIz/L**2
    d = 12*DEIy/L**3
    e = 6*DEIy/L**2
    f = DGK/L
    g = 2*DEIy/L
    h = 2*DEIz/L

    Kle = np.array([
        [a, 0, 0, 0, 0,   0,  -a, 0, 0, 0, 0,   0],
        [0, b, 0, 0, 0,   c,   0, -b, 0, 0, 0,   c],
        [0, 0, d, 0, -e,   0,   0, 0, -d, 0, -e,   0],
        [0, 0, 0, f, 0,   0,   0, 0, 0, -f, 0,   0],
        [0, 0, -e, 0, 2*g, 0,   0, 0, e, 0, g,   0],
        [0, c, 0, 0, 0,   2*h, 0, -c, 0, 0, 0,   h],
        [-a, 0, 0, 0, 0,   0,   a, 0, 0, 0, 0,   0],
        [0, -b, 0, 0, 0,  -c,   0, b, 0, 0, 0,  -c],
        [0, 0, -d, 0, e,   0,   0, 0, d, 0, e,   0],
        [0, 0, 0, -f, 0,   0,   0, 0, 0, f, 0,   0],
        [0, 0, -e, 0, g,   0,   0, 0, e, 0, 2*g, 0],
        [0, c, 0, 0, 0,   h,   0, -c, 0, 0, 0,   2*h]
    ])

    fle = L/2*np.array([qX, qY, qZ, qW, -qZ*L/6, qY*L/6,
                     qX, qY, qZ, qW, qZ*L/6, -qY*L/6]).reshape(12,1)

    n1 = np.array([dx, dy, dz])/L
    
    lc = np.sqrt(eo @ eo.T)
    eo1= np.array(eo/lc)
    n2a = np.cross(eo1,n1)
    n2al = np.sqrt(n2a @ n2a.T)
    n2=n2a/n2al
    
    n3 = np.cross(n1,n2)
    
    G = np.array([
        [n1[0], n1[1], n1[2], 0,     0,     0,
            0,     0,     0,     0,     0,     0],
        [n2[0], n2[1], n2[2], 0,     0,     0,
            0,     0,     0,     0,     0,     0],
        [n3[0], n3[1], n3[2], 0,     0,     0,
            0,     0,     0,     0,     0,     0],
        [0,     0,     0,     n1[0], n1[1], n1[2],
            0,     0,     0,     0,     0,     0],
        [0,     0,     0,     n2[0], n2[1], n2[2],
            0,     0,     0,     0,     0,     0],
        [0,     0,     0,     n3[0], n3[1], n3[2],
            0,     0,     0,     0,     0,     0],
        [0,     0,     0,     0,     0,     0,
            n1[0], n1[1], n1[2], 0,     0,     0],
        [0,     0,     0,     0,     0,     0,
            n2[0], n2[1], n2[2], 0,     0,     0],
        [0,     0,     0,     0,     0,     0,
            n3[0], n3[1], n3[2], 0,     0,     0],
        [0,     0,     0,     0,     0,     0,     0,
            0,     0,     n1[0], n1[1], n1[2]],
        [0,     0,     0,     0,     0,     0,     0,
            0,     0,     n2[0], n2[1], n2[2]],
        [0,     0,     0,     0,     0,     0,
            0,     0,     0,     n3[0], n3[1], n3[2]]
    ])

    Ke = G.T @ Kle @ G
    fe = G.T @ fle

    if eq is None:
        return Ke
    else:
        return Ke, fe


def beam3s(ex, ey, ez, eo, ep, ed, eq=None, nep=None):
    """
    Calculate the variation of section forces and displacements along a 3D beam element.

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2].
    ey : array_like
        Element node y-coordinates [y1, y2].
    ez : array_like
        Element node z-coordinates [z1, z2].
    eo : array_like
        Orientation of local z-axis [xz, yz, zz].
    ep : array_like
        Element properties [E, G, A, Iy, Iz, Kv], where E is Young's modulus, G is shear modulus,
        A is cross-sectional area, Iy is moment of inertia about local y-axis, 
        Iz is moment of inertia about local z-axis, and Kv is Saint-Venant's torsion constant.
    ed : array_like
        Element displacement vector [u1, ..., u12].
    eq : array_like, optional
        Distributed loads in local directions [qX, qY, qZ, qW].
    nep : int, optional
        Number of evaluation points (default is 2).

    Returns
    -------
    es : ndarray
        Section forces at evaluation points, shape (nep, 6).
        Each row contains [N, Vy, Vz, T, My, Mz] (normal force, shear forces, torque, moments).
    edi : ndarray, optional
        Element displacements at evaluation points, shape (nep, 4), if nep is given.
        Each row contains [u, v, w, phi] (displacements and rotation).
    eci : ndarray, optional
        Local x-coordinates of the evaluation points, shape (nep, 1), if nep is given.

    History
    -------
    LAST MODIFIED: O Dahlblom   2015-10-19
                   O Dahlblom   2022-11-23 (Python version)
    """
    E, Gs, A, Iy, Iz, Kv = ep
    DEA = E*A 
    DEIz = E*Iz
    DEIy = E*Iy
    DGK = Gs*Kv

    qX = 0.
    qY = 0.
    qZ = 0.
    qW = 0.
    if not eq is None:
        qX, qY, qZ, qW = eq

    ne = 2
    if nep != None:
        ne = nep

    x1, x2 = ex
    y1, y2 = ey
    z1, z2 = ez
    dx = x2-x1
    dy = y2-y1
    dz = z2-z1
    L = np.sqrt(dx*dx+dy*dy+dz*dz)
    n1 = np.array([dx, dy, dz])/L
    
    lc = np.sqrt(eo @ eo.T)
    eo1= np.array(eo/lc)
    n2a = np.cross(eo1,n1)
    n2al = np.sqrt(n2a @ n2a.T)
    n2=n2a/n2al
    
    n3 = np.cross(n1,n2)
        
    G = np.array([
        [n1[0], n1[1], n1[2], 0,     0,     0,
            0,     0,     0,     0,     0,     0],
        [n2[0], n2[1], n2[2], 0,     0,     0,
            0,     0,     0,     0,     0,     0],
        [n3[0], n3[1], n3[2], 0,     0,     0,
            0,     0,     0,     0,     0,     0],
        [0,     0,     0,     n1[0], n1[1], n1[2],
            0,     0,     0,     0,     0,     0],
        [0,     0,     0,     n2[0], n2[1], n2[2],
            0,     0,     0,     0,     0,     0],
        [0,     0,     0,     n3[0], n3[1], n3[2],
            0,     0,     0,     0,     0,     0],
        [0,     0,     0,     0,     0,     0,
            n1[0], n1[1], n1[2], 0,     0,     0],
        [0,     0,     0,     0,     0,     0,
            n2[0], n2[1], n2[2], 0,     0,     0],
        [0,     0,     0,     0,     0,     0,
            n3[0], n3[1], n3[2], 0,     0,     0],
        [0,     0,     0,     0,     0,     0,     0,
            0,     0,     n1[0], n1[1], n1[2]],
        [0,     0,     0,     0,     0,     0,     0,
            0,     0,     n2[0], n2[1], n2[2]],
        [0,     0,     0,     0,     0,     0,
            0,     0,     0,     n3[0], n3[1], n3[2]]
    ])

    edl = G @ ed.reshape(12,1)

    a1 = np.array([
        edl[0],
        edl[6]
    ])
    C1 = np.array([
        [1.,      0.],
        [-1/L,   1/L]
    ]) 
    C1a = C1 @ a1

    a2 = np.array([
        edl[1],
        edl[5],
        edl[7],
        edl[11]
    ])
    C2 = np.array([
        [1.,      0.,    0.,     0.],
        [0.,      1.,    0.,     0.],
        [-3/L**2, -2/L,  3/L**2, -1/L],
        [2/L**3,  1/L**2, -2/L**3, 1/L**2]
    ]) 
    C2a = C2 @ a2

    a3 = np.array([
        edl[2],
        -edl[4],
        edl[8],
        -edl[10]
    ])
    C3 = C2
    C3a = C3 @ a3
  
    a4 = np.array([
        edl[3],
        edl[9]
    ])
    C4 = C1
    C4a = C4 @ a4
  
    X = np.linspace(0., L+L/(ne-1), ne).reshape(ne,1) 
    #X = np.arange(0., L+L/(ne-1), L/(ne-1)).reshape(ne,1) 
    zero = np.zeros(ne).reshape(ne,1)    
    one = np.ones(ne).reshape(ne,1)
  
    u = np.concatenate((one,  X), 1) @ C1a
    du = np.concatenate((zero,  one), 1) @ C1a
    if DEA != 0:
       u = u-(X**2-L*X)*qX/(2*DEA)
       du = du-(2*X-L)*qX/(2*DEA)

    v = np.concatenate((one,  X, X**2, X**3), 1) @ C2a
    d2v=np.concatenate((zero, zero, 2*one, 6*X), 1) @ C2a
    d3v = np.concatenate((zero, zero, zero, 6*one), 1) @ C2a
    if DEIz != 0:
       v = v+(X**4 - 2*L*X**3 + L**2*X**2)*qY/(24*DEIz)
       d2v = d2v+(6*X**2 - 6*L*X + L**2*one)*qY/(12*DEIz)
       d3v = d3v+(2*X - L*one)*qY/(2*DEIz)
 
    w = np.concatenate((one,  X, X**2, X**3), 1) @ C3a
    d2w = np.concatenate((zero, zero, 2*one, 6*X), 1) @ C3a
    d3w = np.concatenate((zero, zero, zero, 6*one), 1) @ C3a
    if DEIy != 0:
       w = w+(X**4 - 2*L*X**3 + L**2*X**2)*qZ/(24*DEIy)
       d2w = d2w+(6*X**2 - 6*L*X + L**2*one)*qZ/(12*DEIy)
       d3w = d3w+(2*X - L*one)*qZ/(2*DEIy)
 
    fi = np.concatenate((one,  X), 1) @ C4a
    dfi = np.concatenate((zero,  one), 1) @ C4a
    if DGK != 0:
       fi = fi-(X**2-L*X)*qW/(2*DGK)
       dfi = dfi-(2*X-L)*qW/(2*DGK)
    N = DEA*du
    Mz = DEIz*d2v
    Vy = -DEIz*d3v 
    My = -DEIy*d2w
    Vz = -DEIy*d3w
    T = DGK*dfi
    es = np.concatenate((N, Vy, Vz, T, My, Mz), 1)
    edi = np.concatenate((u, v, w, fi), 1)
    eci = X

    if nep is None:
        return es
    else:
        return es, edi, eci
 

def flw2te(ex, ey, ep, D, eq=None):
    """
    Compute element stiffness (conductivity) matrix for a triangular field element.

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2, x3].
    ey : array_like
        Element node y-coordinates [y1, y2, y3].
    ep : array_like
        Element properties [t], where t is the element thickness.
    D : array_like
        Constitutive matrix [[kxx, kxy], [kyx, kyy]].
    eq : float, optional
        Heat supply per unit volume.

    Returns
    -------
    Ke : ndarray
        Element 'stiffness' matrix, shape (3, 3).
    fe : ndarray
        Element load vector, shape (3, 1).
    """
    t = ep[0]
    if eq is None:
        eq = 0.

    exm = np.asmatrix(ex)
    eym = np.asmatrix(ey)
    C = np.asmatrix(np.hstack([np.ones((3, 1)), exm.T, eym.T]))
    B = np.matrix([
        [0., 1., 0.],
        [0., 0., 1.]
    ])*C.I
    A = 0.5*np.linalg.det(C)

    Ke = B.T*D*B*t*A
    fe = np.matrix([[1., 1., 1.]]).T*eq*A*t/3

    if eq == 0.:
        return Ke
    else:
        return Ke, fe


def flw2ts(ex, ey, D, ed):
    """
    Compute flows or corresponding quantities in the triangular field element.

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2, x3].
    ey : array_like
        Element node y-coordinates [y1, y2, y3].
    D : array_like
        Constitutive matrix [[kxx, kxy], [kyx, kyy]].
    ed : array_like
        Element nodal values [u1, u2, u3], one row per element.

    Returns
    -------
    es : ndarray
        Element flows, shape (n_elem, 2) or (2,) for single element.
    et : ndarray
        Element gradients, shape (n_elem, 2) or (2,) for single element.
    """

    if len(ex.shape) > 1:
        qs = np.zeros([ex.shape[0], 2])
        qt = np.zeros([ex.shape[0], 2])
        row = 0
        for exr, eyr, edr in zip(ex, ey, ed):
            exm = np.asarray(exr).reshape(-1)  # Ensure 1D array
            eym = np.asarray(eyr).reshape(-1)  # Ensure 1D array
            edm = np.asarray(edr).reshape(-1)  # Ensure 1D array
            
            # Create C matrix with consistent dimensions
            ones_col = np.ones((3, 1))
            # Reshape column vectors properly
            ex_col = exm.reshape(-1, 1)
            ey_col = eym.reshape(-1, 1)
            
            # Stack horizontally with consistent 2D shapes
            C = np.hstack([ones_col, ex_col, ey_col])
            
            B = np.array([
                [0., 1., 0.],
                [0., 0., 1.]
            ]) @ np.linalg.inv(C)

            qs[row, :] = (-D @ B @ edm.T).T
            qt[row, :] = (B @ edm.T).T
            row += 1

        return qs, qt
    else:
        exm = np.asarray(ex).reshape(-1)  # Ensure 1D array
        eym = np.asarray(ey).reshape(-1)  # Ensure 1D array
        edm = np.asarray(ed).reshape(-1)  # Ensure 1D array
        
        # Create C matrix with consistent dimensions
        ones_col = np.ones((3, 1))
        # Reshape column vectors properly
        ex_col = exm.reshape(-1, 1)
        ey_col = eym.reshape(-1, 1)
        
        # Stack horizontally with consistent 2D shapes
        C = np.hstack([ones_col, ex_col, ey_col])
        
        B = np.array([
            [0., 1., 0.],
            [0., 0., 1.]
        ]) @ np.linalg.inv(C)

        qs = -D @ B @ edm.T
        qt = B @ edm.T

        return qs.T, qt.T
    
def flw2qe(ex, ey, ep, D, eq=None):
    """
    Compute element stiffness (conductivity) matrix for a quadrilateral field element.
    This function calculates the element stiffness matrix and optionally the load vector
    for a 4-node quadrilateral element used in 2D heat flow analysis. The quadrilateral
    is divided into 4 triangular sub-elements using the centroid as a common node.

    Parameters
    ----------
    ex : array_like
        Element x-coordinates [x1, x2, x3, x4] for the 4 corner nodes.
    ey : array_like
        Element y-coordinates [y1, y2, y3, y4] for the 4 corner nodes.
    ep : array_like
        Element properties [t] where t is the element thickness.
    D : array_like
        Constitutive matrix (2x2) for heat conductivity:
        [[kxx, kxy],
         [kyx, kyy]]
        where kxx, kyy are conductivities in x and y directions,
        and kxy, kyx are cross-conductivities.
    eq : float, optional
        Heat supply per unit volume. If None, only stiffness matrix is computed.
        Default is None.

    Returns
    -------
    Ke : ndarray
        Element stiffness matrix (4x4) for the quadrilateral element.
    fe : ndarray, optional
        Element load vector (4x1). Only returned when eq is provided.
    """
    xc = sum(ex)/4.
    yc = sum(ey)/4.

    K = np.zeros((5, 5))
    f = np.zeros((5, 1))

    if eq is None:
        k1 = flw2te([ex[0], ex[1], xc], [ey[0], ey[1], yc], ep, D)
        K = assem(np.array([1, 2, 5]), K, k1)
        k1 = flw2te([ex[1], ex[2], xc], [ey[1], ey[2], yc], ep, D)
        K = assem(np.array([2, 3, 5]), K, k1)
        k1 = flw2te([ex[2], ex[3], xc], [ey[2], ey[3], yc], ep, D)
        K = assem(np.array([3, 4, 5]), K, k1)
        k1 = flw2te([ex[3], ex[0], xc], [ey[3], ey[0], yc], ep, D)
        K = assem(np.array([4, 1, 5]), K, k1)
    else:
        k1, f1 = flw2te([ex[0], ex[1], xc], [ey[0], ey[1], yc], ep, D, eq)
        K, f = assem(np.array([1, 2, 5]), K, k1, f, f1)
        k1, f1 = flw2te([ex[1], ex[2], xc], [ey[1], ey[2], yc], ep, D, eq)
        K, f = assem(np.array([2, 3, 5]), K, k1, f, f1)
        k1, f1 = flw2te([ex[2], ex[3], xc], [ey[2], ey[3], yc], ep, D, eq)
        K, f = assem(np.array([3, 4, 5]), K, k1, f, f1)
        k1, f1 = flw2te([ex[3], ex[0], xc], [ey[3], ey[0], yc], ep, D, eq)
        K, f = assem(np.array([4, 1, 5]), K, k1, f, f1)
    Ke1, fe1 = statcon(K, f, np.array([5]))

    Ke = Ke1
    fe = fe1

    if eq is None:
        return Ke
    else:
        return Ke, fe


def flw2qs(ex, ey, ep, D, ed, eq=None):
    """
    Compute flows or corresponding quantities in the quadrilateral field element.
    
    Parameters
    ----------

    ex : array_like
        Element node x-coordinates [x1, x2, x3, x4].
    ey : array_like
        Element node y-coordinates [y1, y2, y3, y4].
    ep : array_like
        Element properties [t], where t is element thickness.
    D : array_like
        Constitutive matrix [[kxx, kxy], [kyx, kyy]].
    ed : array_like
        Element nodal values [[u1, u2, u3, u4], [.., .., .., ..]], where u1,u2,u3,u4 are nodal values.
    
    Returns
    -------

    es : ndarray
        Element flows [[qx, qy], [.., ..]].
    et : ndarray
        Element gradients [[gx, gy], [.., ..]].
    """
    K = np.zeros((5, 5))
    f = np.zeros((5, 1))

    xm = sum(ex)/4
    ym = sum(ey)/4

    if eq is None:
        q = 0
    else:
        q = eq

    En = np.array([
        [1, 2, 5],
        [2, 3, 5],
        [3, 4, 5],
        [4, 1, 5]
    ])
    ex1 = np.array([ex[0], ex[1], xm])
    ey1 = np.array([ey[0], ey[1], ym])
    ex2 = np.array([ex[1], ex[2], xm])
    ey2 = np.array([ey[1], ey[2], ym])
    ex3 = np.array([ex[2], ex[3], xm])
    ey3 = np.array([ey[2], ey[3], ym])
    ex4 = np.array([ex[3], ex[0], xm])
    ey4 = np.array([ey[3], ey[0], ym])

    if eq is None:
        k1 = flw2te(ex1, ey1, ep, D)
        K = assem(En[0], K, k1)
        k1 = flw2te(ex2, ey2, ep, D)
        K = assem(En[1], K, k1)
        k1 = flw2te(ex3, ey3, ep, D)
        K = assem(En[2], K, k1)
        k1 = flw2te(ex4, ey4, ep, D)
        K = assem(En[3], K, k1)
    else:
        k1, f1 = flw2te(ex1, ey1, ep, D, q)
        K, f = assem(En[0], K, k1, f, f1)
        k1, f1 = flw2te(ex2, ey2, ep, D, q)
        K, f = assem(En[1], K, k1, f, f1)
        k1, f1 = flw2te(ex3, ey3, ep, D, q)
        K, f = assem(En[2], K, k1, f, f1)
        k1, f1 = flw2te(ex4, ey4, ep, D, q)
        K, f = assem(En[3], K, k1, f, f1)

    if ed.ndim == 1:
        ed = np.array([ed])

    ni, nj = np.shape(ed)

    a = np.zeros((5, ni))
    for i in range(ni):
        a[np.ix_(range(5), [i])], r = np.asarray(
            solveq(K, f, np.arange(1, 5), ed[i]))

    s1, t1 = flw2ts(ex1, ey1, D, a[np.ix_(En[0, :]-1, np.arange(ni))].T)
    s2, t2 = flw2ts(ex2, ey2, D, a[np.ix_(En[1, :]-1, np.arange(ni))].T)
    s3, t3 = flw2ts(ex3, ey3, D, a[np.ix_(En[2, :]-1, np.arange(ni))].T)
    s4, t4 = flw2ts(ex4, ey4, D, a[np.ix_(En[3, :]-1, np.arange(ni))].T)

    es = (s1+s2+s3+s4)/4.
    et = (t1+t2+t3+t4)/4.

    return es, et

def flw2i4e(ex, ey, ep, D, eq=None):
    """
    Compute element stiffness (conductivity) matrix for 4 node isoparametric field element.

    Parameters
    ----------
    ex : array_like
        Element coordinates [x1, x2, x3, x4].
    ey : array_like
        Element coordinates [y1, y2, y3, y4].
    ep : array_like
        Element properties [t, ir], where t is thickness and ir is integration rule.
    D : array_like
        Constitutive matrix [[kxx, kxy], [kyx, kyy]].
    eq : float, optional
        Heat supply per unit volume.

    Returns
    -------
    Ke : ndarray
        Element 'stiffness' matrix, shape (4, 4).
    fe : ndarray, optional
        Element load vector, shape (4, 1), if eq is not None.
    """
    t = ep[0]
    ir = ep[1]
    ngp = ir*ir

    if eq is None:
        q = 0
    else:
        q = eq

    if ir == 1:
        g1 = 0.0
        w1 = 2.0
        gp = np.array([[g1, g1]])  # Make this explicitly 2D
        w = np.array([[w1, w1]])   # Make this explicitly 2D
    elif ir == 2:
        g1 = 0.577350269189626
        w1 = 1
        gp = np.array([
            [-g1, -g1],
            [g1, -g1],
            [-g1, g1],
            [g1, g1]
        ])
        w = np.array([
            [w1, w1],
            [w1, w1],
            [w1, w1],
            [w1, w1]
        ])
    elif ir == 3:
        g1 = 0.774596669241483
        g2 = 0.
        w1 = 0.555555555555555
        w2 = 0.888888888888888
        gp = np.array([
            [-g1, -g1],
            [-g2, -g1],
            [g1, -g1],
            [-g1, g2],
            [g2, g2],
            [g1, g2],
            [-g1, g1],
            [g2, g1],
            [g1, g1]
        ])
        w = np.array([
            [w1, w1],
            [w2, w1],
            [w1, w1],
            [w1, w2],
            [w2, w2],
            [w1, w2],
            [w1, w1],
            [w2, w1],
            [w1, w1]
        ])
    else:
        info("Used number of integration points not implemented")
        return
    
    # Make sure w is a 2D array
    w = np.asarray(w)
    if w.ndim == 1:
        w = w.reshape(-1, 1)
        
    # Compute the weight products - handle both matrix and array cases safely
    if w.shape[1] >= 2:
        wp = w[:, 0] * w[:, 1]
    else:
        # Handle the case where w might only have one column
        wp = w[:, 0]

    xsi = gp[:, 0]
    eta = gp[:, 1]
    r2 = ngp*2

    # Calculate shape functions and derivatives
    N = np.multiply((1-xsi), (1-eta))/4.
    N = np.column_stack((N, np.multiply((1+xsi), (1-eta))/4.))
    N = np.column_stack((N, np.multiply((1+xsi), (1+eta))/4.))
    N = np.column_stack((N, np.multiply((1-xsi), (1+eta))/4.))

    dNr = np.zeros((r2, 4))
    dNr[0:r2:2, 0] = -(1-eta)/4.
    dNr[0:r2:2, 1] = (1-eta)/4.
    dNr[0:r2:2, 2] = (1+eta)/4.
    dNr[0:r2:2, 3] = -(1+eta)/4.
    dNr[1:r2+1:2, 0] = -(1-xsi)/4.
    dNr[1:r2+1:2, 1] = -(1+xsi)/4.
    dNr[1:r2+1:2, 2] = (1+xsi)/4.
    dNr[1:r2+1:2, 3] = (1-xsi)/4.

    Ke1 = np.zeros((4, 4))
    fe1 = np.zeros((4, 1))
    
    # Convert ex and ey to arrays and combine for JT calculation
    ex_array = np.asarray(ex).reshape(-1, 1)
    ey_array = np.asarray(ey).reshape(-1, 1)
    coords = np.hstack([ex_array, ey_array])
    
    JT = dNr @ coords

    for i in range(ngp):
        indx = np.array([2*(i+1)-1, 2*(i+1)])
        detJ = np.linalg.det(JT[indx-1, :])
        if detJ < 10*np.finfo(float).eps:
            info("Jacobi determinant == 0")
        JTinv = np.linalg.inv(JT[indx-1, :])
        B = JTinv @ dNr[indx-1, :]

        Ke1 = Ke1 + B.T @ D @ B * detJ * wp[i]
        fe1 = fe1 + N[i, :].reshape(-1, 1) * detJ * wp[i]

    if eq is None:
        return Ke1 * t
    else:
        return Ke1 * t, fe1 * t * q

def flw2i4s(ex, ey, ep, D, ed):
    """
    Compute flows or corresponding quantities in the 4 node isoparametric element.

    Parameters
    ----------
    
    ex : array_like
        Element coordinates [x1, x2, x3, x4].
    ey : array_like
        Element coordinates [y1, y2, y3, y4].
    ep : array_like
        Element properties [t, ir], where t is thickness and ir is integration rule.
    D : array_like
        Constitutive matrix [[kxx, kxy], [kyx, kyy]].
    ed : array_like
        Element nodal values [u1, u2, u3, u4].

    Returns
    -------

    es : ndarray
        Element flows [[qx, qy], [.., ..]].
    et : ndarray
        Element gradients [[qx, qy], [..., ..]].
    eci : ndarray
        Gauss point location vector [[ix1, iy1], [..., ...], [ix(nint), iy(nint)]].
    """
    t = ep[0]
    ir = ep[1]
    ngp = ir*ir

    if ir == 1:
        g1 = 0.0
        w1 = 2.0
        gp = np.array([[g1, g1]])  # Make this explicitly 2D
        w = np.array([[w1, w1]])   # Make this explicitly 2D
    elif ir == 2:
        g1 = 0.577350269189626
        w1 = 1
        gp = np.array([
            [-g1, -g1],
            [g1, -g1],
            [-g1, g1],
            [g1, g1]
        ])
        w = np.array([
            [w1, w1],
            [w1, w1],
            [w1, w1],
            [w1, w1]
        ])
    elif ir == 3:
        g1 = 0.774596669241483
        g2 = 0.
        w1 = 0.555555555555555
        w2 = 0.888888888888888
        gp = np.array([
            [-g1, -g1],
            [-g2, -g1],
            [g1, -g1],
            [-g1, g2],
            [g2, g2],
            [g1, g2],
            [-g1, g1],
            [g2, g1],
            [g1, g1]
        ])
        w = np.array([
            [w1, w1],
            [w2, w1],
            [w1, w1],
            [w1, w2],
            [w2, w2],
            [w1, w2],
            [w1, w1],
            [w2, w1],
            [w1, w1]
        ])
    else:
        info("Used number of integration points not implemented")
    
    # Make sure w is a 2D array
    w = np.asarray(w)
    if w.ndim == 1:
        w = w.reshape(-1, 1)
        
    # We don't need wp for this function, so no need to compute it

    xsi = gp[:, 0]
    eta = gp[:, 1]
    r2 = ngp*2

    N = np.multiply((1-xsi), (1-eta))/4.
    N = np.column_stack((N, np.multiply((1+xsi), (1-eta))/4.))
    N = np.column_stack((N, np.multiply((1+xsi), (1+eta))/4.))
    N = np.column_stack((N, np.multiply((1-xsi), (1+eta))/4.))

    dNr = np.zeros((r2, 4))
    dNr[0:r2:2, 0] = -(1-eta)/4.
    dNr[0:r2:2, 1] = (1-eta)/4.
    dNr[0:r2:2, 2] = (1+eta)/4.
    dNr[0:r2:2, 3] = -(1+eta)/4.
    dNr[1:r2+1:2, 0] = -(1-xsi)/4.
    dNr[1:r2+1:2, 1] = -(1+xsi)/4.
    dNr[1:r2+1:2, 2] = (1+xsi)/4.
    dNr[1:r2+1:2, 3] = (1-xsi)/4.

    # Calculate Gauss point locations (for output)
    # Use numpy arrays instead of matrix multiplication
    ex_array = np.asarray(ex).reshape(-1, 1)
    ey_array = np.asarray(ey).reshape(-1, 1)
    coords = np.hstack([ex_array, ey_array])
    
    eci = N @ coords
    
    # Ensure ed is properly shaped for calculations
    if np.ndim(ed) == 1:
        ed = np.array([ed])

    # Get the number of rows in ed
    red = ed.shape[0]
    
    JT = dNr @ coords

    es = np.zeros((ngp*red, 2))
    et = np.zeros((ngp*red, 2))
    
    for i in range(ngp):
        indx = np.array([2*(i+1)-1, 2*(i+1)])
        detJ = np.linalg.det(JT[indx-1, :])
        if detJ < 10*np.finfo(float).eps:
            info("Jacobi determinant == 0")
        JTinv = np.linalg.inv(JT[indx-1, :])
        B = JTinv @ dNr[indx-1, :]

        # Process each row of ed
        for j in range(red):
            p1 = -D @ B @ ed[j].T
            p2 = B @ ed[j].T
            es[i + j*ngp, :] = p1.T
            et[i + j*ngp, :] = p2.T

    return es, et, eci

def flw2i8e(ex, ey, ep, D, eq=None):
    """
    Compute element stiffness (conductivity) matrix for 8 node isoparametric field element.

    Parameters
    ----------

    ex : array_like
        Element coordinates [x1, ..., x8].
    ey : array_like
        Element coordinates [y1, ..., y8].
    ep : array_like
        Element properties [t, ir], where t is thickness and ir is integration rule.
    D : array_like
        Constitutive matrix [[kxx, kxy], [kyx, kyy]].
    eq : float, optional
        Heat supply per unit volume.

    Returns
    -------

    Ke : ndarray
        Element 'stiffness' matrix, shape (8, 8).
    fe : ndarray, optional
        Element load vector, shape (8, 1), if eq is not None.
    """
    t = ep[0]
    ir = ep[1]
    ngp = ir*ir

    if eq is None:
        q = 0
    else:
        q = eq

    if ir == 1:
        g1 = 0.0
        w1 = 2.0
        gp = np.array([[g1, g1]])  # Make this explicitly 2D
        w = np.array([[w1, w1]])   # Make this explicitly 2D
    elif ir == 2:
        g1 = 0.577350269189626
        w1 = 1
        gp = np.array([
            [-g1, -g1],
            [g1, -g1],
            [-g1, g1],
            [g1, g1]
        ])
        w = np.array([
            [w1, w1],
            [w1, w1],
            [w1, w1],
            [w1, w1]
        ])
    elif ir == 3:
        g1 = 0.774596669241483
        g2 = 0.
        w1 = 0.555555555555555
        w2 = 0.888888888888888
        gp = np.array([
            [-g1, -g1],
            [-g2, -g1],
            [g1, -g1],
            [-g1, g2],
            [g2, g2],
            [g1, g2],
            [-g1, g1],
            [g2, g1],
            [g1, g1]
        ])
        w = np.array([
            [w1, w1],
            [w2, w1],
            [w1, w1],
            [w1, w2],
            [w2, w2],
            [w1, w2],
            [w1, w1],
            [w2, w1],
            [w1, w1]
        ])
    else:
        info("Used number of integration points not implemented")
        return
    
    # Make sure w is a 2D array
    w = np.asarray(w)
    if w.ndim == 1:
        w = w.reshape(-1, 1)
        
    # Compute the weight products - handle both matrix and array cases safely
    if w.shape[1] >= 2:
        wp = w[:, 0] * w[:, 1]
    else:
        # Handle the case where w might only have one column
        wp = w[:, 0]

    xsi = gp[:, 0]
    eta = gp[:, 1]
    r2 = ngp*2

    N = np.multiply(np.multiply(-(1-xsi), (1-eta)), (1+xsi+eta))/4.
    N = np.column_stack((
        N, 
        np.multiply(np.multiply(-(1+xsi), (1-eta)), (1-xsi+eta))/4.
    ))
    N = np.column_stack((
        N,
        np.multiply(np.multiply(-(1+xsi), (1+eta)), (1-xsi-eta))/4.
    ))
    N = np.column_stack((
        N,
        np.multiply(np.multiply(-(1-xsi), (1+eta)), (1+xsi-eta))/4.
    ))
    N = np.column_stack((
        N,
        np.multiply((1-np.multiply(xsi, xsi)), (1-eta))/2.
    ))
    N = np.column_stack((
        N,
        np.multiply((1+xsi), (1-np.multiply(eta, eta)))/2.
    ))
    N = np.column_stack((
        N,
        np.multiply((1-np.multiply(xsi, xsi)), (1+eta))/2.
    ))
    N = np.column_stack((
        N,
        np.multiply((1-xsi), (1-np.multiply(eta, eta)))/2.
    ))

    dNr = np.zeros((r2, 8))
    dNr[0:r2:2, 0] = -(-np.multiply((1-eta), (1+xsi+eta)) +
                       np.multiply((1-xsi), (1-eta)))/4.
    dNr[0:r2:2, 1] = -(np.multiply((1-eta), (1-xsi+eta)) -
                       np.multiply((1+xsi), (1-eta)))/4.
    dNr[0:r2:2, 2] = -(np.multiply((1+eta), (1-xsi-eta)) -
                       np.multiply((1+xsi), (1+eta)))/4.
    dNr[0:r2:2, 3] = -(-np.multiply((1+eta), (1+xsi-eta)) +
                       np.multiply((1-xsi), (1+eta)))/4.
    dNr[0:r2:2, 4] = -np.multiply(xsi, (1-eta))
    dNr[0:r2:2, 5] = (1-np.multiply(eta, eta))/2.
    dNr[0:r2:2, 6] = -np.multiply(xsi, (1+eta))
    dNr[0:r2:2, 7] = -(1-np.multiply(eta, eta))/2.
    dNr[1:r2+1:2, 0] = -(-np.multiply((1-xsi), (1+xsi+eta)) +
                         np.multiply((1-xsi), (1-eta)))/4.
    dNr[1:r2+1:2, 1] = -(-np.multiply((1+xsi), (1-xsi+eta)) +
                         np.multiply((1+xsi), (1-eta)))/4.
    dNr[1:r2+1:2, 2] = -(np.multiply((1+xsi), (1-xsi-eta)) -
                         np.multiply((1+xsi), (1+eta)))/4.
    dNr[1:r2+1:2, 3] = -(np.multiply((1-xsi), (1+xsi-eta)) -
                         np.multiply((1-xsi), (1+eta)))/4.
    dNr[1:r2+1:2, 4] = -(1-np.multiply(xsi, xsi))/2.
    dNr[1:r2+1:2, 5] = -np.multiply(eta, (1+xsi))
    dNr[1:r2+1:2, 6] = (1-np.multiply(xsi, xsi))/2.
    dNr[1:r2+1:2, 7] = -np.multiply(eta, (1-xsi))

    Ke1 = np.zeros((8, 8))
    fe1 = np.zeros((8, 1))
    
    # Convert ex and ey to arrays and combine for JT calculation
    ex_array = np.asarray(ex).reshape(-1, 1)
    ey_array = np.asarray(ey).reshape(-1, 1)
    coords = np.hstack([ex_array, ey_array])
    
    JT = dNr @ coords

    for i in range(ngp):
        indx = np.array([2*(i+1)-1, 2*(i+1)])
        detJ = np.linalg.det(JT[indx-1, :])
        if detJ < 10*np.finfo(float).eps:
            info("Jacobideterminanten lika med noll!")
        JTinv = np.linalg.inv(JT[indx-1, :])
        B = JTinv @ dNr[indx-1, :]
        Ke1 = Ke1 + B.T @ D @ B * detJ * wp[i]
        fe1 = fe1 + N[i, :].reshape(-1, 1) * detJ * wp[i]

    if eq != None:
        return Ke1*t, fe1*t*q
    else:
        return Ke1*t

def flw2i8s(ex, ey, ep, D, ed):
    """
    Compute flows or corresponding quantities in the 8 node isoparametric element.
    
    Parameters
    ----------

    ex : array_like
        Element coordinates [x1, x2, x3, ..., x8].
    ey : array_like
        Element coordinates [y1, y2, y3, ..., y8].
    ep : array_like
        Element properties [t, ir], where t is thickness and ir is integration rule.
    D : array_like
        Constitutive matrix [[kxx, kxy], [kyx, kyy]].
    ed : array_like
        Element nodal values [u1, ..., u8].

    Returns
    -------

    es : ndarray
        Element flows [[qx, qy], [.., ..]].
    et : ndarray
        Element gradients [[qx, qy], [.., ..]].
    eci : ndarray
        Gauss point location vector [[ix1, iy1], [..., ...], [ix(nint), iy(nint)]].
    """
    t = ep[0]
    ir = ep[1]
    ngp = ir*ir

    if ir == 1:
        g1 = 0.0
        w1 = 2.0
        gp = np.array([[g1, g1]])  # Make this explicitly 2D
        w = np.array([[w1, w1]])   # Make this explicitly 2D
    elif ir == 2:
        g1 = 0.577350269189626
        w1 = 1
        gp = np.array([
            [-g1, -g1],
            [g1, -g1],
            [-g1, g1],
            [g1, g1]
        ])
        w = np.array([
            [w1, w1],
            [w1, w1],
            [w1, w1],
            [w1, w1]
        ])
    elif ir == 3:
        g1 = 0.774596669241483
        g2 = 0.
        w1 = 0.555555555555555
        w2 = 0.888888888888888
        gp = np.array([
            [-g1, -g1],
            [-g2, -g1],
            [g1, -g1],
            [-g1, g2],
            [g2, g2],
            [g1, g2],
            [-g1, g1],
            [g2, g1],
            [g1, g1]
        ])
        w = np.array([
            [w1, w1],
            [w2, w1],
            [w1, w1],
            [w1, w2],
            [w2, w2],
            [w1, w2],
            [w1, w1],
            [w2, w1],
            [w1, w1]
        ])
    else:
        info("Used number of integration points not implemented")
        return
    
    # We don't need wp for this function, so no need to compute it

    xsi = gp[:, 0]
    eta = gp[:, 1]
    r2 = ngp*2

    N = np.multiply(np.multiply(-(1-xsi), (1-eta)), (1+xsi+eta))/4.
    N = np.column_stack((
        N, 
        np.multiply(np.multiply(-(1+xsi), (1-eta)), (1-xsi+eta))/4.
    ))
    N = np.column_stack((
        N,
        np.multiply(np.multiply(-(1+xsi), (1+eta)), (1-xsi-eta))/4.
    ))
    N = np.column_stack((
        N,
        np.multiply(np.multiply(-(1-xsi), (1+eta)), (1+xsi-eta))/4.
    ))
    N = np.column_stack((
        N,
        np.multiply((1-np.multiply(xsi, xsi)), (1-eta))/2.
    ))
    N = np.column_stack((
        N,
        np.multiply((1+xsi), (1-np.multiply(eta, eta)))/2.
    ))
    N = np.column_stack((
        N,
        np.multiply((1-np.multiply(xsi, xsi)), (1+eta))/2.
    ))
    N = np.column_stack((
        N,
        np.multiply((1-xsi), (1-np.multiply(eta, eta)))/2.
    ))

    dNr = np.zeros((r2, 8))
    dNr[0:r2:2, 0] = -(-np.multiply((1-eta), (1+xsi+eta)) +
                       np.multiply((1-xsi), (1-eta)))/4.
    dNr[0:r2:2, 1] = -(np.multiply((1-eta), (1-xsi+eta)) -
                       np.multiply((1+xsi), (1-eta)))/4.
    dNr[0:r2:2, 2] = -(np.multiply((1+eta), (1-xsi-eta)) -
                       np.multiply((1+xsi), (1+eta)))/4.
    dNr[0:r2:2, 3] = -(-np.multiply((1+eta), (1+xsi-eta)) +
                       np.multiply((1-xsi), (1+eta)))/4.
    dNr[0:r2:2, 4] = -np.multiply(xsi, (1-eta))
    dNr[0:r2:2, 5] = (1-np.multiply(eta, eta))/2.
    dNr[0:r2:2, 6] = -np.multiply(xsi, (1+eta))
    dNr[0:r2:2, 7] = -(1-np.multiply(eta, eta))/2.
    dNr[1:r2+1:2, 0] = -(-np.multiply((1-xsi), (1+xsi+eta)) +
                         np.multiply((1-xsi), (1-eta)))/4.
    dNr[1:r2+1:2, 1] = -(-np.multiply((1+xsi), (1-xsi+eta)) +
                         np.multiply((1+xsi), (1-eta)))/4.
    dNr[1:r2+1:2, 2] = -(np.multiply((1+xsi), (1-xsi-eta)) -
                         np.multiply((1+xsi), (1+eta)))/4.
    dNr[1:r2+1:2, 3] = -(np.multiply((1-xsi), (1+xsi-eta)) -
                         np.multiply((1-xsi), (1+eta)))/4.
    dNr[1:r2+1:2, 4] = -(1-np.multiply(xsi, xsi))/2.
    dNr[1:r2+1:2, 5] = -np.multiply(eta, (1+xsi))
    dNr[1:r2+1:2, 6] = (1-np.multiply(xsi, xsi))/2.
    dNr[1:r2+1:2, 7] = -np.multiply(eta, (1-xsi))

    # Calculate Gauss point locations (for output)
    # Use numpy arrays instead of matrix multiplication
    ex_array = np.asarray(ex).reshape(-1, 1)
    ey_array = np.asarray(ey).reshape(-1, 1)
    coords = np.hstack([ex_array, ey_array])
    
    eci = N @ coords
    
    # Ensure ed is properly shaped for calculations
    if np.ndim(ed) == 1:
        ed = np.array([ed])

    # Get the number of rows in ed
    red = ed.shape[0]
    
    JT = dNr @ coords

    es = np.zeros((ngp*red, 2))
    et = np.zeros((ngp*red, 2))
    
    for i in range(ngp):
        indx = np.array([2*(i+1)-1, 2*(i+1)])
        detJ = np.linalg.det(JT[indx-1, :])
        if detJ < 10*np.finfo(float).eps:
            info("Jacobi determinant == 0")
        JTinv = np.linalg.inv(JT[indx-1, :])
        B = JTinv @ dNr[indx-1, :]

        # Process each row of ed
        for j in range(red):
            p1 = -D @ B @ ed[j].T
            p2 = B @ ed[j].T
            es[i + j*ngp, :] = p1.T
            et[i + j*ngp, :] = p2.T

    return es, et, eci

def flw3i8e(ex, ey, ez, ep, D, eq=None):
    """
    Compute element stiffness (conductivity) matrix for 8 node isoparametric field element.
    
    Parameters
    ----------

    ex : array_like
        Element node x-coordinates [x1, x2, x3, ..., x8].
    ey : array_like
        Element node y-coordinates [y1, y2, y3, ..., y8].
    ez : array_like
        Element node z-coordinates [z1, z2, z3, ..., z8].
    ep : array_like
        Element properties [ir], where ir is integration rule.
    D : array_like
        Constitutive matrix [[kxx, kxy, kxz], [kyx, kyy, kyz], [kzx, kzy, kzz]].
    eq : float, optional
        Heat supply per unit volume.

    Returns
    -------

    Ke : ndarray
        Element 'stiffness' matrix, shape (8, 8).
    fe : ndarray, optional
        Element load vector, shape (8, 1), if eq is not None.
    """
    ir = ep[0]
    ngp = ir*ir*ir

    if eq is None:
        q = 0
    else:
        q = eq

    if ir == 2:
        g1 = 0.577350269189626
        w1 = 1
        gp = np.array([
            [-1, -1, -1],
            [1, -1, -1],
            [1, 1, -1],
            [-1, 1, -1],
            [-1, -1, 1],
            [1, -1, 1],
            [1, 1, 1],
            [-1, 1, 1]
        ])*g1
        w = np.ones((8, 3))*w1
    elif ir == 3:
        g1 = 0.774596669241483
        g2 = 0.
        w1 = 0.555555555555555
        w2 = 0.888888888888888
        gp = np.zeros((27, 3))
        w = np.zeros((27, 3))
        I1 = np.array([-1, 0, 1, -1, 0, 1, -1, 0, 1])
        I2 = np.array([0, -1, 0, 0, 1, 0, 0, 1, 0])
        gp[:, 0] = np.repeat(I1, 3).reshape(27, 1).flatten()
        I2_tiled = np.tile(I2, 3).reshape(27, 1)
        gp[:, 0] = gp[:, 0]*g1 + I2_tiled.flatten()*g2
        I1 = np.abs(I1)
        I2 = np.abs(I2)
        w[:, 0] = np.repeat(I1, 3).reshape(27, 1).flatten()*w1
        w[:, 0] = w[:, 0] + np.tile(I2, 3).reshape(27, 1).flatten()*w2
        I1 = np.repeat([-1, -1, -1, 0, 0, 0, 1, 1, 1], 3)
        I2 = np.repeat([0, 0, 0, 1, 1, 1, 0, 0, 0], 3)
        gp[:, 1] = I1*g1
        gp[:, 1] = gp[:, 1] + I2*g2
        I1 = np.abs(I1)
        I2 = np.abs(I2)
        w[:, 1] = I1*w1
        w[:, 1] = w[:, 1] + I2*w2
        I1 = np.tile([-1, -1, -1, -1, -1, -1, -1, -1, -1], 3)
        I2 = np.tile([0, 0, 0, 0, 0, 0, 0, 0, 0], 3)
        I3 = np.abs(I1)
        gp[:, 2] = I1*g1
        gp[:, 2] = gp[:, 2] + I2*g2
        w[:, 2] = I3*w1
        w[:, 2] = w[:, 2] + I2*w2
    else:
        info("Used number of integration points not implemented")
        return

    wp = w[:, 0]*w[:, 1]*w[:, 2]

    xsi = gp[:, 0]
    eta = gp[:, 1]
    zet = gp[:, 2]
    r2 = ngp*3

    N = np.zeros((ngp, 8))
    dNr = np.zeros((r2, 8))

    N[:, 0] = (1-xsi)*(1-eta)*(1-zet)/8
    N[:, 1] = (1+xsi)*(1-eta)*(1-zet)/8
    N[:, 2] = (1+xsi)*(1+eta)*(1-zet)/8
    N[:, 3] = (1-xsi)*(1+eta)*(1-zet)/8
    N[:, 4] = (1-xsi)*(1-eta)*(1+zet)/8
    N[:, 5] = (1+xsi)*(1-eta)*(1+zet)/8
    N[:, 6] = (1+xsi)*(1+eta)*(1+zet)/8
    N[:, 7] = (1-xsi)*(1+eta)*(1+zet)/8

    dNr[0:r2+1:3, 0] = -(1-eta)*(1-zet)
    dNr[0:r2+1:3, 1] = (1-eta)*(1-zet)
    dNr[0:r2+1:3, 2] = (1+eta)*(1-zet)
    dNr[0:r2+1:3, 3] = -(1+eta)*(1-zet)
    dNr[0:r2+1:3, 4] = -(1-eta)*(1+zet)
    dNr[0:r2+1:3, 5] = (1-eta)*(1+zet)
    dNr[0:r2+1:3, 6] = (1+eta)*(1+zet)
    dNr[0:r2+1:3, 7] = -(1+eta)*(1+zet)
    dNr[1:r2+2:3, 0] = -(1-xsi)*(1-zet)
    dNr[1:r2+2:3, 1] = -(1+xsi)*(1-zet)
    dNr[1:r2+2:3, 2] = (1+xsi)*(1-zet)
    dNr[1:r2+2:3, 3] = (1-xsi)*(1-zet)
    dNr[1:r2+2:3, 4] = -(1-xsi)*(1+zet)
    dNr[1:r2+2:3, 5] = -(1+xsi)*(1+zet)
    dNr[1:r2+2:3, 6] = (1+xsi)*(1+zet)
    dNr[1:r2+2:3, 7] = (1-xsi)*(1+zet)
    dNr[2:r2+3:3, 0] = -(1-xsi)*(1-eta)
    dNr[2:r2+3:3, 1] = -(1+xsi)*(1-eta)
    dNr[2:r2+3:3, 2] = -(1+xsi)*(1+eta)
    dNr[2:r2+3:3, 3] = -(1-xsi)*(1+eta)
    dNr[2:r2+3:3, 4] = (1-xsi)*(1-eta)
    dNr[2:r2+3:3, 5] = (1+xsi)*(1-eta)
    dNr[2:r2+3:3, 6] = (1+xsi)*(1+eta)
    dNr[2:r2+3:3, 7] = (1-xsi)*(1+eta)

    dNr = dNr/8.0

    Ke = np.zeros((8, 8))
    fe = np.zeros((8, 1))

    ex = np.asarray(ex).reshape((8, 1))
    ey = np.asarray(ey).reshape((8, 1))
    ez = np.asarray(ez).reshape((8, 1))
    coords = np.hstack((ex, ey, ez))

    JT = dNr @ coords

    eps = np.finfo(float).eps

    for i in range(ngp):
        indx = [i*3, i*3+1, i*3+2]
        detJ = np.linalg.det(JT[indx, :])
        if detJ < 10*eps:
            info("Jacobi determinant equal or less than zero!")
        JTinv = np.linalg.inv(JT[indx, :])
        dNx = JTinv @ dNr[indx, :]

        B = np.zeros((3, 8))
        B[0, :] = dNx[0, :]
        B[1, :] = dNx[1, :]
        B[2, :] = dNx[2, :]

        Ke = Ke + B.T @ D @ B * detJ * wp[i]
        fe = fe + (N[i, :].reshape(-1, 1) * detJ * wp[i]).reshape(-1, 1)

    if eq is not None:
        return Ke, fe*q
    else:
        return Ke

def flw3i8s(ex, ey, ez, ep, D, ed):
    """
    Compute flows or corresponding quantities in the 8 node (3-dim) isoparametric field element.
    
    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2, x3, ..., x8].
    ey : array_like
        Element node y-coordinates [y1, y2, y3, ..., y8].
    ez : array_like
        Element node z-coordinates [z1, z2, z3, ..., z8].
    ep : array_like
        Element properties [ir], where ir is integration rule.
    D : array_like
        Constitutive matrix [[kxx, kxy, kxz], [kyx, kyy, kyz], [kzx, kzy, kzz]].
    ed : array_like
        Element nodal values [[u1, ..., u8], [.., ..., ..]].

    Returns
    -------
    es : ndarray
        Element flows [[qx, qy, qz], [.., .., ..]].
    et : ndarray
        Element gradients [[qx, qy, qz], [.., .., ..]].
    eci : ndarray
        Gauss point location vector [[ix1, iy1, iz1], [..., ..., ...], [ix(nint), iy(nint), iz(nint)]].
    """
    ir = ep[0]
    ngp = ir*ir*ir

    if ir == 2:
        g1 = 0.577350269189626
        w1 = 1
        gp = np.array([
            [-1, -1, -1],
            [1, -1, -1],
            [1, 1, -1],
            [-1, 1, -1],
            [-1, -1, 1],
            [1, -1, 1],
            [1, 1, 1],
            [-1, 1, 1]
        ])*g1
        w = np.ones((8, 3))*w1
    elif ir == 3:
        g1 = 0.774596669241483
        g2 = 0.
        w1 = 0.555555555555555
        w2 = 0.888888888888888
        gp = np.zeros((27, 3))
        w = np.zeros((27, 3))
        I1 = np.array([-1, 0, 1, -1, 0, 1, -1, 0, 1])
        I2 = np.array([0, -1, 0, 0, 1, 0, 0, 1, 0])
        gp[:, 0] = np.repeat(I1, 3).reshape(27, 1).flatten()
        I2_tiled = np.tile(I2, 3).reshape(27, 1)
        gp[:, 0] = gp[:, 0]*g1 + I2_tiled.flatten()*g2
        I1 = np.abs(I1)
        I2 = np.abs(I2)
        w[:, 0] = np.repeat(I1, 3).reshape(27, 1).flatten()*w1
        w[:, 0] = w[:, 0] + np.tile(I2, 3).reshape(27, 1).flatten()*w2
        I1 = np.repeat([-1, -1, -1, 0, 0, 0, 1, 1, 1], 3)
        I2 = np.repeat([0, 0, 0, 1, 1, 1, 0, 0, 0], 3)
        gp[:, 1] = I1*g1
        gp[:, 1] = gp[:, 1] + I2*g2
        I1 = np.abs(I1)
        I2 = np.abs(I2)
        w[:, 1] = I1*w1
        w[:, 1] = w[:, 1] + I2*w2
        I1 = np.tile([-1, -1, -1, -1, -1, -1, -1, -1, -1], 3)
        I2 = np.tile([0, 0, 0, 0, 0, 0, 0, 0, 0], 3)
        I3 = np.abs(I1)
        gp[:, 2] = I1*g1
        gp[:, 2] = gp[:, 2] + I2*g2
        w[:, 2] = I3*w1
        w[:, 2] = w[:, 2] + I2*w2
    else:
        info("Used number of integration points not implemented")
        return

    wp = w[:, 0]*w[:, 1]*w[:, 2]

    xsi = gp[:, 0]
    eta = gp[:, 1]
    zet = gp[:, 2]
    r2 = ngp*3

    N = np.zeros((ngp, 8))
    dNr = np.zeros((r2, 8))

    N[:, 0] = (1-xsi)*(1-eta)*(1-zet)/8
    N[:, 1] = (1+xsi)*(1-eta)*(1-zet)/8
    N[:, 2] = (1+xsi)*(1+eta)*(1-zet)/8
    N[:, 3] = (1-xsi)*(1+eta)*(1-zet)/8
    N[:, 4] = (1-xsi)*(1-eta)*(1+zet)/8
    N[:, 5] = (1+xsi)*(1-eta)*(1+zet)/8
    N[:, 6] = (1+xsi)*(1+eta)*(1+zet)/8
    N[:, 7] = (1-xsi)*(1+eta)*(1+zet)/8

    dNr[0:r2+1:3, 0] = -(1-eta)*(1-zet)
    dNr[0:r2+1:3, 1] = (1-eta)*(1-zet)
    dNr[0:r2+1:3, 2] = (1+eta)*(1-zet)
    dNr[0:r2+1:3, 3] = -(1+eta)*(1-zet)
    dNr[0:r2+1:3, 4] = -(1-eta)*(1+zet)
    dNr[0:r2+1:3, 5] = (1-eta)*(1+zet)
    dNr[0:r2+1:3, 6] = (1+eta)*(1+zet)
    dNr[0:r2+1:3, 7] = -(1+eta)*(1+zet)
    dNr[1:r2+2:3, 0] = -(1-xsi)*(1-zet)
    dNr[1:r2+2:3, 1] = -(1+xsi)*(1-zet)
    dNr[1:r2+2:3, 2] = (1+xsi)*(1-zet)
    dNr[1:r2+2:3, 3] = (1-xsi)*(1-zet)
    dNr[1:r2+2:3, 4] = -(1-xsi)*(1+zet)
    dNr[1:r2+2:3, 5] = -(1+xsi)*(1+zet)
    dNr[1:r2+2:3, 6] = (1+xsi)*(1+zet)
    dNr[1:r2+2:3, 7] = (1-xsi)*(1+zet)
    dNr[2:r2+3:3, 0] = -(1-xsi)*(1-eta)
    dNr[2:r2+3:3, 1] = -(1+xsi)*(1-eta)
    dNr[2:r2+3:3, 2] = -(1+xsi)*(1+eta)
    dNr[2:r2+3:3, 3] = -(1-xsi)*(1+eta)
    dNr[2:r2+3:3, 4] = (1-xsi)*(1-eta)
    dNr[2:r2+3:3, 5] = (1+xsi)*(1-eta)
    dNr[2:r2+3:3, 6] = (1+xsi)*(1+eta)
    dNr[2:r2+3:3, 7] = (1-xsi)*(1+eta)

    dNr = dNr/8.0

    # Calculate Gauss point locations (for output)
    ex = np.asarray(ex).reshape((8, 1))
    ey = np.asarray(ey).reshape((8, 1))
    ez = np.asarray(ez).reshape((8, 1))
    coords = np.hstack((ex, ey, ez))
    
    eci = N @ coords
    
    # Ensure ed is properly shaped for calculations
    if np.ndim(ed) == 1:
        ed = np.array([ed])

    # Get the number of rows in ed
    red = ed.shape[0]
    
    JT = dNr @ coords

    es = np.zeros((ngp*red, 3))
    et = np.zeros((ngp*red, 3))
    
    eps = np.finfo(float).eps
    
    for i in range(ngp):
        indx = [i*3, i*3+1, i*3+2]
        detJ = np.linalg.det(JT[indx, :])
        if detJ < 10*eps:
            info("Jacobi determinant equal or less than zero!")
        JTinv = np.linalg.inv(JT[indx, :])
        dNx = JTinv @ dNr[indx, :]

        B = np.zeros((3, 8))
        B[0, :] = dNx[0, :]
        B[1, :] = dNx[1, :]
        B[2, :] = dNx[2, :]

        # Process each row of ed
        for j in range(red):
            p1 = -D @ B @ ed[j].T
            p2 = B @ ed[j].T
            es[i + j*ngp, :] = p1.T
            et[i + j*ngp, :] = p2.T

    return es, et, eci

def plante(ex, ey, ep, D, eq=None):
    """
    Calculate the stiffness matrix for a triangular plane stress or plane strain element.
    
    Parameters
    ----------
    ex : array_like
        Element coordinates [x1, x2, x3].
    ey : array_like
        Element coordinates [y1, y2, y3].
    ep : array_like
        Element properties [ptype, t], where ptype is analysis type and t is thickness.
    D : array_like
        Constitutive matrix.
    eq : array_like, optional
        Body force vector [bx, by], where bx, by are body forces in x, y directions.
              
    Returns
    -------
    Ke : ndarray
        Element stiffness matrix, shape (6, 6).
    fe : ndarray, optional
        Equivalent nodal forces, shape (6, 1), if eq is given.
    """

    ptype, t = ep

    bx = 0.0
    by = 0.0

    if not eq is None:
        bx = eq[0]
        by = eq[1]

    # Ensure arrays have correct dimensions for operations
    ex_array = np.asarray(ex).reshape(-1, 1)
    ey_array = np.asarray(ey).reshape(-1, 1)
    ones_col = np.ones((3, 1))
    
    # Construct C matrix with proper dimensions
    C = np.array([
        [1, ex[0], ey[0], 0,     0,     0],
        [0,     0,     0, 1, ex[0], ey[0]],
        [1, ex[1], ey[1], 0,     0,     0],
        [0,     0,     0, 1, ex[1], ey[1]],
        [1, ex[2], ey[2], 0,     0,     0],
        [0,     0,     0, 1, ex[2], ey[2]]
    ])

    # Calculate area using a properly shaped array
    A_matrix = np.hstack([ones_col, ex_array, ey_array])
    A = 0.5 * np.linalg.det(A_matrix)

    # --------- plane stress --------------------------------------

    if ptype == 1:
        B = np.array([
            [0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 1],
            [0, 0, 1, 0, 1, 0]
        ]) @ np.linalg.inv(C)

        colD = D.shape[1]

        if colD > 3:
            Cm = np.linalg.inv(D)
            Dm = np.linalg.inv(Cm[np.ix_((0, 1, 3), (0, 1, 3))])
        else:
            Dm = D

        Ke = B.T @ Dm @ B * A * t
        fe = A/3 * np.array([bx, by, bx, by, bx, by]).reshape(-1, 1) * t

        if eq is None:
            return Ke
        else:
            return Ke, fe.T

    #--------- plane strain --------------------------------------

    elif ptype == 2:
        B = np.array([
            [0, 1, 0, 0, 0, 0, ],
            [0, 0, 0, 0, 0, 1, ],
            [0, 0, 1, 0, 1, 0, ]
        ]) @ np.linalg.inv(C)

        colD = D.shape[1]

        if colD > 3:
            Dm = D[np.ix_((0, 1, 3), (0, 1, 3))]
        else:
            Dm = D

        Ke = B.T @ Dm @ B * A * t
        fe = A/3 * np.array([bx, by, bx, by, bx, by]).reshape(-1, 1) * t

        if eq is None:
            return Ke
        else:
            return Ke, fe.T

    else:
        info("Error ! Check first argument, ptype=1 or 2 allowed")
        if eq is None:
            return None
        else:
            return None, None

def plants(ex, ey, ep, D, ed):
    """
    Calculate element normal and shear stress for a triangular plane stress or plane strain element.

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2, x3].
    ey : array_like
        Element node y-coordinates [y1, y2, y3].
    ep : array_like
        Element properties [ptype, t], where ptype is analysis type and t is thickness.
    D : array_like
        Constitutive matrix.
    ed : array_like
        Element displacement vector [u1, u2, ..., u6], one row for each element.

    Returns
    -------
    es : ndarray
        Element stress matrix, one row for each element.
        Each row contains [sigx, sigy, [sigz], tauxy].
    et : ndarray
        Element strain matrix, one row for each element.
        Each row contains [epsx, epsy, [epsz], gamxy].
    """

    ptype = ep[0]

    # Ensure all inputs are proper arrays
    ex = np.asarray(ex)
    ey = np.asarray(ey)
    ed = np.asarray(ed)

    if np.ndim(ex) == 1:
        ex = np.array([ex])
    if np.ndim(ey) == 1:
        ey = np.array([ey])
    if np.ndim(ed) == 1:
        ed = np.array([ed])

    rowed = ed.shape[0]
    rowex = ex.shape[0]

    # --------- plane stress --------------------------------------

    if ptype == 1:

        colD = D.shape[1]

        if colD > 3:
            Cm = np.linalg.inv(D)
            Dm = np.linalg.inv(Cm[np.ix_((0, 1, 3), (0, 1, 3))])
        else:
            Dm = D

        incie = 0

        if rowex == 1:
            incie = 0
        else:
            incie = 1

        et = np.zeros([rowed, colD])
        es = np.zeros([rowed, colD])

        ie = 0

        for i in range(rowed):
            # Create C matrix for element
            ex1 = ex[ie, :].flatten()
            ey1 = ey[ie, :].flatten()
            
            # Construct the C matrix properly
            C = np.array([
                [1, ex1[0], ey1[0], 0, 0, 0],
                [0, 0, 0, 1, ex1[0], ey1[0]],
                [1, ex1[1], ey1[1], 0, 0, 0],
                [0, 0, 0, 1, ex1[1], ey1[1]],
                [1, ex1[2], ey1[2], 0, 0, 0],
                [0, 0, 0, 1, ex1[2], ey1[2]]
            ])
            
            # Create the element B matrix
            B = np.array([
                [0, 1, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 1],
                [0, 0, 1, 0, 1, 0]]) @ np.linalg.inv(C)

            # Make sure element displacements are the right shape
            element_disp = ed[ie, :].flatten()
            
            # Calculate strains
            ee = B @ element_disp.reshape(-1, 1)

            if colD > 3:
                ss = np.zeros([colD, 1])
                ss[[0, 1, 3]] = Dm @ ee
                ee = Cm @ ss
            else:
                ss = Dm @ ee

            et[ie, :] = ee.T.flatten()
            es[ie, :] = ss.T.flatten()

            ie = ie + incie

        return es, et

    # --------- plane strain --------------------------------------
    elif ptype == 2:  # Implementation by LAPM
        colD = D.shape[1]
        incie = 0

        if rowex == 1:
            incie = 0
        else:
            incie = 1

        et = np.zeros([rowed, colD])
        es = np.zeros([rowed, colD])

        ie = 0

        ee = np.zeros([colD, 1])

        for i in range(rowed):
            # Create C matrix for element
            ex1 = ex[ie, :].flatten()
            ey1 = ey[ie, :].flatten()
            
            # Construct the C matrix properly
            C = np.array([
                [1, ex1[0], ey1[0], 0, 0, 0],
                [0, 0, 0, 1, ex1[0], ey1[0]],
                [1, ex1[1], ey1[1], 0, 0, 0],
                [0, 0, 0, 1, ex1[1], ey1[1]],
                [1, ex1[2], ey1[2], 0, 0, 0],
                [0, 0, 0, 1, ex1[2], ey1[2]]
            ])

            # Create the element B matrix
            B = np.array([
                [0, 1, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 1],
                [0, 0, 1, 0, 1, 0]]) @ np.linalg.inv(C)

            # Make sure element displacements are the right shape
            element_disp = ed[ie, :].flatten()
            
            # Calculate strains
            e = B @ element_disp.reshape(-1, 1)

            if colD > 3:
                ee[[0, 1, 3]] = e
            else:
                ee = e

            et[ie, :] = ee.T.flatten()
            es[ie, :] = (D @ ee).T.flatten()

            ie = ie + incie

        return es, et

    else:
        print("Error ! Check first argument, ptype=1 or 2 allowed")
        return None
    
def plantf(ex, ey, ep, es):
    """
    Compute internal element force vector in a triangular element in plane stress or plane strain.

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2, x3].
    ey : array_like
        Element node y-coordinates [y1, y2, y3].
    ep : array_like
        Element properties [ptype, t], where ptype is analysis type and t is thickness.
    es : array_like
        Element stress matrix [[sigx, sigy, [sigz], tauxy], [...]], one row for each element.

    Returns
    -------
    fe : ndarray
        Internal force vector [[f1], [f2], ..., [f8]].
    """

    ptype, t = ep

    colD = es.shape[1]

    #--------- plane stress --------------------------------------

    if ptype == 1:

        C = np.matrix([
            [1, ex[0], ey[0], 0, 0,     0],
            [0, 0,     0,     1, ex[0], ey[0]],
            [1, ex[1], ey[1], 0, 0,     0],
            [0, 0,     0,     1, ex[1], ey[1]],
            [1, ex[2], ey[2], 0, 0,     0],
            [0, 0,     0,     1, ex[2], ey[2]]
        ])

        A = 0.5*np.linalg.det(np.matrix([
            [1, ex[0], ey[0]],
            [1, ex[1], ey[1]],
            [1, ex[2], ey[2]]
        ]))

        B = np.matrix([
            [0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 1],
            [0, 0, 1, 0, 1, 0]
        ])*np.linalg.inv(C)

        if colD > 3:
            stress = np.asmatrix(es[np.ix_((0, 1, 3))])
        else:
            stress = np.asmatrix(es)

        ef = (A*t*B.T*stress.T).T

        return np.reshape(np.asarray(ef), 6)

    #--------- plane strain --------------------------------------

    elif ptype == 2:

        C = np.matrix([
            [1, ex[0], ey[0], 0, 0,     0],
            [0, 0,     0,     1, ex[0], ey[0]],
            [1, ex[1], ey[1], 0, 0,     0],
            [0, 0,     0,     1, ex[1], ey[1]],
            [1, ex[2], ey[2], 0, 0,     0],
            [0, 0,     0,     1, ex[2], ey[2]]
        ])

        A = 0.5*np.linalg.det(np.matrix([
            [1, ex[0], ey[0]],
            [1, ex[1], ey[1]],
            [1, ex[2], ey[2]]
        ]))

        B = np.matrix([
            [0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 1],
            [0, 0, 1, 0, 1, 0]
        ])*np.linalg.inv(C)

        if colD > 3:
            stress = np.asmatrix(es[np.ix_((0, 1, 3))])
        else:
            stress = np.asmatrix(es)

        ef = (A*t*B.T*stress.T).T

        return np.reshape(np.asarray(ef), (6,1))

    else:
        info("Error ! Check first argument, ptype=1 or 2 allowed")
        return None


def platre(ex, ey, ep, D, eq=None):
    """
    Calculate the stiffness matrix for a rectangular plate element.

    NOTE! Element sides must be parallel to the coordinate axis.

    Parameters
    ----------
    ex : array_like
        Element coordinates [x1, x2, x3, x4].
    ey : array_like
        Element coordinates [y1, y2, y3, y4].
    ep : array_like
        Element properties [t], where t is thickness.
    D : array_like
        Constitutive matrix for plane stress.
    eq : array_like, optional
        Load per unit area [qz].

    Returns
    -------
    Ke : ndarray
        Element stiffness matrix, shape (12, 12).
    fe : ndarray, optional
        Equivalent nodal forces, shape (12, 1), if eq is not None.
    """
    Lx = (ex[2]-ex[0]).astype(float)
    Ly = (ey[2]-ey[0]).astype(float)
    t = ep[0]

    D = t**3/12.*D

    A1 = Ly/(Lx**3)
    A2 = Lx/(Ly**3)
    A3 = 1/Lx/Ly
    A4 = Ly/(Lx**2)
    A5 = Lx/(Ly**2)
    A6 = 1/Lx
    A7 = 1/Ly
    A8 = Ly/Lx
    A9 = Lx/Ly

    C1 = 4*A1*D[0, 0]+4*A2*D[1, 1]+2*A3*D[0, 1]+5.6*A3*D[2, 2]
    C2 = -4*A1*D[0, 0]+2*A2*D[1, 1]-2*A3*D[0, 1]-5.6*A3*D[2, 2]
    C3 = 2*A1*D[0, 0]-4*A2*D[1, 1]-2*A3*D[0, 1]-5.6*A3*D[2, 2]
    C4 = -2*A1*D[0, 0]-2*A2*D[1, 1]+2*A3*D[0, 1]+5.6*A3*D[2, 2]
    C5 = 2*A5*D[1, 1]+A6*D[0, 1]+0.4*A6*D[2, 2]
    C6 = 2*A4*D[0, 0]+A7*D[0, 1]+0.4*A7*D[2, 2]

    C7 = 2*A5*D[1, 1]+0.4*A6*D[2, 2]
    C8 = 2*A4*D[0, 0]+0.4*A7*D[2, 2]
    C9 = A5*D[1, 1]-A6*D[0, 1]-0.4*A6*D[2, 2]
    C10 = A4*D[0, 0]-A7*D[0, 1]-0.4*A7*D[2, 2]
    C11 = A5*D[1, 1]-0.4*A6*D[2, 2]
    C12 = A4*D[0, 0]-0.4*A7*D[2, 2]

    C13 = 4/3.*A9*D[1, 1]+8/15.*A8*D[2, 2]
    C14 = 4/3.*A8*D[0, 0]+8/15.*A9*D[2, 2]
    C15 = 2/3.*A9*D[1, 1]-8/15.*A8*D[2, 2]
    C16 = 2/3.*A8*D[0, 0]-8/15.*A9*D[2, 2]
    C17 = 2/3.*A9*D[1, 1]-2/15.*A8*D[2, 2]
    C18 = 2/3.*A8*D[0, 0]-2/15.*A9*D[2, 2]
    C19 = 1/3.*A9*D[1, 1]+2/15.*A8*D[2, 2]
    C20 = 1/3.*A8*D[0, 0]+2/15.*A9*D[2, 2]
    C21 = D[0, 1]

    Keq = np.matrix(np.zeros((12, 12)))
    Keq[0, 0:13] = C1, C5, -C6, C2, C9, -C8, C4, C11, -C12, C3, C7, -C10
    Keq[1, 1:13] = C13, -C21, C9, C15, 0, -C11, C19, 0, -C7, C17, 0
    Keq[2, 2:13] = C14, C8, 0, C18, C12, 0, C20, -C10, 0, C16
    Keq[3, 3:13] = C1, C5, C6, C3, C7, C10, C4, C11, C12
    Keq[4, 4:13] = C13, C21, -C7, C17, 0, -C11, C19, 0
    Keq[5, 5:13] = C14, C10, 0, C16, -C12, 0, C20
    Keq[6, 6:13] = C1, -C5, C6, C2, -C9, C8
    Keq[7, 7:13] = C13, -C21, -C9, C15, 0
    Keq[8, 8:13] = C14, -C8, 0, C18
    Keq[9, 9:13] = C1, -C5, -C6
    Keq[10, 10:13] = C13, C21
    Keq[11, 11] = C14
    Keq = Keq.T+Keq-np.diag(np.diag(Keq))

    if eq != None:
        q = eq
        R1 = q*Lx*Ly/4
        R2 = q*Lx*Ly**2/24
        R3 = q*Ly*Lx**2/24

        feq = np.matrix([R1, R2, -R3, R1, R2, R3, R1, -R2, R3, R1, -R2, -R3])

    if eq != None:
        return Keq, feq
    else:
        return Keq


def planqe(ex, ey, ep, D, eq=None):
    """
    Calculate the stiffness matrix for a quadrilateral plane stress or plane strain element.

    Parameters
    ----------
    ex : array_like
        Element coordinates [x1, x2, x3, x4].
    ey : array_like
        Element coordinates [y1, y2, y3, y4].
    ep : array_like
        Element properties [ptype, t], where ptype is analysis type and t is element thickness.
    D : array_like
        Constitutive matrix.
    eq : array_like, optional
        Body force vector [bx, by], where bx, by are body forces in x, y directions.

    Returns
    -------
    Ke : ndarray
        Element stiffness matrix, shape (8, 8).
    fe : ndarray, optional
        Equivalent nodal forces, if eq is provided.
    """
    K = np.zeros((10, 10))
    f = np.zeros((10, 1))

    xm = sum(ex)/4.
    ym = sum(ey)/4.

    b1 = eq if eq is not None else np.array([[0], [0]])

    # Make sure b1 is properly shaped for plante
    b1 = np.asarray(b1).flatten()

    # Create element coordinates for triangular subelements
    ex1 = np.array([ex[0], ex[1], xm])
    ey1 = np.array([ey[0], ey[1], ym])
    ex2 = np.array([ex[1], ex[2], xm])
    ey2 = np.array([ey[1], ey[2], ym])
    ex3 = np.array([ex[2], ex[3], xm])
    ey3 = np.array([ey[2], ey[3], ym])
    ex4 = np.array([ex[3], ex[0], xm])
    ey4 = np.array([ey[3], ey[0], ym])

    # Create element matrices for each subelement
    ke1, fe1 = plante(ex1, ey1, ep, D, b1)
    # Convert fe1 to column vector if needed
    fe1 = np.asarray(fe1).reshape(-1, 1)
    K, f = assem(np.array([1, 2, 3, 4, 9, 10]), K, ke1, f, fe1)

    ke1, fe1 = plante(ex2, ey2, ep, D, b1)
    # Convert fe1 to column vector if needed
    fe1 = np.asarray(fe1).reshape(-1, 1)
    K, f = assem(np.array([3, 4, 5, 6, 9, 10]), K, ke1, f, fe1)

    ke1, fe1 = plante(ex3, ey3, ep, D, b1)
    # Convert fe1 to column vector if needed
    fe1 = np.asarray(fe1).reshape(-1, 1)
    K, f = assem(np.array([5, 6, 7, 8, 9, 10]), K, ke1, f, fe1)

    ke1, fe1 = plante(ex4, ey4, ep, D, b1)
    # Convert fe1 to column vector if needed
    fe1 = np.asarray(fe1).reshape(-1, 1)
    K, f = assem(np.array([7, 8, 1, 2, 9, 10]), K, ke1, f, fe1)

    # Static condensation
    Ke, fe = statcon(K, f, np.array([[9], [10]]))

    if eq is None:
        return Ke
    else:
        return Ke, fe

def planqs(ex, ey, ep, D, ed, eq=None):
    """
    Calculate element normal and shear stress for a quadrilateral plane stress or plane strain element.

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2, x3, x4].
    ey : array_like
        Element node y-coordinates [y1, y2, y3, y4].
    ep : array_like
        Element properties [ptype, t], where ptype is analysis type and t is thickness.
    D : array_like
        Constitutive matrix.
    ed : array_like
        Element displacement vector [u1, u2, ..., u8].
    eq : array_like, optional
        Body force vector [bx, by], where bx, by are body forces in x, y directions.

    Returns
    -------
    es : ndarray
        Element stress array [sigx, sigy, (sigz), tauxy].
    et : ndarray
        Element strain array [epsx, epsy, (epsz), gamxy].
    """

    # Convert inputs to arrays for consistency
    ex = np.asarray(ex).flatten()
    ey = np.asarray(ey).flatten()
    ed = np.asarray(ed).flatten()

    if len(ex) != 4 or len(ey) != 4 or len(ed) != 8:
        raise ValueError(
            'Error ! PLANQS: only one element at the time (ex, ey, ed must be arrays with 4, 4, and 8 elements)')

    K = np.zeros((10, 10))
    f = np.zeros((10, 1))

    xm = sum(ex)/4.
    ym = sum(ey)/4.

    b1 = eq if eq is not None else np.array([[0], [0]])

    ex1 = np.array([ex[0], ex[1], xm])
    ey1 = np.array([ey[0], ey[1], ym])
    ex2 = np.array([ex[1], ex[2], xm])
    ey2 = np.array([ey[1], ey[2], ym])
    ex3 = np.array([ex[2], ex[3], xm])
    ey3 = np.array([ey[2], ey[3], ym])
    ex4 = np.array([ex[3], ex[0], xm])
    ey4 = np.array([ey[3], ey[0], ym])

    ke1, fe1 = plante(ex1, ey1, ep, D, b1)
    K, f = assem(np.array([1, 2, 3, 4, 9, 10]), K, ke1, f, fe1)
    ke1, fe1 = plante(ex2, ey2, ep, D, b1)
    K, f = assem(np.array([3, 4, 5, 6, 9, 10]), K, ke1, f, fe1)
    ke1, fe1 = plante(ex3, ey3, ep, D, b1)
    K, f = assem(np.array([5, 6, 7, 8, 9, 10]), K, ke1, f, fe1)
    ke1, fe1 = plante(ex4, ey4, ep, D, b1)
    K, f = assem(np.array([7, 8, 1, 2, 9, 10]), K, ke1, f, fe1)

    A1 = 0.5 * np.linalg.det(np.vstack([
        np.ones(3),
        ex1,
        ey1
    ]).T)
    
    A2 = 0.5 * np.linalg.det(np.vstack([
        np.ones(3),
        ex2,
        ey2
    ]).T)
    
    A3 = 0.5 * np.linalg.det(np.vstack([
        np.ones(3),
        ex3,
        ey3
    ]).T)
    
    A4 = 0.5 * np.linalg.det(np.vstack([
        np.ones(3),
        ex4,
        ey4
    ]).T)
    
    Atot = A1+A2+A3+A4

    a, _ = solveq(K, f, np.array(range(1, 9)), ed)

    # Create proper element displacement arrays for each triangular sub-element
    ed1 = np.array([a[0,0], a[1,0], a[2,0], a[3,0], a[8,0], a[9,0]])
    ed2 = np.array([a[2,0], a[3,0], a[4,0], a[5,0], a[8,0], a[9,0]])
    ed3 = np.array([a[4,0], a[5,0], a[6,0], a[7,0], a[8,0], a[9,0]])
    ed4 = np.array([a[6,0], a[7,0], a[0,0], a[1,0], a[8,0], a[9,0]])

    s1, t1 = plants(ex1, ey1, ep, D, ed1)
    s2, t2 = plants(ex2, ey2, ep, D, ed2)
    s3, t3 = plants(ex3, ey3, ep, D, ed3)
    s4, t4 = plants(ex4, ey4, ep, D, ed4)

    es = (s1*A1+s2*A2+s3*A3+s4*A4)/Atot
    et = (t1*A1+t2*A2+t3*A3+t4*A4)/Atot

    return es[0], et[0]

def plani4e(ex, ey, ep, D, eq=None):
    """
    Calculate the stiffness matrix for a 4 node isoparametric element in plane strain or plane stress.
    
    Parameters
    ----------
    ex : array_like
        Element coordinates [x1, x2, x3, x4].
    ey : array_like
        Element coordinates [y1, y2, y3, y4].
    ep : array_like
        Element properties [ptype, t, ir], where ptype is analysis type, t is thickness, and ir is integration rule.
    D : array_like
        Constitutive matrix.
    eq : array_like, optional
        Body force vector [bx, by], where bx, by are body forces in x, y directions.
    
    Returns
    -------
    Ke : ndarray
        Element stiffness matrix, shape (8, 8).
    fe : ndarray, optional
        Equivalent nodal forces, shape (8, 1), if eq is provided.
    """
    ptype = ep[0]
    t = ep[1]
    ir = ep[2]
    ngp = ir*ir
    if eq is None:
        q = np.zeros((2, 1))
    else:
        q = np.reshape(eq, (2, 1))
#--------- gauss points --------------------------------------
    if ir == 1:
        g1 = 0.0
        w1 = 2.0
        gp = np.matrix([g1, g1])
        w = np.matrix([w1, w1])
    elif ir == 2:
        g1 = 0.577350269189626
        w1 = 1
        gp = np.matrix([
            [-g1, -g1],
            [g1, -g1],
            [-g1, g1],
            [g1, g1]])
        w = np.matrix([
            [w1, w1],
            [w1, w1],
            [w1, w1],
            [w1, w1]])
    elif ir == 3:
        g1 = 0.774596669241483
        g2 = 0.
        w1 = 0.555555555555555
        w2 = 0.888888888888888
        gp = np.matrix([
            [-g1, -g1],
            [-g2, -g1],
            [g1, -g1],
            [-g1, g2],
            [g2, g2],
            [g1, g2],
            [-g1, g1],
            [g2, g1],
            [g1, g1]])
        w = np.matrix([
            [w1, w1],
            [w2, w1],
            [w1, w1],
            [w1, w2],
            [w2, w2],
            [w1, w2],
            [w1, w1],
            [w2, w1],
            [w1, w1]])
    else:
        info("Used number of integrat     ion points not implemented")
    wp = np.multiply(w[:, 0], w[:, 1])
    xsi = gp[:, 0]
    eta = gp[:, 1]
    r2 = ngp*2
    # Shape Functions
    N = np.multiply((1-xsi), (1-eta))/4.
    N = np.append(N, np.multiply((1+xsi), (1-eta))/4., axis=1)
    N = np.append(N, np.multiply((1+xsi), (1+eta))/4., axis=1)
    N = np.append(N, np.multiply((1-xsi), (1+eta))/4., axis=1)

    dNr = np.matrix(np.zeros((r2, 4)))
    dNr[0:r2:2, 0] = -(1-eta)/4.
    dNr[0:r2:2, 1] = (1-eta)/4.
    dNr[0:r2:2, 2] = (1+eta)/4.
    dNr[0:r2:2, 3] = -(1+eta)/4.
    dNr[1:r2+1:2, 0] = -(1-xsi)/4.
    dNr[1:r2+1:2, 1] = -(1+xsi)/4.
    dNr[1:r2+1:2, 2] = (1+xsi)/4.
    dNr[1:r2+1:2, 3] = (1-xsi)/4.

#
    Ke1 = np.matrix(np.zeros((8, 8)))
    fe1 = np.matrix(np.zeros((8, 1)))
    JT = dNr*np.matrix([ex, ey]).T
    # --------- plane stress --------------------------------------
    if ptype == 1:
        colD = np.shape(D)[0]
        if colD > 3:
            Cm = np.linalg.inv(D)
            Dm = np.linalg.inv(Cm[np.ix_([0, 1, 3], [0, 1, 3])])
        else:
            Dm = D
#
        B = np.matrix(np.zeros((3, 8)))
        N2 = np.matrix(np.zeros((2, 8)))
        for i in range(ngp):
            indx = np.array([2*(i+1)-1, 2*(i+1)])
            detJ = np.linalg.det(JT[indx-1, :])
            if detJ < 10*np.finfo(float).eps:
                info("Jacobi determinant equal or less than zero!")
            JTinv = np.linalg.inv(JT[indx-1, :])
            dNx = JTinv*dNr[indx-1, :]
#
            index_array_even = np.array([0, 2, 4, 6])
            index_array_odd = np.array([1, 3, 5, 7])
#
            counter = 0
            for index in index_array_even:
                B[0, index] = dNx[0, counter]
                B[2, index] = dNx[1, counter]
                N2[0, index] = N[i, counter]
                counter = counter+1
#
            counter = 0
            for index in index_array_odd:
                B[1, index] = dNx[1, counter]
                B[2, index] = dNx[0, counter]
                N2[1, index] = N[i, counter]
                counter = counter+1
#
            Ke1 = Ke1+B.T*Dm*B*detJ*wp[i].item()*t
            fe1 = fe1 + N2.T * q * detJ * wp[i].item() * t

        return Ke1, fe1
#--------- plane strain --------------------------------------
    elif ptype == 2:
        #
        colD = np.shape(D)[0]
        if colD > 3:
            Dm = D[np.ix_([0, 1, 3], [0, 1, 3])]
        else:
            Dm = D
#
        B = np.matrix(np.zeros((3, 8)))
        N2 = np.matrix(np.zeros((2, 8)))
        for i in range(ngp):
            indx = np.array([2*(i+1)-1, 2*(i+1)])
            detJ = np.linalg.det(JT[indx-1, :])
            if detJ < 10*np.finfo(float).eps:
                info("Jacobideterminant equal or less than zero!")
            JTinv = np.linalg.inv(JT[indx-1, :])
            dNx = JTinv*dNr[indx-1, :]
#
            index_array_even = np.array([0, 2, 4, 6])
            index_array_odd = np.array([1, 3, 5, 7])
#
            counter = 0
            for index in index_array_even:
                #
                B[0, index] = dNx[0, counter]
                B[2, index] = dNx[1, counter]
                N2[0, index] = N[i, counter]
#
                counter = counter+1
#
            counter = 0
            for index in index_array_odd:
                B[1, index] = dNx[1, counter]
                B[2, index] = dNx[0, counter]
                N2[1, index] = N[i, counter]
                counter = counter+1
#
            Ke1 = Ke1 + B.T * Dm * B * detJ * wp[i].item() * t
            fe1 = fe1+N2.T*q*detJ*wp[i].item()*t
        return Ke1, fe1
    else:
        info("Error ! Check first argument, ptype=1 or 2 allowed")


def soli8e(ex, ey, ez, ep, D, eqp=None):
    """
    Calculate the stiffness matrix (and optionally load vector) for an 8-node (brick) isoparametric element.

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2, x3, ..., x8].
    ey : array_like
        Element node y-coordinates [y1, y2, y3, ..., y8].
    ez : array_like
        Element node z-coordinates [z1, z2, z3, ..., z8].
    ep : array_like
        Element properties [ir], where ir is the integration rule.
    D : array_like
        Constitutive matrix.
    eqp : array_like, optional
        Body force vector [bx, by, bz], where bx, by, bz are body forces in x, y, z directions.

    Returns
    -------
    Ke : ndarray
        Element stiffness matrix.
    fe : ndarray, optional
        Equivalent nodal forces, if eqp is not None.

    History
    -------
    LAST MODIFIED: M Ristinmaa   1995-10-25
                   J Lindemann   2022-01-24 (Python version)
    """
    ir = ep[0]
    ngp = ir*ir*ir

    if eqp is None:
        eq = np.zeros((3, 1))
    else:
        eq = eqp

    if ir == 1:
        g1 = 0.0
        w1 = 2.0
        gp = np.array([g1, g1, g1]).reshape(1, 3)
        w = np.array([w1, w1, w1]).reshape(1, 3)
    elif ir == 2:
        g1 = 0.577350269189626
        w1 = 1
        gp = np.zeros((8, 3))
        w = np.zeros((8, 3))
        gp[:, 0] = np.array([-1, 1, 1, -1, -1, 1, 1, -1])*g1
        w[:, 0] = np.array([1, 1, 1, 1, 1, 1, 1, 1])*w1
        gp[:, 1] = np.array([-1, -1, 1, 1, -1, -1, 1, 1])*g1
        w[:, 1] = np.array([1, 1, 1, 1, 1, 1, 1, 1])*w1
        gp[:, 2] = np.array([-1, -1, -1, -1, 1, 1, 1, 1])*g1
        w[:, 2] = np.array([1, 1, 1, 1, 1, 1, 1, 1])*w1
    else:
        g1 = 0.774596669241483,
        g2 = 0.0
        w1 = 0.555555555555555
        w2 = 0.888888888888888

        gp = np.zeros((27, 3))
        w = np.zeros((27, 3))

        I1 = np.array([-1, 0, 1, -1, 0, 1, -1, 0, 1]).reshape(1, 9)
        I2 = np.array([0, -1, 0, 0, 1, 0, 0, 1, 0]).reshape(1, 9)

        gp[:, 0] = np.concatenate((I1, I1, I1), axis=1)*g1
        gp[:, 0] = np.concatenate((I2, I2, I2), axis=1)*g2 + gp[:, 0]

        I1 = np.abs(I1)
        I2 = np.abs(I2)

        w[:, 0] = np.concatenate((I1, I1, I1), axis=1)*w1
        w[:, 0] = np.concatenate((I2, I2, I2), axis=1)*w2 + w[:, 0]

        I1 = np.array([-1, -1, -1, 0, 0, 0, 1, 1, 1]).reshape(1, 9)
        I2 = np.array([0, 0, 0, 1, 1, 1, 0, 0, 0]).reshape(1, 9)

        gp[:, 1] = np.concatenate((I1, I1, I1), axis=1)*g1
        gp[:, 1] = np.concatenate((I2, I2, I2), axis=1)*g2 + gp[:, 1]

        I1 = np.abs(I1)
        I2 = np.abs(I2)

        w[:, 1] = np.concatenate((I1, I1, I1), axis=1)*w1
        w[:, 1] = np.concatenate((I2, I2, I2), axis=1)*w2 + w[:, 1]

        I1 = np.array([-1, -1, -1, -1, -1, -1, -1, -1, -1]).reshape(1, 9)
        I2 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0]).reshape(1, 9)
        I3 = np.abs(I1)

        gp[:, 2] = np.concatenate((I1, I2, I3), axis=1)*g1
        gp[:, 2] = np.concatenate((I2, I3, I2), axis=1)*g2 + gp[:, 2]

        w[:, 2] = np.concatenate((I3, I2, I3), axis=1)*w1
        w[:, 2] = np.concatenate((I2, I3, I2), axis=1)*w2 + w[:, 2]

    wp = w[:, 0]*w[:, 1]*w[:, 2]

    xsi = gp[:, 0]
    eta = gp[:, 1]
    zet = gp[:, 2]
    r2 = ngp*3

    N = np.zeros((ngp, 8))
    dNr = np.zeros((r2, 8))

    N[:, 0] = (1-xsi)*(1-eta)*(1-zet)/8
    N[:, 1] = (1+xsi)*(1-eta)*(1-zet)/8
    N[:, 2] = (1+xsi)*(1+eta)*(1-zet)/8
    N[:, 3] = (1-xsi)*(1+eta)*(1-zet)/8
    N[:, 4] = (1-xsi)*(1-eta)*(1+zet)/8
    N[:, 5] = (1+xsi)*(1-eta)*(1+zet)/8
    N[:, 6] = (1+xsi)*(1+eta)*(1+zet)/8
    N[:, 7] = (1-xsi)*(1+eta)*(1+zet)/8

    dNr[0:r2+1:3, 0] = -(1-eta)*(1-zet)
    dNr[0:r2+1:3, 1] = (1-eta)*(1-zet)
    dNr[0:r2+1:3, 2] = (1+eta)*(1-zet)
    dNr[0:r2+1:3, 3] = -(1+eta)*(1-zet)
    dNr[0:r2+1:3, 4] = -(1-eta)*(1+zet)
    dNr[0:r2+1:3, 5] = (1-eta)*(1+zet)
    dNr[0:r2+1:3, 6] = (1+eta)*(1+zet)
    dNr[0:r2+1:3, 7] = -(1+eta)*(1+zet)
    dNr[1:r2+2:3, 0] = -(1-xsi)*(1-zet)
    dNr[1:r2+2:3, 1] = -(1+xsi)*(1-zet)
    dNr[1:r2+2:3, 2] = (1+xsi)*(1-zet)
    dNr[1:r2+2:3, 3] = (1-xsi)*(1-zet)
    dNr[1:r2+2:3, 4] = -(1-xsi)*(1+zet)
    dNr[1:r2+2:3, 5] = -(1+xsi)*(1+zet)
    dNr[1:r2+2:3, 6] = (1+xsi)*(1+zet)
    dNr[1:r2+2:3, 7] = (1-xsi)*(1+zet)
    dNr[2:r2+3:3, 0] = -(1-xsi)*(1-eta)
    dNr[2:r2+3:3, 1] = -(1+xsi)*(1-eta)
    dNr[2:r2+3:3, 2] = -(1+xsi)*(1+eta)
    dNr[2:r2+3:3, 3] = -(1-xsi)*(1+eta)
    dNr[2:r2+3:3, 4] = (1-xsi)*(1-eta)
    dNr[2:r2+3:3, 5] = (1+xsi)*(1-eta)
    dNr[2:r2+3:3, 6] = (1+xsi)*(1+eta)
    dNr[2:r2+3:3, 7] = (1-xsi)*(1+eta)

    dNr = dNr/8.0

    Ke = np.zeros((24, 24))
    fe = np.zeros((24, 1))

    ex = np.asarray(ex).reshape((8, 1))
    ey = np.asarray(ey).reshape((8, 1))
    ez = np.asarray(ez).reshape((8, 1))

    JT = dNr@np.concatenate((ex, ey, ez), axis=1)

    eps = np.finfo(float).eps

    for i in range(ngp):
        indx = [i*3, i*3+1, i*3+2]
        detJ = np.linalg.det(JT[indx, :])
        if detJ < 10*eps:
            print('Jacobideterminant equal or less than zero!')
        JTinv = np.linalg.inv(JT[indx, :])
        dNx = JTinv@dNr[indx, :]

        B = np.zeros((6, 24))
        N2 = np.zeros((3, 24))

        B[0, 0:24:3] = dNx[0, :]
        B[1, 1:25:3] = dNx[1, :]
        B[2, 2:26:3] = dNx[2, :]
        B[3, 0:24:3] = dNx[1, :]
        B[3, 1:25:3] = dNx[0, :]
        B[4, 0:24:3] = dNx[2, :]
        B[4, 2:26:3] = dNx[0, :]
        B[5, 1:25:3] = dNx[2, :]
        B[5, 2:26:3] = dNx[1, :]

        N2[0, 0:24:3] = N[i, :]
        N2[1, 1:25:3] = N[i, :]
        N2[2, 2:26:3] = N[i, :]

        Ke = Ke + (np.transpose(B)@D@B)*detJ*wp[i]
        fe = fe + (np.transpose(N2)@eq)*detJ*wp[i]

    if eqp != None:
        return Ke, fe
    else:
        return Ke


def soli8s(ex, ey, ez, ep, D, ed):
    """
    Calculate element normal and shear stress for an 8-node (brick) isoparametric element.

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates [x1, x2, x3, ..., x8].
    ey : array_like
        Element node y-coordinates [y1, y2, y3, ..., y8].
    ez : array_like
        Element node z-coordinates [z1, z2, z3, ..., z8].
    ep : array_like
        Element properties [Ir], where Ir is the integration rule.
    D : array_like
        Constitutive matrix.
    ed : array_like
        Element displacement vector [u1, u2, ..., u24].

    Returns
    -------
    es : ndarray
        Element stress matrix, one row for each integration point.
        Each row contains [sigx, sigy, sigz, sigxy, sigyz, sigxz].
    et : ndarray
        Element strain matrix, one row for each integration point.
        Each row contains [epsx, epsy, epsz, epsxy, epsyz, epsxz].

    History
    -------
    LAST MODIFIED: M Ristinmaa   1995-10-25
                   J Lindemann   2022-02-23 (Python version)
    """

    ir = ep[0]
    ngp = ir*ir*ir

    ir = ep[0]
    ngp = ir*ir*ir

    if ir == 1:
        g1 = 0.0
        w1 = 2.0
        gp = np.array([g1, g1, g1]).reshape(1, 3)
        w = np.array([w1, w1, w1]).reshape(1, 3)
    elif ir == 2:
        g1 = 0.577350269189626
        w1 = 1
        gp = np.zeros((8, 3))
        w = np.zeros((8, 3))
        gp[:, 0] = np.array([-1, 1, 1, -1, -1, 1, 1, -1])*g1
        w[:, 0] = np.array([1, 1, 1, 1, 1, 1, 1, 1])*w1
        gp[:, 1] = np.array([-1, -1, 1, 1, -1, -1, 1, 1])*g1
        w[:, 1] = np.array([1, 1, 1, 1, 1, 1, 1, 1])*w1
        gp[:, 2] = np.array([-1, -1, -1, -1, 1, 1, 1, 1])*g1
        w[:, 2] = np.array([1, 1, 1, 1, 1, 1, 1, 1])*w1
    else:
        g1 = 0.774596669241483,
        g2 = 0.0
        w1 = 0.555555555555555
        w2 = 0.888888888888888

        gp = np.zeros((27, 3))
        w = np.zeros((27, 3))

        I1 = np.array([-1, 0, 1, -1, 0, 1, -1, 0, 1]).reshape(1, 9)
        I2 = np.array([0, -1, 0, 0, 1, 0, 0, 1, 0]).reshape(1, 9)

        gp[:, 0] = np.concatenate((I1, I1, I1), axis=1)*g1
        gp[:, 0] = np.concatenate((I2, I2, I2), axis=1)*g2 + gp[:, 0]

        I1 = np.abs(I1)
        I2 = np.abs(I2)

        w[:, 0] = np.concatenate((I1, I1, I1), axis=1)*w1
        w[:, 0] = np.concatenate((I2, I2, I2), axis=1)*w2 + w[:, 0]

        I1 = np.array([-1, -1, -1, 0, 0, 0, 1, 1, 1]).reshape(1, 9)
        I2 = np.array([0, 0, 0, 1, 1, 1, 0, 0, 0]).reshape(1, 9)

        gp[:, 1] = np.concatenate((I1, I1, I1), axis=1)*g1
        gp[:, 1] = np.concatenate((I2, I2, I2), axis=1)*g2 + gp[:, 1]

        I1 = np.abs(I1)
        I2 = np.abs(I2)

        w[:, 1] = np.concatenate((I1, I1, I1), axis=1)*w1
        w[:, 1] = np.concatenate((I2, I2, I2), axis=1)*w2 + w[:, 1]

        I1 = np.array([-1, -1, -1, -1, -1, -1, -1, -1, -1]).reshape(1, 9)
        I2 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0]).reshape(1, 9)
        I3 = np.abs(I1)

        gp[:, 2] = np.concatenate((I1, I2, I3), axis=1)*g1
        gp[:, 2] = np.concatenate((I2, I3, I2), axis=1)*g2 + gp[:, 2]

        w[:, 2] = np.concatenate((I3, I2, I3), axis=1)*w1
        w[:, 2] = np.concatenate((I2, I3, I2), axis=1)*w2 + w[:, 2]

    wp = w[:, 0]*w[:, 1]*w[:, 2]

    xsi = gp[:, 0]
    eta = gp[:, 1]
    zet = gp[:, 2]
    r2 = ngp*3

    N = np.zeros((ngp, 8))
    dNr = np.zeros((r2, 8))

    N[:, 0] = (1-xsi)*(1-eta)*(1-zet)/8
    N[:, 1] = (1+xsi)*(1-eta)*(1-zet)/8
    N[:, 2] = (1+xsi)*(1+eta)*(1-zet)/8
    N[:, 3] = (1-xsi)*(1+eta)*(1-zet)/8
    N[:, 4] = (1-xsi)*(1-eta)*(1+zet)/8
    N[:, 5] = (1+xsi)*(1-eta)*(1+zet)/8
    N[:, 6] = (1+xsi)*(1+eta)*(1+zet)/8
    N[:, 7] = (1-xsi)*(1+eta)*(1+zet)/8

    dNr[0:r2+1:3, 0] = -(1-eta)*(1-zet)
    dNr[0:r2+1:3, 1] = (1-eta)*(1-zet)
    dNr[0:r2+1:3, 2] = (1+eta)*(1-zet)
    dNr[0:r2+1:3, 3] = -(1+eta)*(1-zet)
    dNr[0:r2+1:3, 4] = -(1-eta)*(1+zet)
    dNr[0:r2+1:3, 5] = (1-eta)*(1+zet)
    dNr[0:r2+1:3, 6] = (1+eta)*(1+zet)
    dNr[0:r2+1:3, 7] = -(1+eta)*(1+zet)
    dNr[1:r2+2:3, 0] = -(1-xsi)*(1-zet)
    dNr[1:r2+2:3, 1] = -(1+xsi)*(1-zet)
    dNr[1:r2+2:3, 2] = (1+xsi)*(1-zet)
    dNr[1:r2+2:3, 3] = (1-xsi)*(1-zet)
    dNr[1:r2+2:3, 4] = -(1-xsi)*(1+zet)
    dNr[1:r2+2:3, 5] = -(1+xsi)*(1+zet)
    dNr[1:r2+2:3, 6] = (1+xsi)*(1+zet)
    dNr[1:r2+2:3, 7] = (1-xsi)*(1+zet)
    dNr[2:r2+3:3, 0] = -(1-xsi)*(1-eta)
    dNr[2:r2+3:3, 1] = -(1+xsi)*(1-eta)
    dNr[2:r2+3:3, 2] = -(1+xsi)*(1+eta)
    dNr[2:r2+3:3, 3] = -(1-xsi)*(1+eta)
    dNr[2:r2+3:3, 4] = (1-xsi)*(1-eta)
    dNr[2:r2+3:3, 5] = (1+xsi)*(1-eta)
    dNr[2:r2+3:3, 6] = (1+xsi)*(1+eta)
    dNr[2:r2+3:3, 7] = (1-xsi)*(1+eta)

    dNr = dNr/8.0

    ex = np.asarray(ex).reshape((8, 1))
    ey = np.asarray(ey).reshape((8, 1))
    ez = np.asarray(ez).reshape((8, 1))

    JT = dNr@np.concatenate((ex, ey, ez), axis=1)

    eps = np.finfo(float).eps
    
    eci = N@np.concatenate((ex, ey, ez), axis=1)
    et = np.zeros((ngp, 6))
    es = np.zeros((ngp, 6))

    ed = ed.reshape(1, 24)

    for i in range(ngp):
        indx = [i*3, i*3+1, i*3+2]
        detJ = np.linalg.det(JT[indx, :])
        if detJ < 10*eps:
            print('Jacobideterminant equal or less than zero!')
        JTinv = np.linalg.inv(JT[indx, :])
        dNx = JTinv@dNr[indx, :]

        B = np.zeros((6, 24))
        N2 = np.zeros((3, 24))

        B[0, 0:24:3] = dNx[0, :]
        B[1, 1:25:3] = dNx[1, :]
        B[2, 2:26:3] = dNx[2, :]
        B[3, 0:24:3] = dNx[1, :]
        B[3, 1:25:3] = dNx[0, :]
        B[4, 0:24:3] = dNx[2, :]
        B[4, 2:26:3] = dNx[0, :]
        B[5, 1:25:3] = dNx[2, :]
        B[5, 2:26:3] = dNx[1, :]

        N2[0, 0:24:3] = N[i, :]
        N2[1, 1:25:3] = N[i, :]
        N2[2, 2:26:3] = N[i, :]

        # [6x24] x [24,1]
        ee = B@np.transpose(ed)

        et[i, :] = ee.reshape(6,)
        es[i, :] = (D@ee).reshape(6,)

    return et, es, eci


def assem(edof: ArrayLike, K: Union[NDArray[np.floating], csr_matrix, csc_matrix, lil_matrix], Ke: ArrayLike, f: Optional[ArrayLike] = None, fe: Optional[ArrayLike] = None) -> Union[NDArray[np.floating], Tuple[NDArray[np.floating], NDArray[np.floating]]]:
    """
    Assemble element matrices Ke (and fe) into the global stiffness matrix K (and global force vector f).

    Parameters
    ----------
    edof : array_like
        DOF topology array.
    K : array_like
        The global stiffness matrix.
    Ke : array_like
        Element stiffness matrix.
    f : array_like, optional
        The global force vector.
    fe : array_like, optional
        Element force vector.

    Returns
    -------
    K : ndarray
        The updated global stiffness matrix.
    f : ndarray, optional
        The updated global force vector, if f and fe are provided.
    """

    import scipy.sparse as sp
    import numpy as np
    
    # Handle sparse matrices case
    if isinstance(K, sp.lil_matrix) or isinstance(K, sp.csr_matrix) or isinstance(K, sp.csc_matrix):
        if edof.ndim == 1:
            idx = edof-1
            # Convert to COO format for efficient modification
            if not isinstance(K, sp.lil_matrix):
                K = K.tolil()
                
            # Manually update each element
            for i in range(len(idx)):
                for j in range(len(idx)):
                    K[idx[i], idx[j]] += Ke[i, j]
            
            if (not f is None) and (not fe is None):
                fe_array = np.asarray(fe)
                if fe_array.ndim == 2 and fe_array.shape[1] > 1:
                    fe_shaped = fe_array.reshape(-1, 1)
                else:
                    fe_shaped = fe_array
                f[idx] += fe_shaped.flatten()
        else:
            for row in edof:
                idx = row-1
                # Convert to LIL format for efficient modification
                if not isinstance(K, sp.lil_matrix):
                    K = K.tolil()
                    
                # Manually update each element    
                for i in range(len(idx)):
                    for j in range(len(idx)):
                        K[idx[i], idx[j]] += Ke[i, j]
                
                if (not f is None) and (not fe is None):
                    fe_array = np.asarray(fe)
                    if fe_array.ndim == 2 and fe_array.shape[1] > 1:
                        fe_shaped = fe_array.reshape(-1, 1)
                    else:
                        fe_shaped = fe_array
                    f[idx] += fe_shaped.flatten()
    else:
        # Original code for dense matrices
        if edof.ndim == 1:
            idx = edof-1
            K[np.ix_(idx, idx)] = K[np.ix_(idx, idx)] + Ke
            
            if (not f is None) and (not fe is None):
                fe_array = np.asarray(fe)
                if fe_array.ndim == 2 and fe_array.shape[1] > 1:
                    fe_shaped = fe_array.reshape(-1, 1)
                else:
                    fe_shaped = fe_array
                f[np.ix_(idx)] = f[np.ix_(idx)] + fe_shaped
        else:
            for row in edof:
                idx = row-1
                K[np.ix_(idx, idx)] = K[np.ix_(idx, idx)] + Ke
                if (not f is None) and (not fe is None):
                    fe_array = np.asarray(fe)
                    if fe_array.ndim == 2 and fe_array.shape[1] > 1:
                        fe_shaped = fe_array.reshape(-1, 1)
                    else:
                        fe_shaped = fe_array
                    f[np.ix_(idx)] = f[np.ix_(idx)] + fe_shaped

    if f is None:
        return K
    else:
        return K, f
    
def solveq(K, f, bcPrescr=None, bcVal=None):
    """
    Solve static FE-equations considering boundary conditions.
    
    Parameters
    ----------
    K : array_like
        Global stiffness matrix, shape (nd, nd).
    f : array_like
        Global load vector, shape (nd, 1).
    bcPrescr : array_like
        1-dim integer array containing prescribed dofs.
    bcVal : array_like, optional
        1-dim float array containing prescribed values.
        If not given all prescribed dofs are assumed 0.
        
    Returns
    -------
    a : ndarray
        Solution including boundary values, shape (nd, 1).
    Q : ndarray
        Reaction force vector, shape (nd, 1).
    """

    nDofs = K.shape[0]
    nPdofs = bcPrescr.shape[0]

    if bcVal is None:
        bcVal = np.zeros([nPdofs], 'd')

    if bcPrescr is None:
        return np.asmatrix(np.linalg.solve(K, f))

    bc = np.ones(nDofs, 'bool')
    bcDofs = np.arange(nDofs)

    bc[np.ix_(bcPrescr-1)] = False
    bcDofs = bcDofs[bc]

    fsys = f[bcDofs]-K[np.ix_((bcDofs), (bcPrescr-1))] * \
        np.asmatrix(bcVal).reshape(nPdofs, 1)
    asys = np.linalg.solve(K[np.ix_((bcDofs), (bcDofs))], fsys)

    a = np.zeros([nDofs, 1])
    a[np.ix_(bcPrescr-1)] = np.asmatrix(bcVal).reshape(nPdofs, 1)
    a[np.ix_(bcDofs)] = asys

    Q = K*np.asmatrix(a)-f

    return (np.asmatrix(a), Q)


def spsolveq(K, f, bcPrescr, bcVal=None):
    """
    Solve static FE-equations considering boundary conditions.
    
    Parameters
    ----------
    K : array_like
        Global stiffness matrix, shape (nd, nd).
    f : array_like
        Global load vector, shape (nd, 1).
    bcPrescr : array_like
        1-dim integer array containing prescribed dofs.
    bcVal : array_like, optional
        1-dim float array containing prescribed values.
        If not given all prescribed dofs are assumed 0.
        
    Returns
    -------
    a : ndarray
        Solution including boundary values, shape (nd, 1).
    Q : ndarray
        Reaction force vector, shape (nd, 1).
    """

    nDofs = K.shape[0]
    nPdofs = bcPrescr.shape[0]

    if bcVal is None:
        bcVal = np.zeros([nPdofs], 'd')

    bc = np.ones(nDofs, 'bool')
    bcDofs = np.arange(nDofs)

    bc[np.ix_(bcPrescr-1)] = False
    bcDofs = bcDofs[bc]

    bcVal_m = np.asmatrix(bcVal).reshape(nPdofs, 1)

    info("Preparing system matrix...")

    mask = np.ones(K.shape[0], dtype=bool)
    mask[bcDofs] = False

    info("step 1... converting K->CSR")
    Kcsr = K.asformat("csr")
    info("step 2... Kt")
    #Kt1 = K[bcDofs]
    #Kt = Kt1[:,bcPrescr]
    Kt = K[np.ix_((bcDofs), (bcPrescr-1))]
    info("step 3... fsys")
    fsys = f[bcDofs]-Kt*bcVal_m
    info("step 4... Ksys")
    Ksys1 = Kcsr[bcDofs]
    Ksys = Ksys1[:, bcDofs]
    #Ksys = Kcsr[np.ix_((bcDofs),(bcDofs))]
    info("done...")

    info("Solving system...")
    asys = dsolve.spsolve(Ksys, fsys)

    info("Reconstructing full a...")
    a = np.zeros([nDofs, 1])
    a[np.ix_(bcPrescr-1)] = bcVal_m
    a[np.ix_(bcDofs)] = np.asmatrix(asys).transpose()

    a_m = np.asmatrix(a)
    Q = K*a_m-f
    info("done...")
    return (a_m, Q)


def eigen(K: ArrayLike, M: ArrayLike, b: Optional[ArrayLike] = None) -> Tuple[NDArray[np.floating], NDArray[np.floating]]:
    """
    Solve the generalized eigenvalue problem |K-LM|X = 0, considering boundary conditions.

    Parameters
    ----------
    K : array_like
        Global stiffness matrix, shape (ndof, ndof).
    M : array_like
        Global mass matrix, shape (ndof, ndof).
    b : array_like, optional
        Boundary condition vector, shape (nbc, 1).

    Returns
    -------
    L : ndarray
        Eigenvalue vector, shape (ndof-nbc, 1).
    X : ndarray
        Eigenvectors, shape (ndof, ndof-nbc).
    """
    nd, _ = K.shape
    if b is not None:
        fdof = np.setdiff1d(np.arange(nd), b-1)
        D, X1 = eig(K[np.ix_(fdof,fdof)], M[np.ix_(fdof,fdof)])
        D = np.real(D)
        nfdof, _ = X1.shape
        for j in range(nfdof):
            mnorm = np.sqrt(X1[:,j].T@M[np.ix_(fdof,fdof)]@X1[:,j])
            X1[:,j] /= mnorm
        s_order = np.argsort(D)
        L = np.sort(D)
        X2 = np.zeros(X1.shape)
        for ind,j in enumerate(s_order):
            X2[:,ind] = X1[:,j]
        X = np.zeros((nd,nfdof))
        X[fdof,:] = X2
        return L, X
    else:
        D, X1 = eig(K, M)
        D = np.real(D)
        for j in range(nd):
            mnorm = np.sqrt(X1[:,j].T@M@X1[:,j])
            X1[:,j] /= mnorm
        s_order = np.argsort(D)
        L = np.sort(D)
        X = np.copy(X1)
        for ind,j in enumerate(s_order):
            X[:,ind] = X1[:,j]
        return L, X


def gfunc(G,dt):
    """
    Form vector with function values at equally spaced points by linear interpolation.

    Parameters
    ----------
    G : array_like
        Time-function pairs [t_i, g_i], where t_i is time i and g_i is g(t_i), shape (np, 2).
    dt : float
        Time step.

    Returns
    -------
    t : ndarray
        1-D vector with equally spaced time points.
    g : ndarray
        1-D vector with corresponding function values.
    """
    ti = np.arange(G[0,0],G[-1,0]+dt,dt)
    g1 = np.interp(ti,G[:,0],G[:,1])
    return ti, g1


def step1(K,C,f,a0,bc,ip,times,dofs):
    """
    Algorithm for dynamic solution of first-order FE equations considering boundary conditions.

    Parameters
    ----------
    K : array_like
        Conductivity matrix, shape (ndof, ndof).
    C : array_like
        Capacity matrix, shape (ndof, ndof).
    f : array_like
        Load vector, shape (ndof, nstep + 1).
        If shape (ndof, 1), the values are kept constant during time integration.
    a0 : array_like
        Initial vector a(0), shape (ndof, 1).
    bc : array_like
        Boundary condition matrix, shape (nbc, nstep + 2).
        where nbc = number of prescribed degrees of freedom (either constant or time-dependent).
        The first column contains the numbers of the prescribed degrees of freedom
        and the subsequent columns contain the time history.
        If shape (nbc, 2), the values from the second column are kept constant
        during time integration.
    ip : array_like
        Array [dt, tottime, alpha], where
        dt is the size of the time increment,
        tottime is the total time,
        alpha is time integration constant.
        Frequently used values of alpha are:
        alpha=0: forward difference; forward Euler,
        alpha=1/2: trapezoidal rule; Crank-Nicholson
        alpha=1: backward difference; backward Euler
    times : array_like
        Array [t(i) ...] of times at which output should be written to a and da.
    dofs : array_like
        Array [dof(i) ...] of degree of freedom numbers for which history output
        should be written to ahist and dahist.

    Returns
    -------
    modelhist : dict
        Dictionary containing solution history for the whole model with keys:
        
        - 'a' : ndarray
            Values of a at all timesteps, alternatively at times specified in 'times',
            shape (ndof, nstep + 1) or (ndof, ntimes).
        - 'da' : ndarray
            Values of da at all timesteps, alternatively at times specified in 'times',
            shape (ndof, nstep + 1) or (ndof, ntimes).
    dofhist : dict
        Dictionary containing solution history for the degrees of freedom selected in 'dofs' with keys:
        
        - 'a' : ndarray
            Time history of a at the dofs specified in 'dofs',
            shape (ndof, nstep + 1).
        - 'da' : ndarray
            Time history of da at the dofs specified in 'dofs',
            shape (ndof, nstep + 1).
    """
    ndof, _ = K.shape
    dt, tottime, alpha = ip
    a1 = (1-alpha)*dt
    a2 = alpha*dt

    nstep = 1
    if np.array(f).any():
        _, ncf = f.shape
        if ncf>1:
            nstep = ncf-1

    if np.array(bc).any():
        _, ncb = bc.shape
        if ncb>2:
            nstep = ncb-2
        bound = 1
    if not np.array(bc).any():
        bound = 0

    ns = int(tottime/dt)
    if (ns < nstep or nstep==1):
        nstep=ns

    tf = np.zeros((ndof,nstep+1))
    if np.array(f).any():
        if ncf==1:
            tf = f[:,0].reshape(-1,1)@np.ones((1,nstep+1))
        if ncf>1:
            tf = np.copy(f)

    modelhist = {}
    sa=0
    if not np.array(times).any():
        ntimes=0
        sa=1
        modelhist['a'] = np.zeros((ndof,nstep+1))
        modelhist['da'] = np.zeros((ndof,nstep+1))
    else:
        ntimes = len(times)
        if ntimes:
            sa=2
            modelhist['a'] = np.zeros((ndof,ntimes))
            modelhist['da'] = np.zeros((ndof,ntimes))

    dofhist = {}
    if np.array(dofs).all():
        ndofs = len(dofs)
        if ndofs:
            dofhist['a'] = np.zeros((ndofs,nstep+1))
            dofhist['da'] = np.zeros((ndofs,nstep+1))
    else:
        ndofs=0

    itime = 0

    # Calculate initial second time derivative d2a0
    da0 = np.linalg.solve(C,tf[:,0].reshape(-1,1) - K@a0)
    # Save initial values
    if sa==1:
        modelhist['a'][:,0] = a0.ravel()
        modelhist['da'][:,0] = da0.ravel()
    elif sa==2:
        if times[itime]==0:
            modelhist['a'][:,itime] = a0.ravel()
            modelhist['da'][:,itime] = da0.ravel()
            itime += 1

    if ndofs:
        dofhist['a'][:,0] = a0[np.ix_(dofs-1)].ravel()
        dofhist['da'][:,0] = da0[np.ix_(dofs-1)].ravel()

    # Reduce matrices due to bcs
    tempa = np.zeros((ndof,1))
    tempda = np.zeros((ndof,1))
    fdof=np.arange(1,ndof+1).astype(int)
    if bound:
        nrb, ncb = bc.shape
        if ncb==2:
            pa = bc[:,1].reshape(-1,1)@np.ones((1,nstep+1))
            pda = np.zeros((nrb,nstep+1))
        elif ncb>2:
            pa = np.copy(bc[:,1:])
            pda1 = (pa[:,1]-pa[:,0])/dt
            pdarest = (pa[:,1:] - pa[:,0:-1])/dt
            pda = np.hstack((pda1.reshape(-1,1),pdarest))
        pdof = np.copy(bc[:,0]).astype(int)
        fdof = np.setdiff1d(fdof,pdof).astype(int) - 1
        pdof -= 1 #adjusting for indexing starting from 0
        Keff = C[np.ix_(fdof,fdof)] + a2*K[np.ix_(fdof,fdof)]
    else:
        fdof -= 1 #adjusting for indexing starting from 0
        Keff = C + a2*K

    L, U = lu(Keff,permute_l=True)
    anew = a0[np.ix_(fdof)]
    danew = da0[np.ix_(fdof)]

    # Iterate over time steps
    for j in range(1,nstep+1):
        time = dt*j
        aold = np.copy(anew)
        daold = np.copy(danew)
        apred = aold + a1*daold
        if not bound:
            reff = tf[:,j].reshape(-1,1) - K@apred
        else:
            pdeff = C[np.ix_(fdof,pdof)]@pda[:,j].reshape(-1,1) + K[np.ix_(fdof,pdof)]@pa[:,j].reshape(-1,1)
            reff = tf[np.ix_(fdof),j].reshape(-1,1) - K[np.ix_(fdof,fdof)]@apred - pdeff
        y = np.linalg.solve(L,reff)
        danew = np.linalg.solve(U,y)
        anew = apred + a2*danew
        # Save to modelhist and dofhist
        if bound:
            tempa[np.ix_(pdof)] = pa[:,j].reshape(-1,1)
            tempda[np.ix_(pdof)] = pda[:,j].reshape(-1,1)
        tempa[np.ix_(fdof)] = anew
        tempda[np.ix_(fdof)] = danew
        if sa==1:
            modelhist['a'][:,j] = tempa.ravel()
            modelhist['da'][:,j] = tempda.ravel()
        elif sa==2:
            if ntimes and itime < ntimes:
                if time >= times[itime]:
                    modelhist['a'][:,itime] = tempa.ravel()
                    modelhist['da'][:,itime] = tempda.ravel()
                    itime += 1
            if ndofs:
                dofhist['a'][:,j] = tempa[np.ix_(dofs-1)].ravel()
                dofhist['da'][:,j] = tempda[np.ix_(dofs-1)].ravel()

    return modelhist, dofhist


def step2(K,C,M,f,a0,da0,bc,ip,times,dofs):
    """
    Algorithm for dynamic solution of second-order FE equations considering boundary conditions.

    Parameters
    ----------
    K : array_like
        Global stiffness matrix, shape (ndof, ndof).
    C : array_like
        Global damping matrix, shape (ndof, ndof).
        If there is no damping in the system, simply set C=[].
    M : array_like
        Global mass matrix, shape (ndof, ndof).
    f : array_like
        Global load vector, shape (ndof, nstep + 1).
        If shape (ndof, 1), the values are kept constant during time integration.
    a0 : array_like
        Initial displacement vector a(0), shape (ndof, 1).
    da0 : array_like
        Initial velocity vector v(0), shape (ndof, 1).
    bc : array_like
        Boundary condition matrix, shape (nbc, nstep + 2).
        where nbc = number of prescribed degrees of freedom (either constant or time-dependent).
        The first column contains the numbers of the prescribed degrees of freedom
        and the subsequent columns contain the time history.
        If shape (nbc, 2), the values from the second column are kept constant
        during time integration.
    ip : array_like
        Array [dt, tottime, alpha, delta], where
        dt is the size of the time increment,
        tottime is the total time,
        alpha and delta are time integration constants for the Newmark family of methods.
        Frequently used values of alpha and delta are:
        alpha=1/4, delta=1/2: average acceleration (trapezoidal) rule,
        alpha=1/6, delta=1/2: linear acceleration,
        alpha=0, delta=1/2: central difference.
    times : array_like
        Array [t(i) ...] of times at which output should be written to a, da and d2a.
    dofs : array_like
        Array [dof(i) ...] of degree of freedom numbers for which history output
        should be written to ahist, dahist and d2ahist.

    Returns
    -------
    modelhist : dict
        Dictionary containing solution history for the whole model with keys:
        
        - 'a' : ndarray
            Displacement values at all timesteps, alternatively at times specified in 'times',
            shape (ndof, nstep + 1) or (ndof, ntimes).
        - 'da' : ndarray
            Velocity values at all timesteps, alternatively at times specified in 'times',
            shape (ndof, nstep + 1) or (ndof, ntimes).
        - 'd2a' : ndarray
            Acceleration values at all timesteps, alternatively at times specified in 'times',
            shape (ndof, nstep + 1) or (ndof, ntimes).
    dofhist : dict
        Dictionary containing solution history for the degrees of freedom selected in 'dofs' with keys:
        
        - 'a' : ndarray
            Displacement time history at the dofs specified in 'dofs',
            shape (ndof, nstep + 1).
        - 'da' : ndarray
            Velocity time history at the dofs specified in 'dofs',
            shape (ndof, nstep + 1).
        - 'd2a' : ndarray
            Acceleration time history at the dofs specified in 'dofs',
            shape (ndof, nstep + 1).
    """
    ndof, _ = K.shape
    if not np.array(C).any():
        C = np.zeros((ndof,ndof))
    dt, tottime, alpha, delta = ip
    b1 = dt*dt*0.5*(1-2*alpha)
    b2 = (1-delta)*dt
    b3 = delta*dt
    b4 = alpha*dt*dt

    nstep = 1
    if np.array(f).any():
        _, ncf = f.shape
        if ncf>1:
            nstep = ncf-1

    if np.array(bc).any():
        _, ncb = bc.shape
        if ncb>2:
            nstep = ncb-2
        bound = 1
    if not np.array(bc).any():
        bound = 0

    ns = int(tottime/dt)
    if (ns < nstep or nstep==1):
        nstep=ns

    tf = np.zeros((ndof,nstep+1))
    if np.array(f).any():
        if ncf==1:
            tf = f[:,0].reshape(-1,1)@np.ones((1,nstep+1))
        if ncf>1:
            tf = np.copy(f)

    modelhist = {}
    sa=0
    if not np.array(times).any():
        ntimes=0
        sa=1
        modelhist['a'] = np.zeros((ndof,nstep+1))
        modelhist['da'] = np.zeros((ndof,nstep+1))
        modelhist['d2a'] = np.zeros((ndof,nstep+1))
    else:
        ntimes = len(times)
        if ntimes:
            sa=2
            modelhist['a'] = np.zeros((ndof,ntimes))
            modelhist['da'] = np.zeros((ndof,ntimes))
            modelhist['d2a'] = np.zeros((ndof,ntimes))

    dofhist = {}
    if np.array(dofs).all():
        ndofs = len(dofs)
        if ndofs:
            dofhist['a'] = np.zeros((ndofs,nstep+1))
            dofhist['da'] = np.zeros((ndofs,nstep+1))
            dofhist['d2a'] = np.zeros((ndofs,nstep+1))
    else:
        ndofs=0

    itime = 0

    # Calculate initial second time derivative d2a0
    d2a0 = np.linalg.solve(M,tf[:,0].reshape(-1,1) - C@da0 - K@a0)
    # Save initial values
    if sa==1:
        modelhist['a'][:,0] = a0.ravel()
        modelhist['da'][:,0] = da0.ravel()
        modelhist['d2a'][:,0] = d2a0.ravel()
    elif sa==2:
        if times[itime]==0:
            modelhist['a'][:,itime] = a0.ravel()
            modelhist['da'][:,itime] = da0.ravel()
            modelhist['d2a'][:,itime] = d2a0.ravel()
            itime += 1

    if ndofs:
        dofhist['a'][:,0] = a0[np.ix_(dofs-1)].ravel()
        dofhist['da'][:,0] = da0[np.ix_(dofs-1)].ravel()
        dofhist['d2a'][:,0] = d2a0[np.ix_(dofs-1)].ravel()

    # Reduce matrices due to bcs
    tempa = np.zeros((ndof,1))
    tempda = np.zeros((ndof,1))
    tempd2a = np.zeros((ndof,1))
    fdof=np.arange(1,ndof+1).astype(int)
    if bound:
        nrb, ncb = bc.shape
        if ncb==2:
            pa = bc[:,1].reshape(-1,1)@np.ones((1,nstep+1))
            pda = np.zeros((nrb,nstep+1))
        elif ncb>2:
            pa = np.copy(bc[:,1:])
            pda1 = (pa[:,1]-pa[:,0])/dt
            pdarest = (pa[:,1:] - pa[:,0:-1])/dt
            pda = np.hstack((pda1.reshape(-1,1),pdarest))
        pdof = np.copy(bc[:,0]).astype(int)
        fdof = np.setdiff1d(fdof,pdof).astype(int) - 1
        pdof -= 1 #adjusting for indexing starting from 0
        Keff = M[np.ix_(fdof,fdof)] + b3*C[np.ix_(fdof,fdof)] +b4*K[np.ix_(fdof,fdof)]
    else:
        fdof -= 1 #adjusting for indexing starting from 0
        Keff = M + b3*C + b4*K

    L, U = lu(Keff,permute_l=True)
    anew = a0[np.ix_(fdof)]
    danew = da0[np.ix_(fdof)]
    d2anew = d2a0[np.ix_(fdof)]

    # Iterate over time steps
    for j in range(1,nstep+1):
        time = dt*j
        aold = np.copy(anew)
        daold = np.copy(danew)
        d2aold = np.copy(d2anew)
        apred = aold + dt*daold + b1*d2aold
        dapred = daold + b2*d2aold
        if not bound:
            reff = tf[:,j].reshape(-1,1) - C@dapred - K@apred
        else:
            pdeff = C[np.ix_(fdof,pdof)]@pda[:,j].reshape(-1,1) + K[np.ix_(fdof,pdof)]@pa[:,j].reshape(-1,1)
            reff = tf[np.ix_(fdof),j].reshape(-1,1) - C[np.ix_(fdof,fdof)]@dapred - K[np.ix_(fdof,fdof)]@apred - pdeff
        y = np.linalg.solve(L,reff)
        d2anew = np.linalg.solve(U,y)
        anew = apred + b4*d2anew
        danew = dapred + b3*d2anew
        # Save to modelhist and dofhist
        if bound:
            tempa[np.ix_(pdof)] = pa[:,j].reshape(-1,1)
            tempda[np.ix_(pdof)] = pda[:,j].reshape(-1,1)
        tempa[np.ix_(fdof)] = anew
        tempda[np.ix_(fdof)] = danew
        tempd2a[np.ix_(fdof)] = d2anew
        if sa==1:
            modelhist['a'][:,j] = tempa.ravel()
            modelhist['da'][:,j] = tempda.ravel()
            modelhist['d2a'][:,j] = tempd2a.ravel()
        elif sa==2:
            if ntimes and itime < ntimes:
                if time >= times[itime]:
                    modelhist['a'][:,itime] = tempa.ravel()
                    modelhist['da'][:,itime] = tempda.ravel()
                    modelhist['d2a'][:,itime] = tempd2a.ravel()
                    itime += 1
            if ndofs:
                dofhist['a'][:,j] = tempa[np.ix_(dofs-1)].ravel()
                dofhist['da'][:,j] = tempda[np.ix_(dofs-1)].ravel()
                dofhist['d2a'][:,j] = tempd2a[np.ix_(dofs-1)].ravel()

    return modelhist, dofhist


def extract_eldisp(edof, a):
    """
    Extract element displacements from the global displacement vector according to the topology matrix edof.
    
    Parameters
    ----------
    a : array_like
        The global displacement vector.
    edof : array_like
        DOF topology array.
    
    Returns
    -------
    ed : ndarray
        Element displacement array.
    """

    ed = None

    if edof.ndim == 1:
        nDofs = len(edof)
        ed = np.zeros([nDofs])
        idx = edof-1
        ed[:] = a[np.ix_(idx)].T
    else:
        nElements = edof.shape[0]
        nDofs = edof.shape[1]
        ed = np.zeros([nElements, nDofs])
        i = 0
        for row in edof:
            idx = row-1
            ed[i, :] = a[np.ix_(idx)].T
            i += 1

    return ed


extractEldisp = extract_eldisp
extract_ed = extract_eldisp


def statcon(K, f, cd):
    """
    Condensation of static FE-equations according to the vector cd.

    Parameters
    ----------
    K : array_like
        Global stiffness matrix, shape (nd, nd).
    f : array_like
        Global load vector, shape (nd, 1).
    cd : array_like
        Vector containing dof's to be eliminated, shape (nc, 1), where nc is number of condensed dof's.

    Returns
    -------
    K1 : ndarray
        Condensed stiffness matrix, shape (nd-nc, nd-nc).
    f1 : ndarray
        Condensed load vector, shape (nd-nc, 1).
    """
    nd = K.shape[0]
    
    # Ensure cd is properly shaped (flatten to 1D array)
    cd = np.asarray(cd).flatten()
    # Adjust for 0-based indexing
    cd = (cd-1).astype(int)

    # Create indices for free (a) and constrained (b) dofs
    aindx = np.setdiff1d(np.arange(nd), cd)
    bindx = cd

    # Extract submatrices
    Kaa = K[np.ix_(aindx, aindx)]
    Kab = K[np.ix_(aindx, bindx)]
    Kbb = K[np.ix_(bindx, bindx)]

    # Extract portions of force vector
    f = np.asarray(f)
    if f.ndim > 1 and f.shape[1] > 1:
        # Handle matrix-like f
        f = f.reshape(-1, 1)
    
    fa = f[aindx]
    fb = f[bindx]

    # Calculate the condensed matrices
    Kbb_inv = np.linalg.inv(Kbb)
    K1 = Kaa - Kab @ Kbb_inv @ Kab.T
    f1 = fa - Kab @ Kbb_inv @ fb

    return K1, f1

def c_mul(a, b):
    """
    Multiply two integers with 32-bit overflow handling.

    Parameters
    ----------
    a : int
        First integer.
    b : int
        Second integer.

    Returns
    -------
    int
        Product of a and b with 32-bit overflow handling.
    """
    return eval(hex((int(a) * b) & 0xFFFFFFFF)[:-1])


def dof_hash(dof):
    """
    Compute a hash value for a degree of freedom array.

    Parameters
    ----------
    dof : array_like
        Degree of freedom array.

    Returns
    -------
    int
        Hash value for the DOF array.
    """
    if len(dof) == 1:
        return dof[0]
    value = 0x345678
    for item in dof:
        value = c_mul(1000003, value) ^ hash(item)
    value = value ^ len(dof)
    if value == -1:
        value = -2
    return value

dofHash = dof_hash

def create_dofs(nCoords: int, nDof: int) -> NDArray[np.integer]:
    """
    Create degree of freedom (DOF) array.

    Parameters
    ----------
    nCoords : int
        Number of coordinates (nodes).
    nDof : int
        Number of degrees of freedom per coordinate.

    Returns
    -------
    ndarray
        DOF array, shape (nCoords, nDof), with sequential DOF numbering starting from 1.
    """
    return np.arange(nCoords*nDof).reshape(nCoords, nDof)+1


createdofs = create_dofs


def coordxtr(edof, coords, dofs, nen=-1):
    """
    Create element coordinate matrices ex, ey, ez from edof coord and dofs matrices.

    Parameters
    ----------
    edof : array_like
        Element topology array, shape (nel, nen * nnd), where nel is number of elements,
        nen is number of element nodes, and nnd is number of node DOFs.
    coords : array_like
        Node coordinates array, shape (ncoords, ndims), where ncoords is number of coordinates
        and ndims is node dimensions.
    dofs : array_like
        DOF array, shape (ncoords, nnd), where nnd is number of node DOFs.
    nen : int, optional
        Number of element nodes. If -1, calculated from edof and dofs.

    Returns
    -------
    ex : ndarray
        Element x-coordinates, returned if ndims = 1.
    ex, ey : tuple of ndarray
        Element x and y coordinates, returned if ndims = 2.
    ex, ey, ez : tuple of ndarray
        Element x, y and z coordinates, returned if ndims = 3.
    """

    # Create dictionary with dof indices

    dofDict = {}
    nDofs = np.size(dofs, 1)
    nElements = np.size(edof, 0)
    n_element_dofs = np.size(edof, 1)
    nDimensions = np.size(coords, 1)
    nElementDofs = np.size(edof, 1)

    if nen == -1:
        nElementNodes = int(nElementDofs/nDofs)
    else:
        nElementNodes = nen

    if nElementNodes*nDofs != n_element_dofs:
        nDofs = nElementNodes*nDofs - n_element_dofs
        user_warning(
            "dofs/edof mismatch. Using %d dofs per node when indexing." % nDofs)

    idx = 0
    for dof in dofs:
        #dofDict[dofHash(dof)] = idx
        dofDict[hash(tuple(dof[0:nDofs]))] = idx
        idx += 1

    # Loop over edof and extract element coords

    ex = np.zeros((nElements, nElementNodes))
    ey = np.zeros((nElements, nElementNodes))
    ez = np.zeros((nElements, nElementNodes))

    elementIdx = 0
    for etopo in edof:
        for i in range(nElementNodes):
            i0 = i*nDofs
            i1 = i*nDofs+nDofs-1
            dof = []
            if i0 == i1:
                dof = [etopo[i*nDofs]]
            else:
                dof = etopo[i*nDofs:(i*nDofs+nDofs)]

            nodeCoord = coords[dofDict[hash(tuple(dof[0:nDofs]))]]

            if nDimensions >= 1:
                ex[elementIdx, i] = nodeCoord[0]
            if nDimensions >= 2:
                ey[elementIdx, i] = nodeCoord[1]
            if nDimensions >= 3:
                ez[elementIdx, i] = nodeCoord[2]

        elementIdx += 1

    if nDimensions == 1:
        return ex

    if nDimensions == 2:
        return ex, ey

    if nDimensions == 3:
        return ex, ey, ez


coord_extract = coordxtr


def hooke(ptype, E, v):
    """
    Calculate the material matrix for a linear elastic and isotropic material.

    Parameters
    ----------
    ptype : int
        Analysis type:
        1 : plane stress
        2 : plane strain
        3 : axisymmetry
        4 : three dimensional
    E : float
        Young's modulus.
    v : float
        Poisson's ratio.

    Returns
    -------
    D : ndarray
        Material matrix.
    """

    if ptype == 1:
        D = E*np.matrix(
            [[1, v, 0],
             [v, 1, 0],
             [0, 0, (1-v)/2]]
        )/(1-v**2)
    elif ptype == 2:
        D = E/(1+v)*np.matrix(
            [[1-v, v, v, 0],
             [v, 1-v, v, 0],
             [v, v, 1-v, 0],
             [0, 0, 0, (1-2*v)/2]]
        )/(1-2*v)
    elif ptype == 3:
        D = E/(1+v)*np.matrix(
            [[1-v, v, v, 0],
             [v, 1-v, v, 0],
             [v, v, 1-v, 0],
             [0, 0, 0, (1-2*v)/2]]
        )/(1-2*v)
    elif ptype == 4:
        D = E*np.matrix(
            [[1-v, v, v, 0, 0, 0],
             [v, 1-v, v, 0, 0, 0],
             [v, v, 1-v, 0, 0, 0],
             [0, 0, 0, (1-2*v)/2, 0, 0],
             [0, 0, 0, 0, (1-2*v)/2, 0],
             [0, 0, 0, 0, 0, (1-2*v)/2]]
        )/(1+v)/(1-2*v)
    else:
        info("ptype not supported.")

    return D


def effmises(es, ptype):
    """
    Calculate effective von mises stresses.

    Parameters
    ----------
    es : array_like
        Element stress matrix [[sigx, sigy, [sigz], tauxy], [...]], one row for each element.
    ptype : int
        Analysis type:
        1 : plane stress
        2 : plane strain
        3 : axisymmetry
        4 : three dimensional

    Returns
    -------
    eseff : ndarray
        Effective stress array [eseff_0, ..., eseff_nel-1].
    """

    nel = np.size(es, 0)
    escomps = np.size(es, 1)

    eseff = np.zeros([nel])

    if ptype == 1:
        sigxx = es[:, 0]
        sigyy = es[:, 1]
        sigxy = es[:, 2]
        eseff = np.sqrt(sigxx*sigxx+sigyy*sigyy-sigxx*sigyy+3*sigxy*sigxy)
        return eseff


def stress2nodal(eseff, edof):
    """
    Convert element effective stresses to nodal effective
    stresses.
    
    Parameters:
        
        eseff  = [eseff_0 .. eseff_nel-1] 
        edof   = [dof topology array]
    
    Returns:
    
        ev:     element value array [[ev_0_0 ev_0_1 ev_0_nen-1 ]
                                      ..
                                      ev_nel-1_0 ev_nel-1_1 ev_nel-1_nen-1]
                      
    """

    values = np.zeros(edof.max())
    elnodes = int(np.size(edof, 1) / 2)

    for etopo, eleseff in zip(edof, eseff):
        values[etopo-1] = values[etopo-1] + eleseff / elnodes

    evtemp = extractEldisp(edof, values)
    ev = evtemp[:, range(0, elnodes*2, 2)]

    return ev


def beam2crd_old(ex, ey, ed, mag):
    """
    Calculate the element continuous displacements for multiple identical 2D Bernoulli beam elements.

    Parameters
    ----------
    ex : array_like
        Element node x-coordinates.
    ey : array_like
        Element node y-coordinates.
    ed : array_like
        Element displacement matrix.
    mag : float
        Magnification factor.

    Returns
    -------
    excd : ndarray
        Continuous x-coordinates.
    eycd : ndarray
        Continuous y-coordinates.

    History
    -------
    LAST MODIFIED: P-E AUSTRELL 1993-10-15
                   J Lindemann 2021-12-30 (Python)
    """
    nie, ned = ed.shape

    excd = np.zeros([nie, 20])
    eycd = np.zeros([nie, 20])

    for i in range(nie):

        b = np.array([ex[i, 1]-ex[i, 0], ey[i, 1]-ey[i, 0]])
        L = np.sqrt(b@np.transpose(b)).item()
        n = b/L

        G = np.array([
            [n[0], n[1],  0,    0,    0,   0],
            [-n[1], n[0],  0,    0,    0,   0],
            [0,    0,     1,    0,    0,   0],
            [0,    0,     0,   n[0], n[1], 0],
            [0,    0,     0,  -n[1], n[0], 0],
            [0,    0,     0,    0,    0,   1]
        ])

        d = ed[i, :]
        dl = G @ d

        xl = np.linspace(0.0, L, 20)
        one = np.ones(xl.shape)

        Cis = np.array([
            [-1.0, 1.0],
            [L, 0.0]
        ]) / L

        ds = np.array([dl[0], dl[3]]).reshape(2, 1)

        xl_one = np.transpose(np.vstack((xl, one)))

        ul = np.transpose(xl_one@Cis@ds)  # [20x1][2]

        Cib = np.array([
            [12,      6*L,  -12,    6*L],
            [-6*L, -4*L**2,  6*L, -2*L**2],
            [0,      L**3,    0,     0],
            [L**3,      0,     0,     0]
        ])/L**3

        db = np.array([dl[1], dl[2], dl[4], dl[5]]).reshape(4, 1)
        vl = np.transpose(np.transpose(
            np.vstack((xl**3/6, xl**2/2, xl, one)))@Cib@db)

        cld = np.vstack((ul, vl))
        A = np.array([
            [n[0], -n[1]],
            [n[1], n[0]]
        ])
        cd = A@cld

        # [2,1] x [1,20] + [2 x 1] x [1 x 20]
        #    [2 x 20]    +      [2 x 20]

        AA = A[:, 0].reshape(2, 1)
        XL = xl.reshape(1, 20)
        xyc = AA@XL + np.array([[ex[i, 0]], [ey[i, 0]]])@one.reshape(1, 20)

        excd[i, :] = xyc[0, :]+mag*cd[0, :]
        eycd[i, :] = xyc[1, :]+mag*cd[1, :]

    return excd, eycd


def beam2crd(ex=None, ey=None, ed=None, mag=None):
    """
    Calculate the element continuous displacements for multiple identical 2D Bernoulli beam elements.

    Parameters
    ----------
    ex : array_like, optional
        Element node x-coordinates.
    ey : array_like, optional
        Element node y-coordinates.
    ed : array_like, optional
        Element displacement matrix.
    mag : float, optional
        Magnification factor.

    Returns
    -------
    excd : ndarray
        Continuous x-coordinates.
    eycd : ndarray
        Continuous y-coordinates.

    History
    -------
    LAST MODIFIED: P-E AUSTRELL 1993-10-15
    """

    nie, ned = ed.shape

    n_coords = 21

    excd = np.zeros([nie, n_coords])
    eycd = np.zeros([nie, n_coords])

    for i in range(nie):
        b = np.array([ex[i, 1] - ex[i, 0], ey[i, 1] - ey[i, 0]])
        L = np.sqrt(b @ np.transpose(b))
        n = b / L

        G = np.array([
            [n[0], n[1], 0, 0, 0, 0],
            [-n[1], n[0], 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, n[0], n[1], 0],
            [0, 0, 0, - n[1], n[0], 0],
            [0, 0, 0, 0, 0, 1]
        ])

        d = np.transpose(ed[i, :])
        dl = G @ d
        xl = np.transpose(np.linspace(0, L, n_coords))
        one = np.ones(xl.shape)

        Cis = np.array([
            [-1, 1],
            [L, 0]
        ]) / L

        ds = np.array([dl[0], dl[3]]).reshape(2, 1)
        xl_one = np.transpose(np.vstack((xl, one)))
        ul = np.transpose(xl_one@Cis@ds)  # [20x1][2]

        Cib = np.array([
            [12, 6 * L, - 12, 6 * L],
            [- 6 * L, - 4 * L ** 2, 6 * L, - 2 * L ** 2],
            [0, L ** 3, 0, 0],
            [L ** 3, 0, 0, 0]
        ]) / L ** 3

        db = np.array([dl[1], dl[2], dl[4], dl[5]]).reshape(4, 1)
        vl = np.transpose(np.transpose(
            np.vstack((xl**3/6, xl**2/2, xl, one)))@Cib@db)

        cld = np.vstack((ul, vl))
        A = np.array([
            [n[0], -n[1]],
            [n[1], n[0]]
        ])
        cd = A@cld

        # [2,1] x [1,20] + [2 x 1] x [1 x 20]
        #    [2 x 20]    +      [2 x 20]

        AA = A[:, 0].reshape(2, 1)
        XL = xl.reshape(1, n_coords)
        xyc = AA@XL + np.array([[ex[i, 0]], [ey[i, 0]]]
                               )@one.reshape(1, n_coords)

        excd[i, :] = xyc[0, :]+mag*cd[0, :]
        eycd[i, :] = xyc[1, :]+mag*cd[1, :]

    return excd, eycd
