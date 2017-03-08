#!/bin/env python
# -*- coding: iso-8859-15 -*-

import os, sys

import numpy as np
import calfem.core as cfc
import logging as cflog

haveMatplotLib = True
haveMlab = True

def error(msg):
    cflog.error(" "+msg)

def info(msg):
    cflog.info(" "+msg)

#haveWx = True
#haveQt = True
    
#globalWxApp = None
#globalQtApp = None
#globalWindows = []
#
#try:
#    from PyQt import QtGui
#    from PyQt.QtOpenGL import *
#    from calfem.classes_qt4 import ElementView
#    globalQtApp = QtGui.QApplication(["PyCalfem"])
#except:
#    haveQt = False    
#      
#if not haveQt:
#    try:
#        import wx
#        from calfem.classes_wx import ElementView
#        globalWxApp = wx.App(0)
#    except: 
#        haveWx = False

class ElementProperties(object):
    def __init__(self):
        self.ep = {}
        self.attributes = {}
        
    def add(self, markerId, ep=[]):
        if not markerId in self.ep:
            self.ep[markerId] = ep
            
    def addAttribute(self, markerId, name, value):
        if not markerId in self.attributes:
            self.attributes[markerId] = {}
            self.attributes[markerId][name] = value

def enableLogging():
    cflog.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s', level=cflog.DEBUG)
    
def disableLogging():
    cflog.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s', level=cflog.NOTSET)
    


def readInt(f):
    """
    Read a row from file, f, and return a list of integers.
    """
    return list(map(int, f.readline().split()))
    
def readFloat(f):
    """
    Read a row from file, f, and return a list of floats.
    """
    return list(map(float, f.readline().split()))
    
def readSingleInt(f):
    """
    Read a single integer from a row in file f. All other values on row are discarded.
    """
    return readInt(f)[0] 

def readSingleFloat(f):
    """
    Read a single float from a row in file f. All other values on row are discarded.
    """
    return readFloat(f)[0]
    
def writeSingleFloat(f, floatValue):
    f.write("%g\n" % floatValue)
    
def writeSingleInt(f, intValue):
    f.write("%d\n" % intValue)

def writeFloatList(f, floatList):
    for floatValue in floatList:
        f.write("%g " % floatValue)
    f.write("\n")
    
def writeIntList(f, intList):
    for intValue in intList:
        f.write("%d " % intValue)
    f.write("\n")
    
def which(filename):
    """
    Return complete path to executable given by filename.
    """
    if not ('PATH' in os.environ) or os.environ['PATH'] == '':
        p = os.defpath
    else:
        p = os.environ['PATH']
                
    pathlist = p.split (os.pathsep)
    pathlist.append(".")
    pathlist.append("/bin")
    pathlist.append("/usr/bin")
    pathlist.append("/usr/local/bin")
    pathlist.append("/opt/local/bin")
    
    for path in pathlist:
        f = os.path.join(path, filename)
        if os.access(f, os.X_OK):
            return f
    return None

def applybc(boundaryDofs, bcPrescr, bcVal, marker, value=0.0, dimension=0):
    """
    Apply boundary condition to bcPresc and bcVal matrices.
    
    Parameters:
    
        boundaryDofs        Dictionary with boundary dofs.
        bcPresc             1-dim integer array containing prescribed dofs.
        bcVal               1-dim float array containing prescribed values.
        marker              Boundary marker to assign boundary condition.
        value               Value to assign boundary condition.
                            If not giben 0.0 is assigned.
        dimension           dimension to apply bc. 0 - all, 1 - x, 2 - y
                            
    """

    if marker in boundaryDofs:
        if (dimension==0):
            bcAdd = np.array(boundaryDofs[marker])
            bcAddVal = np.ones([np.size(bcAdd)])*value
        elif (dimension==1):
            bcAdd = boundaryDofs[marker][(dimension-1)::2]
            bcAddVal = np.ones([np.size(bcAdd)])*value
        else:
            bcAdd = boundaryDofs[marker][(dimension-1)::2]
            bcAddVal = np.ones([np.size(bcAdd)])*value

        return np.hstack([bcPrescr,bcAdd]), np.hstack([bcVal,bcAddVal])
    else:
        print("Error: Boundary marker", marker, "does not exist.")
        
def applybcnode(nodeIdx, dofs, bcPrescr, bcVal, value=0.0, dimension=0):
    
    if (dimension==0):
        bcAdd = np.asarray(dofs[nodeIdx])
        bcAddVal = np.ones([np.size(bcAdd)])*value
    elif (dimension==1):
        bcAdd = np.asarray(dofs[nodeIdx,dimension-1])
        bcAddVal = np.ones([np.size(bcAdd)])*value
    else:
        bcAdd = np.asarray(dofs[nodeIdx,dimension-1])
        bcAddVal = np.ones([np.size(bcAdd)])*value

    return np.hstack([bcPrescr,bcAdd]), np.hstack([bcVal,bcAddVal])
    

def applyforcenode(nodeIdx, dofs, f, value=0.0, dimension=0):

    if (dimension==0):
        f[dofs[nodeIdx]]+=value
    elif (dimension==1):
        f[dofs[nodeIdx,dimension-1]]+=value
    else:
        f[dofs[nodeIdx,dimension-1]]+=value
        
def applyforce(boundaryDofs, f, marker, value=0.0, dimension=0):
    """
    Apply boundary condition to bcPresc and bcVal matrices.
    
    Parameters:
    
        boundaryDofs        Dictionary with boundary dofs.
        f                   force matrix.
        marker              Boundary marker to assign boundary condition.
        value               Value to assign boundary condition.
                            If not giben 0.0 is assigned.
        dimension           dimension to apply force. 0 - all, 1 - x, 2 - y
                            
    """

    if marker in boundaryDofs:
        if dimension == 0:
            f[np.asarray(boundaryDofs[marker])-1] += value
        elif dimension == 1:
            f[np.asarray(boundaryDofs[marker][(dimension-1)::2])-1] += value
        elif dimension == 2:
            f[np.asarray(boundaryDofs[marker][(dimension-1)::2])-1] += value            
    else:
        print("Error: Boundary marker", marker, "does not exist.")

def applyforcetotal(boundaryDofs, f, marker, value=0.0, dimension=0):
    """
    Apply boundary force to f matrix. Total force, value, is
    distributed over all boundaryDofs defined by marker.
    
    Parameters:
    
        boundaryDofs        Dictionary with boundary dofs.
        f                   force matrix.
        marker              Boundary marker to assign boundary condition.
        value               Total force value to assign boundary condition.
                            If not giben 0.0 is assigned.
        dimension           dimension to apply force. 0 - all, 1 - x, 2 - y
                            
    """

    if marker in boundaryDofs:
        if dimension == 0:
            nDofs = len(boundaryDofs[marker])
            valuePerDof = value / nDofs
            f[np.asarray(boundaryDofs[marker])-1] += valuePerDof
        elif dimension == 1:
            nDofs = len(boundaryDofs[marker][(dimension-1)::2])
            valuePerDof = value / nDofs
            f[np.asarray(boundaryDofs[marker][(dimension-1)::2])-1] += valuePerDof
        elif dimension == 2:
            nDofs = len(boundaryDofs[marker][(dimension-1)::2])
            valuePerDof = value / nDofs
            f[np.asarray(boundaryDofs[marker][(dimension-1)::2])-1] += valuePerDof
    else:
        print("Error: Boundary marker", marker, "does not exist.")
               
def scalfact2(ex,ey,ed,rat=0.2):
    """
    Determine scale factor for drawing computational results, such as 
    displacements, section forces or flux.
    
    Parameters:
    
        ex, ey      element node coordinates
                       
        ed          element displacement matrix or section force matrix
    
        rat         relation between illustrated quantity and element size. 
                    If not specified, 0.2 is used.
        
    """

    nen = -1
    if ex.shape != ey.shape:
        print("ex and ey shapes do not match.")
        return 1.0
    
    dlmax = 0.
    edmax = 1.
    
    if np.rank(ex)==1:
        nen = ex.shape[0]
        nel = 1
        dxmax=ex.T.max()-ex.T.min()
        dymax=ey.T.max()-ey.T.min()
        dlmax=max(dxmax,dymax);
        edmax=abs(ed).max();
    else:
        nen = ex.shape[1]
        nel = ex.shape[0]
        dxmax=ex.T.max()-ex.T.min()
        dymax=ey.T.max()-ey.T.min()
        dlmax=max(dxmax,dymax);
        edmax=abs(ed).max();
        
    k = rat
    return k*dlmax/edmax
