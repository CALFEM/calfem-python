#!/bin/env python
# -*- coding: iso-8859-15 -*-

from OpenGL.GL import *
from OpenGL.GLUT import *
from PyQt import QtGui, QtCore, QtOpenGL

import colorsys

def drawBitmapText(text, font=GLUT_BITMAP_TIMES_ROMAN_24):
    for c in text:
        glutBitmapCharacter(font, ord(c))

def floatRgb(mag, cmin, cmax):
    """
    Return a tuple of floats between 0 and 1 for the red, green and
    blue amplitudes.
    """
    
    try:
           # normalize to [0,1]
           x = float(mag-cmin)/float(cmax-cmin)
    except:
           # cmax = cmin
           x = 0.5
           
    red, green, blue = colorsys.hsv_to_rgb(x*240./360, 1.0, 1.0)
    return (red, green, blue)

class OpenGLFrame(QtOpenGL.QGLWidget):
    def __init__(self, parent = None):
        super(OpenGLFrame, self).__init__(parent)

    def paintGL(self):
        self.onDraw()

    def resizeGL(self, w, h):
        self.onReshape(w, h)

    def initializeGL(self):
        self.onInitGL()
        
    def Show(self):
        self.show()
        
    def SwapBuffers(self):
        pass
        
    #
    # GLFrame OpenGL Event Handlers

    def onInitGL(self):
        """Initialize OpenGL for use in the window."""
        glClearColor(1, 1, 1, 1)
        glDisable(GL_DEPTH_TEST)

    def onReshape(self, width, height):
        """Reshape the OpenGL viewport based on the dimensions of the window."""
        glViewport(0, 0, width, height)

        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glOrtho(-0.5, 0.5, -0.5, 0.5, -1, 1)

        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()

    def onDraw(self):
        "Draw the window."
        glClear(GL_COLOR_BUFFER_BIT)

        # Drawing an example triangle in the middle of the screen
        glBegin(GL_TRIANGLES)
        glColor(1, 0, 0)
        glVertex(-.25, -.25)
        glColor(0, 1, 0)
        glVertex(.25, -.25)
        glColor(0, 0, 1)
        glVertex(0, .25)
        glEnd()
        
        
class ElementView(OpenGLFrame):
    def __init__(self, parent = None, a = 0, b = 0):
        super(ElementView, self).__init__(parent)
        self._initView()
        
    def _initView(self):
        
        self._limits = [1e300, 1e300, -1e300, -1e300]
        self._nodeLimits = [1e300, -1e300]
        self._maxNodeValue = -1e300
        self._minNodeValue = 1e300
        self._maxElementValue = -1e300
        self._minElementValue = 1e300
        self._ex = None
        self._ey = None
        self._ed = None
        self._ev = None
        self._magnfac = 0.1
        self._elementNodes = 3
        self._dofsPerNode = 1
        self._showDimension = 1
        self._width = 100
        self._height = 100
        self._showMesh = True
        self._showNodalValues = True
        self._showDisplacements = False
        self._showElementValues = False
        
        self.drawAnnotations = None
        self.drawCustom = None
        
    def calcLimits(self):
        
        self._limits = [1e300, 1e300, -1e300, -1e300]
        
        if not (self._ex is None) and not (self._ey is None):
            for elCoords in self._ex:
                xmin = min(elCoords)
                xmax = max(elCoords)
                
                if xmin<self._limits[0]:
                    self._limits[0]=xmin
                if xmax>self._limits[2]:
                    self._limits[2]=xmax

            for elCoords in self._ey:
                ymin = min(elCoords)
                ymax = max(elCoords)
                
                if ymin<self._limits[1]:
                    self._limits[1]=ymin
                if ymax>self._limits[3]:
                    self._limits[3]=ymax
                    
    def calcNodeLimits(self):
        if self.dofsPerNode == 1:
            self._maxNodeValue = self._ed.max()
            self._minNodeValue = self._ed.min()
        else:
            self._maxNodeValue = abs(self._ed).max()
            
    def calcElementLimits(self):
        self._maxElementValue = self._ev.max()
        self._minElementValue = self._ev.min()
    
    def calcScaling(self):

        self._k = 0.8
        
        factor1 = self._k*self._width/(self._limits[2]-self._limits[0])
        factor2 = self._k*self._height/(self._limits[3]-self._limits[1])
        
        self._scaleFactor = 1.0
        
        if factor1<factor2:
            self._scaleFactor = factor1
        else:
            self._scaleFactor = factor2
            
        self._x0 = -self._scaleFactor*self._limits[0] + (1-self._k)*self._width/2.0
        self._y0 = self._height - self._scaleFactor*self._limits[1] - (1-self._k)*self._height/2.0

    def worldToScreen(self, x, y):
        return (int(self._x0 + self._scaleFactor*x), int(self._y0 - self._scaleFactor*y))
        
    def drawMesh(self):
             
        # Draw elements
                              
        glBegin(GL_LINES)
              
        for elx, ely in zip(self._ex, self._ey):
            
            if self._elementNodes == 3:
            
                (sx1, sy1) = self.worldToScreen(elx[0], ely[0])
                (sx2, sy2) = self.worldToScreen(elx[1], ely[1])
                (sx3, sy3) = self.worldToScreen(elx[2], ely[2])
                
                glColor(0.5, 0.5, 0.5)
                glVertex(sx1,sy1)
                glVertex(sx2,sy2)
                glVertex(sx2,sy2)
                glVertex(sx3,sy3)
                glVertex(sx3,sy3)
                glVertex(sx1,sy1)
                
            elif self._elementNodes == 2:
            
                (sx1, sy1) = self.worldToScreen(elx[0], ely[0])
                (sx2, sy2) = self.worldToScreen(elx[1], ely[1])

                glColor(0.5, 0.5, 0.5)
                glVertex(sx1,sy1)
                glVertex(sx2,sy2)
            
        glEnd()
        
    def drawNodalValues(self):
                
        # Draw nodal values
        
        if self._elementNodes == 3:
            glBegin(GL_TRIANGLES)
        else:
            return
        
        for elx, ely, eld in zip(self._ex, self._ey, self._ed):
            
            if self._elementNodes == 3:
            
                (sx1, sy1) = self.worldToScreen(elx[0], ely[0])
                (sx2, sy2) = self.worldToScreen(elx[1], ely[1])
                (sx3, sy3) = self.worldToScreen(elx[2], ely[2])
                           
                if self._dofsPerNode == 1:
                    c1 = floatRgb(eld[0], self._maxNodeValue, self._minNodeValue)
                    c2 = floatRgb(eld[1], self._maxNodeValue, self._minNodeValue)
                    c3 = floatRgb(eld[2], self._maxNodeValue, self._minNodeValue)
                else:
                    c1 = floatRgb(eld[0], self._maxNodeValue, self._minNodeValue)
                    c2 = floatRgb(eld[1*self._dofsPerNode+(self._showDimension-1)], self._maxNodeValue, self._minNodeValue)
                    c3 = floatRgb(eld[2*self._dofsPerNode+(self._showDimension-1)], self._maxNodeValue, self._minNodeValue)
                        
                glColor3f(c1[0], c1[1], c1[2])
                glVertex(sx1,sy1)
                glColor3f(c2[0], c2[1], c2[2])
                glVertex(sx2,sy2)
                glColor3f(c3[0], c3[1], c3[2])
                glVertex(sx3,sy3)
            
        glEnd()
        
    def drawElementValues(self):
                
        # Draw element values
        
        if self._elementNodes == 3:
            glBegin(GL_TRIANGLES)
        else:
            return
        
        for elx, ely, elv in zip(self._ex, self._ey, self._ev):
            
            if self._elementNodes == 3:
            
                (sx1, sy1) = self.worldToScreen(elx[0], ely[0])
                (sx2, sy2) = self.worldToScreen(elx[1], ely[1])
                (sx3, sy3) = self.worldToScreen(elx[2], ely[2])
                           
                if self._dofsPerNode == 1:
                    c1 = floatRgb(elv, self._maxElementValue, self._minElementValue)
                else:
                    c1 = floatRgb(elv, self._maxElementValue, self._minElementValue)
                        
                glColor3f(c1[0], c1[1], c1[2])
                glVertex(sx1,sy1)
                glVertex(sx2,sy2)
                glVertex(sx3,sy3)
            
        glEnd()
        
    def drawDisplacements(self):
        
        # Draw elements
        
        scl = self._magnfac*self.modelWidth/self._maxNodeValue
        
        glBegin(GL_LINES)
        
        for elx, ely, eld in zip(self._ex, self._ey, self._ed):
            
            if self._elementNodes == 3:
            
                (sx1, sy1) = self.worldToScreen(elx[0]+scl*eld[0], ely[0]+scl*eld[1])
                (sx2, sy2) = self.worldToScreen(elx[1]+scl*eld[2], ely[1]+scl*eld[3])
                (sx3, sy3) = self.worldToScreen(elx[2]+scl*eld[4], ely[2]+scl*eld[5])
                                           
                glColor(0.3, 0.3, 0.3)
                glVertex(sx1,sy1)
                glVertex(sx2,sy2)
                glVertex(sx2,sy2)
                glVertex(sx3,sy3)
                glVertex(sx3,sy3)
                glVertex(sx1,sy1)

            elif self._elementNodes == 2:
            
                (sx1, sy1) = self.worldToScreen(elx[0]+scl*eld[0], ely[0]+scl*eld[1])
                (sx2, sy2) = self.worldToScreen(elx[1]+scl*eld[2], ely[1]+scl*eld[3])

                glColor(0.3, 0.3, 0.3)
                glVertex(sx1,sy1)
                glVertex(sx2,sy2)

        glEnd()

            
    def onReshape(self, width, height):
        """
        Reshape the OpenGL viewport based on the dimensions of the window.
        """
        glViewport(0, 0, width, height)
        
        self._width = width
        self._height = height
        
        self.calcScaling()

        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glOrtho(0.0, width, height, 0.0, -1, 1)

        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        
        glEnable(GL_POINT_SMOOTH)
        glEnable(GL_LINE_SMOOTH)
        glEnable(GL_POLYGON_SMOOTH)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA)
        glHint(GL_LINE_SMOOTH_HINT,GL_NICEST)
        glLineWidth(1.0)
        glDisable(GL_DEPTH_TEST)
        glDisable(GL_LIGHTING)

    def onDraw(self, *args, **kwargs):
        """
        Draw the window.
        """
              
        glClear(GL_COLOR_BUFFER_BIT)
        
        if self._showElementValues:
            self.drawElementValues()
        if self._showNodalValues:
            self.drawNodalValues()
        if self._showMesh:
            self.drawMesh()
        if self._showDisplacements:
            self.drawDisplacements()
            
        if self.drawCustom!=None:
            self.drawCustom(self, self._width, self._height)
            
        if self.drawAnnotations!=None:
            self.drawAnnotations(self, self._width, self._height)
            
        self.SwapBuffers()
        
    def getModelHeight(self):
        return self._ey.max() - self._ey.min()
        
    def getModelWidth(self):
        return self._ex.max() - self._ex.min()

    def setEx(self, ex):
        self._ex = ex
        self._elementNodes = self._ex.shape[1]
        self.calcLimits()
        self.calcScaling()
        
    def getEx(self):
        return self._ex
    
    def setEy(self, ey):
        self._ey = ey
        self._elementNodes = self._ey.shape[1]
        self.calcLimits()
        self.calcScaling()
        
    def getEy(self):
        return self._ey
    
    def setEd(self, ed):
        self._ed = ed
        self.calcNodeLimits()
        
    def getEd(self):
        return self._ed
    
    def setShowMesh(self, showMesh):
        self._showMesh = showMesh
        
    def getShowMesh(self):
        return self._showMesh
    
    def setShowNodalValues(self, showNodalValues):
        self._showNodalValues = showNodalValues
        
    def getShowNodalValues(self):
        return self._showNodalValues
    
    def setShowElementValues(self, showElementValues):
        self._showElementValues = showElementValues
        
    def getShowElementValues(self):
        return self._showElementValues

    def setShowDisplacements(self, showDisplacements):
        self._showDisplacements = showDisplacements
        
    def getShowDisplacements(self):
        return self._showDisplacements
    
    def setDofsPerNode(self, dofsPerNode):
        self._dofsPerNode = dofsPerNode
        
    def getDofsPerNode(self):
        return self._dofsPerNode
    
    def setElementNodes(self, elementNodes):
        self._elementNodes = elementNodes
        
    def getElementNodes(self):
        return self._elementNodes
    
    def setMagnFac(self, magnfac):
        self._magnfac = magnfac
        
    def getMagnFac(self):
        return self._magnfac
    
    def getEv(self):
        return self._ev
    
    def setEv(self, value):
        self._ev = value
        self.calcElementLimits()
    
    ex = property(getEx, setEx)
    ey = property(getEy, setEy)
    ed = property(getEd, setEd)
    ev = property(getEv, setEv)
    dofsPerNode = property(getDofsPerNode, setDofsPerNode)
    elementNodes = property(getElementNodes, setElementNodes)
    modelWidth = property(getModelWidth)
    modelHeight = property(getModelHeight)
    magnfac = property(getMagnFac, setMagnFac)
    showMesh = property(getShowMesh, setShowMesh)
    showNodalValues = property(getShowNodalValues, setShowNodalValues)
    showDisplacements = property(getShowDisplacements, setShowDisplacements)
    showElementValues = property(getShowElementValues, setShowElementValues)        
        

