# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 09:44:29 2016

@author: lindemann
"""

import os, sys

print("------------------------------------")
print("CALFEM/Python ui module initialising")
print("------------------------------------")
print()

from PyQt5.QtCore import pyqtSlot, pyqtSignal
from PyQt5.QtWidgets import QApplication, QDialog, QWidget, QMainWindow
from PyQt5.QtGui import QPixmap
from PyQt5.uic import loadUi

#from PyQt import QtGui, QtCore, uic


g_inSpyder = False

if any('SPYDER' in name for name in os.environ):
    print('Running in Spyder...')
    g_inSpyder = True
    
g_haveVisVis = True

try:
    import visvis as vv
except:
    g_haveVisVis = False
        
if g_haveVisVis:
    print("VisVis installed...")
else:
    print("VisVis not installed...")
    
def loadUiWidget(uifilename, parent=None):
    """Load user interface file and return object model"""
    ui = loadUi(uifilename, parent)
    return ui    
        
def appInstance(useVisVis=True):
    """Create a suitable application instance"""
    print("Creating application instance...")
    
    global g_haveVisVis
    
    app = None
    
    if g_haveVisVis and useVisVis:
        print("Using VisVis application instance...")
        app = vv.use()
        app.Create()
    else:
        print("Trying Qt application instance...") 
        app = QtGui.QApplication.instance()
        if app is None:
            print("Creating new Qt application instance...")
            app = QtGui.QApplication(sys.argv)
        else:
            print("Reusing existing Qt application instance...")
        
        if app!=None:
            app.Run = app.exec_
            
    if app is None:
        print("No application instance found. Exiting...")
        sys.exit(-1)

    return app    