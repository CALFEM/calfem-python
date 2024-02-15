# -*- coding: utf-8 -*-
"""
CALFEM UI module

Routines for interfacing wirt user interface toolkits.
"""

import os, sys

print("------------------------------------")
print("CALFEM/Python ui module initialising")
print("------------------------------------")
print()

from qtpy.QtCore import Slot, Signal, QThread
from qtpy.QtWidgets import QApplication, QDialog, QWidget, QMainWindow
from qtpy.QtGui import QPixmap
from qtpy.uic import loadUi

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

def init_qt_app():
    app = QApplication.instance()
    
    if app is None:
        print("No QApplication instance found. Creating one.")
        # if it does not exist then a QApplication is created
        app = QApplication(sys.argv)    
    else:
        print("QApplication instance found. Reusing.")

    return app
            
def loadUiWidget(uifilename, parent=None):
    """Load user interface file and return object model"""
    ui = loadUi(uifilename, parent)
    return ui

load_ui_widget = loadUiWidget
        
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

app_instance = appInstance