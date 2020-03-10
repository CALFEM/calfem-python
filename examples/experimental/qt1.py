# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 16:23:37 2017

@author: Jonas Lindemann
"""

import sys
import time

from calfem.qt5 import *

import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis as cfv

class MainWindow(QMainWindow):
    def __init__(self):
        """Constructor"""
        super(MainWindow, self).__init__()       
        
if __name__ == "__main__":
    
    app = QApplication(sys.argv)
    widget = MainWindow()
    widget.show()
    sys.exit(app.exec_())  
