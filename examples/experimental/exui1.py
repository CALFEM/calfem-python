# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 09:44:29 2016

@author: lindemann
"""

import sys

from qtpy.QtCore import Slot, Signal, QMetaObject
from qtpy.QtWidgets import *
from qtpy.QtGui import QPixmap
from qtpy.uic import loadUi

import calfem.vis_mpl as cfv

class MainWindow(QMainWindow):
    """Main window class of our UI"""
    def __init__(self):
        """Constructor"""
        super(MainWindow, self).__init__()

        self.resize(640,480)
        self.move(50,50)
        self.setWindowTitle("MyWindow")

        self.executeButton = QPushButton("Tryck", self)
        self.executeButton.setObjectName("executeButton")
        self.executeButton.move(50,50)
        self.executeButton.resize(100,50)
        
        QMetaObject.connectSlotsByName(self)

    @Slot()
    def on_executeButton_clicked(self):
        print("Button pressed")

def init_app():
    app = QApplication.instance()
    
    if app is None:
        print("No QApplication instance found. Creating one.")
        # if it does not exist then a QApplication is created
        app = QApplication(sys.argv)    
    else:
        print("QApplication instance found.")

    return app
        
def main_loop():
    
    app = QApplication.instance()
    
    if app is None:
        print("No QApplication instance found. Creating one.")
        # if it does not exist then a QApplication is created
        app = QApplication(sys.argv)    
    else:
        print("QApplication instance found.")
        
    # För matplotlib kompatibilitet
    
    # plt.show(block=False)
    
    app.exec_()      
        
        
def show_window():
    app = init_app()
    widget = MainWindow()
    widget.show()
    app.exec_()

if __name__ == '__main__':

    show_window()
    print("Hello")
    show_window()
    print("Hello after...")