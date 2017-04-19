# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 09:44:29 2016

@author: lindemann
"""

import sys

from PyQt5.QtCore import pyqtSlot, pyqtSignal, QMetaObject
from PyQt5.QtWidgets import *
from PyQt5.QtGui import QPixmap
from PyQt5.uic import loadUi

class MainWindow(QMainWindow):
    """Main window class of our UI"""
    def __init__(self):
        """Constructor"""
        super(MainWindow, self).__init__()

        self.resize(640,480)
        self.move(50,50)
        self.setWindowTitle("MyWindow")

        self.executeButton = QPushButton("Tryck", self)
        self.executeButton.move(50,50)
        self.executeButton.resize(100,50)
        #self.button.clicked.connect(self.on_button_clicked)
        
        QMetaObject.connectSlotsByName(self)

    @pyqtSlot()
    def on_executeButton_clicked(self):
        print("Button pressed")


#    def show(self):
#        """Show and raise window"""
#        self.ui.show()
#        self.ui.raise_()

if __name__ == '__main__':

    app = QApplication(sys.argv)
    widget = MainWindow()
    widget.show()
    sys.exit(app.exec_())