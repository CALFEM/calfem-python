# -*- coding: utf-8 -*-

import sys

from calfem.qt5 import *

class GuiWindow(QWidget):
    def __init__(self, var_dict):
        """MyWindow constructor"""

        super().__init__()

        # Skapa gränssnittskontroller

        self.var_dict = var_dict

        self.init_gui()

    def init_gui(self):
        """Initiera gränssnitt"""

        self.setGeometry(300, 300, 600, 600)
        self.setWindowTitle("Parameter window")

        self.grid = QGridLayout(self)

        row = 0

        for key, value in self.var_dict.items():
            var_label = QLabel(key)
            var_edit = QLineEdit(str(value))
            self.grid.addWidget(var_label, row, 0)
            self.grid.addWidget(var_edit, row, 1)

            row += 1

        ok_button = QPushButton("OK")
        cancel_button = QPushButton("Cancel")

        self.grid.addWidget(ok_button, row, 0)
        self.grid.addWidget(cancel_button, row, 1)

        # Visa fönster

        self.show()

def create_window(var_dict):
    w = GuiWindow(var_dict)

    return w

def parse_variables(var_dict):
    """Parse variables in current context"""
    valid_vars = {}

    for key, value in var_dict.items():
        if type(value) is int:
            valid_vars[key] = value
        if type(value) is float:
            valid_vars[key] = value
        if type(value) is list:
            valid_vars[key] = value
        if type(value) is bool:
            valid_vars[key] = value

    return valid_vars

def edit_params(var_dict):
    """Run Qt event loop"""

    valid_vars = parse_variables(var_dict)
    print(valid_vars)

    app = QApplication(sys.argv)

    # Skapa vårt MyWindow objekt

    w = create_window(valid_vars)

    # Starta händelseloop

    app.exec_()

def edit_geometry(geometry):
    pass