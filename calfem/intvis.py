# -*- coding: utf-8 -*-

import sys

from PyQt5.QtWidgets import *
from PyQt5.QtGui import *

class GuiWindow(QWidget):
    def __init__(self, var_dict):
        """MyWindow constructor"""

        super().__init__()

        # Skapa gränssnittskontroller

        self.var_dict = var_dict

        self.init_gui()

    def __parse_variables(self, g):

        # a_edit = 1
        # b_slider = 2.0
        # c_list = [1, 2, 3]
        # d_check = True
        # f_param = 42.0
        # g_float = 84.0
        # g_int = 34

        row = 0

        for key, value in self.var_dict.items():
            if '_edit' in key:
                var_label = QLabel(key)
                var_edit = QLineEdit(str(value))

                g.addWidget(var_label, row, 0)
                g.addWidget(var_edit, row, 1)
            elif '_slider' in key:
                var_label = QLabel(key.split("_")[0])
                var_edit = QSlider()


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

        buttons = QHBoxLayout()

        ok_button = QPushButton("OK")
        cancel_button = QPushButton("Cancel")

        buttons.addStretch()
        buttons.addWidget(ok_button)
        buttons.addWidget(cancel_button)
        buttons.addStretch()

        self.grid.setRowStretch(row, 10)

        row += 1

        self.grid.addLayout(buttons, row, 0, 1, 2)

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