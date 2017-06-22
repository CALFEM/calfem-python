# -*- coding: iso-8859-1 -*-

import sys
import os

default_variant = 'PySide'

env_api = os.environ.get('QT_API', 'pyqt')
if '--pyside' in sys.argv:
    variant = 'PySide'
elif '--pyqt4' in sys.argv:
    variant = 'PyQt4'
elif env_api == 'pyside':
    variant = 'PySide'
elif env_api == 'pyqt':
    variant = 'PyQt4'
else:
    variant = default_variant

if variant == 'PySide':
    from PySide import QtGui, QtCore
    # This will be passed on to new versions of matplotlib
    os.environ['QT_API'] = 'pyside'
    def QtLoadUI(uifile):
        from PySide import QtUiTools
        loader = QtUiTools.QUiLoader()
        uif = QtCore.QFile(uifile)
        uif.open(QtCore.QFile.ReadOnly)
        result = loader.load(uif)
        uif.close()
        return result
elif variant == 'PyQt4':
    #import sip
    api2_classes = [
            'QData', 'QDateTime', 'QString', 'QTextStream',
            'QTime', 'QUrl', 'QVariant',
            ]
    #for cl in api2_classes:
    #    sip.setapi(cl, 2)
    from PyQt4 import QtGui, QtCore
    QtCore.Signal = QtCore.pyqtSignal
    QtCore.QString = str
    os.environ['QT_API'] = 'pyqt'
    def QtLoadUI(uifile):
        from PyQt4 import uic
        return uic.loadUi(uifile)
else:
    raise ImportError("Python Variant not specified")

__all__ = [QtGui, QtCore, QtLoadUI, variant]