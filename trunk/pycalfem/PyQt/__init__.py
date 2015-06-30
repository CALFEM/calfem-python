#!/usr/bin/env python

# Copyright (c) 2011, Dirk Thomas, Dorian Scholz, TU Darmstadt
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#   * Redistributions in binary form must reproduce the above
#     copyright notice, this list of conditions and the following
#     disclaimer in the documentation and/or other materials provided
#     with the distribution.
#   * Neither the name of the TU Darmstadt nor the names of its
#     contributors may be used to endorse or promote products derived
#     from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""
Abstraction for different Python Qt bindings.

Supported Python Qt bindings are PyQt4 and PySide.
The Qt modules can be imported like this:

from python_qt_binding.QtCore import QObject
from python_qt_binding import QtGui, loadUi

The name of the selected binding is available in QT_BINDING.
The version of the selected binding is available in QT_BINDING_VERSION.
All available Qt modules are listed in QT_BINDING_MODULES.

The default binding order ('pyqt', 'pyside') can be overridden with a
SELECT_QT_BINDING_ORDER attribute on sys:
  setattr(sys, 'SELECT_QT_BINDING_ORDER', [FIRST_NAME, NEXT_NAME, ..])

A specific binding can be selected with a SELECT_QT_BINDING attribute on sys:
  setattr(sys, 'SELECT_QT_BINDING', MY_BINDING_NAME)
"""

import sys
from .binding_helper import loadUi, QT_BINDING, QT_BINDING_MODULES, QT_BINDING_VERSION  # @UnusedImport

# register all binding modules as sub modules of this package (python_qt_binding) for easy importing
for module_name, module in QT_BINDING_MODULES.items():
    sys.modules[__name__ + '.' + module_name] = module
    setattr(sys.modules[__name__], module_name, module)
    del module_name
    del module

del sys
