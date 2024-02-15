# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'editor.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1493, 785)
        MainWindow.setStyleSheet("")
        MainWindow.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.centralwidget = QtWidgets.QWidget(MainWindow)
        MainWindow.centralwidget.setStyleSheet("QWidget{\n"
"background: rgb(220, 220, 220);\n"
"}")
        MainWindow.centralwidget.setObjectName("centralwidget")
        MainWindow.verticalLayout_2 = QtWidgets.QVBoxLayout(MainWindow.centralwidget)
        MainWindow.verticalLayout_2.setObjectName("verticalLayout_2")
        MainWindow.tabWidget = QtWidgets.QTabWidget(MainWindow.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.tabWidget.sizePolicy().hasHeightForWidth())
        MainWindow.tabWidget.setSizePolicy(sizePolicy)
        MainWindow.tabWidget.setMinimumSize(QtCore.QSize(0, 120))
        MainWindow.tabWidget.setCursor(QtGui.QCursor(QtCore.Qt.ArrowCursor))
        MainWindow.tabWidget.setStyleSheet("QTabWidget::pane { /* The tab widget frame */\n"
"    border-bottom: 6px solid;\n"
"    border-bottom-color: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, stop:0      rgb(200, 200, 200), stop:1 rgba(220,220,220));\n"
"     border-top: 8px solid rgb(0,0,128);\n"
"    /*border-radius: 4px;*/\n"
"}\n"
"\n"
"\n"
"QTabWidget::tab-bar {\n"
"    left: 20px; /* move to the right by 5px */\n"
"    bottom: -8px;\n"
"}\n"
"\n"
"QTabBar::tab {\n"
"    /*border-bottom-color:   transparent;\n"
"    border-top-left-radius: 4px;*/\n"
"    min-width: 18ex;\n"
"    padding: 4px;\n"
"}\n"
"\n"
"QTabBar::tab:selected {\n"
"    background: rgb(240,240,240);\n"
"     border: 2px solid rgb(200,200,200);\n"
"     border-bottom: 0 px solid;\n"
"     border-top-left-radius: 6px;\n"
"     border-top-right-radius: 6px;\n"
"}\n"
"\n"
"QTabBar::tab:!selected {\n"
"    border: 2px solid rgb(200,200,200);\n"
"    margin-top: 4px; /* make non-selected tabs look smaller */\n"
"    border-top-left-radius: 6px;\n"
"     border-top-right-radius: 6px;\n"
"}")
        MainWindow.tabWidget.setObjectName("tabWidget")
        MainWindow.tab_2 = QtWidgets.QWidget()
        MainWindow.tab_2.setStyleSheet("QWidget{\n"
"    background-color: rgb(240, 240, 240);\n"
"}        ")
        MainWindow.tab_2.setObjectName("tab_2")
        MainWindow.layoutWidget = QtWidgets.QWidget(MainWindow.tab_2)
        MainWindow.layoutWidget.setGeometry(QtCore.QRect(0, 0, 1100, 81))
        MainWindow.layoutWidget.setObjectName("layoutWidget")
        MainWindow.horizontalLayout_2 = QtWidgets.QHBoxLayout(MainWindow.layoutWidget)
        MainWindow.horizontalLayout_2.setContentsMargins(0, 0, 0, 0)
        MainWindow.horizontalLayout_2.setObjectName("horizontalLayout_2")
        spacerItem = QtWidgets.QSpacerItem(10, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum)
        MainWindow.horizontalLayout_2.addItem(spacerItem)
        MainWindow.gridLayout = QtWidgets.QGridLayout()
        MainWindow.gridLayout.setObjectName("gridLayout")
        MainWindow.zoomOutButtonSurface = QtWidgets.QToolButton(MainWindow.layoutWidget)
        MainWindow.zoomOutButtonSurface.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:pressed {\n"
"     background-color: rgb(185, 211, 220);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgba(185, 211, 220,100);\n"
"}\n"
"")
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(":/Images/zoom_out.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.zoomOutButtonSurface.setIcon(icon)
        MainWindow.zoomOutButtonSurface.setIconSize(QtCore.QSize(24, 24))
        MainWindow.zoomOutButtonSurface.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.zoomOutButtonSurface.setObjectName("zoomOutButtonSurface")
        MainWindow.gridLayout.addWidget(MainWindow.zoomOutButtonSurface, 0, 1, 1, 1)
        MainWindow.arrowButtonSurface = QtWidgets.QToolButton(MainWindow.layoutWidget)
        MainWindow.arrowButtonSurface.setCursor(QtGui.QCursor(QtCore.Qt.ArrowCursor))
        MainWindow.arrowButtonSurface.setFocusPolicy(QtCore.Qt.TabFocus)
        MainWindow.arrowButtonSurface.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:pressed {\n"
"     background-color: rgb(185, 211, 220);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgba(185, 211, 220,100);\n"
"}\n"
"")
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(":/Images/select.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.arrowButtonSurface.setIcon(icon1)
        MainWindow.arrowButtonSurface.setIconSize(QtCore.QSize(24, 24))
        MainWindow.arrowButtonSurface.setCheckable(False)
        MainWindow.arrowButtonSurface.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.arrowButtonSurface.setObjectName("arrowButtonSurface")
        MainWindow.gridLayout.addWidget(MainWindow.arrowButtonSurface, 0, 0, 1, 1)
        MainWindow.panningButtonSurface = QtWidgets.QToolButton(MainWindow.layoutWidget)
        MainWindow.panningButtonSurface.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:pressed {\n"
"     background-color: rgb(185, 211, 220);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgba(185, 211, 220,100);\n"
"}\n"
"")
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap(":/Images/pan.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.panningButtonSurface.setIcon(icon2)
        MainWindow.panningButtonSurface.setIconSize(QtCore.QSize(24, 24))
        MainWindow.panningButtonSurface.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.panningButtonSurface.setObjectName("panningButtonSurface")
        MainWindow.gridLayout.addWidget(MainWindow.panningButtonSurface, 1, 0, 1, 1)
        MainWindow.zoomInButtonSurface = QtWidgets.QToolButton(MainWindow.layoutWidget)
        MainWindow.zoomInButtonSurface.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:pressed {\n"
"     background-color: rgb(185, 211, 220);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgba(185, 211, 220,100);\n"
"}\n"
"")
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap(":/Images/zoom_in.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.zoomInButtonSurface.setIcon(icon3)
        MainWindow.zoomInButtonSurface.setIconSize(QtCore.QSize(24, 24))
        MainWindow.zoomInButtonSurface.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.zoomInButtonSurface.setObjectName("zoomInButtonSurface")
        MainWindow.gridLayout.addWidget(MainWindow.zoomInButtonSurface, 1, 1, 1, 1)
        MainWindow.horizontalLayout_2.addLayout(MainWindow.gridLayout)
        MainWindow.line_2 = QtWidgets.QFrame(MainWindow.layoutWidget)
        MainWindow.line_2.setMinimumSize(QtCore.QSize(15, 0))
        MainWindow.line_2.setFrameShadow(QtWidgets.QFrame.Plain)
        MainWindow.line_2.setLineWidth(1)
        MainWindow.line_2.setFrameShape(QtWidgets.QFrame.VLine)
        MainWindow.line_2.setObjectName("line_2")
        MainWindow.horizontalLayout_2.addWidget(MainWindow.line_2)
        MainWindow.verticalLayout_7 = QtWidgets.QVBoxLayout()
        MainWindow.verticalLayout_7.setObjectName("verticalLayout_7")
        MainWindow.polyButton = QtWidgets.QToolButton(MainWindow.layoutWidget)
        MainWindow.polyButton.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:pressed {\n"
"     background-color: rgb(185, 211, 220);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgba(185, 211, 220,100);\n"
"}\n"
"")
        icon4 = QtGui.QIcon()
        icon4.addPixmap(QtGui.QPixmap(":/Images/polygon.svg"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.polyButton.setIcon(icon4)
        MainWindow.polyButton.setIconSize(QtCore.QSize(24, 24))
        MainWindow.polyButton.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.polyButton.setObjectName("polyButton")
        MainWindow.verticalLayout_7.addWidget(MainWindow.polyButton)
        MainWindow.rectButton = QtWidgets.QToolButton(MainWindow.layoutWidget)
        MainWindow.rectButton.setEnabled(True)
        MainWindow.rectButton.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:pressed {\n"
"     background-color: rgb(185, 211, 220);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgba(185, 211, 220,100);\n"
"}\n"
"")
        icon5 = QtGui.QIcon()
        icon5.addPixmap(QtGui.QPixmap(":/Images/rectangle.svg"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.rectButton.setIcon(icon5)
        MainWindow.rectButton.setIconSize(QtCore.QSize(24, 24))
        MainWindow.rectButton.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.rectButton.setObjectName("rectButton")
        MainWindow.verticalLayout_7.addWidget(MainWindow.rectButton)
        MainWindow.horizontalLayout_2.addLayout(MainWindow.verticalLayout_7)
        MainWindow.verticalLayout_10 = QtWidgets.QVBoxLayout()
        MainWindow.verticalLayout_10.setObjectName("verticalLayout_10")
        MainWindow.addPolyHoleButton = QtWidgets.QToolButton(MainWindow.layoutWidget)
        MainWindow.addPolyHoleButton.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:pressed {\n"
"     background-color: rgb(185, 211, 220);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgba(185, 211, 220,100);\n"
"}\n"
"")
        icon6 = QtGui.QIcon()
        icon6.addPixmap(QtGui.QPixmap(":/Images/polygon_hole.svg"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.addPolyHoleButton.setIcon(icon6)
        MainWindow.addPolyHoleButton.setIconSize(QtCore.QSize(24, 24))
        MainWindow.addPolyHoleButton.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.addPolyHoleButton.setObjectName("addPolyHoleButton")
        MainWindow.verticalLayout_10.addWidget(MainWindow.addPolyHoleButton)
        MainWindow.addRectHoleButton = QtWidgets.QToolButton(MainWindow.layoutWidget)
        MainWindow.addRectHoleButton.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:pressed {\n"
"     background-color: rgb(185, 211, 220);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgba(185, 211, 220,100);\n"
"}\n"
"")
        icon7 = QtGui.QIcon()
        icon7.addPixmap(QtGui.QPixmap(":/Images/rectangle_hole.svg"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.addRectHoleButton.setIcon(icon7)
        MainWindow.addRectHoleButton.setIconSize(QtCore.QSize(24, 24))
        MainWindow.addRectHoleButton.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.addRectHoleButton.setObjectName("addRectHoleButton")
        MainWindow.verticalLayout_10.addWidget(MainWindow.addRectHoleButton)
        MainWindow.horizontalLayout_2.addLayout(MainWindow.verticalLayout_10)
        MainWindow.verticalLayout_12 = QtWidgets.QVBoxLayout()
        MainWindow.verticalLayout_12.setObjectName("verticalLayout_12")
        MainWindow.mergeButton = QtWidgets.QToolButton(MainWindow.layoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.mergeButton.sizePolicy().hasHeightForWidth())
        MainWindow.mergeButton.setSizePolicy(sizePolicy)
        MainWindow.mergeButton.setLayoutDirection(QtCore.Qt.LeftToRight)
        MainWindow.mergeButton.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:pressed {\n"
"     background-color: rgb(185, 211, 220);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgba(185, 211, 220,100);\n"
"}\n"
"")
        icon8 = QtGui.QIcon()
        icon8.addPixmap(QtGui.QPixmap(":/Images/merge.svg"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.mergeButton.setIcon(icon8)
        MainWindow.mergeButton.setIconSize(QtCore.QSize(24, 24))
        MainWindow.mergeButton.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.mergeButton.setAutoRaise(False)
        MainWindow.mergeButton.setObjectName("mergeButton")
        MainWindow.verticalLayout_12.addWidget(MainWindow.mergeButton)
        MainWindow.deleteButton = QtWidgets.QToolButton(MainWindow.layoutWidget)
        MainWindow.deleteButton.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:pressed {\n"
"     background-color: rgb(185, 211, 220);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgba(185, 211, 220,100);\n"
"}\n"
"")
        icon9 = QtGui.QIcon()
        icon9.addPixmap(QtGui.QPixmap(":/Images/erase.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.deleteButton.setIcon(icon9)
        MainWindow.deleteButton.setIconSize(QtCore.QSize(24, 24))
        MainWindow.deleteButton.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.deleteButton.setObjectName("deleteButton")
        MainWindow.verticalLayout_12.addWidget(MainWindow.deleteButton)
        MainWindow.horizontalLayout_2.addLayout(MainWindow.verticalLayout_12)
        MainWindow.line_5 = QtWidgets.QFrame(MainWindow.layoutWidget)
        MainWindow.line_5.setMinimumSize(QtCore.QSize(15, 0))
        MainWindow.line_5.setFrameShadow(QtWidgets.QFrame.Plain)
        MainWindow.line_5.setLineWidth(1)
        MainWindow.line_5.setFrameShape(QtWidgets.QFrame.VLine)
        MainWindow.line_5.setObjectName("line_5")
        MainWindow.horizontalLayout_2.addWidget(MainWindow.line_5)
        MainWindow.verticalLayout_9 = QtWidgets.QVBoxLayout()
        MainWindow.verticalLayout_9.setObjectName("verticalLayout_9")
        MainWindow.gridSnapButtonSurface = QtWidgets.QToolButton(MainWindow.layoutWidget)
        MainWindow.gridSnapButtonSurface.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:checked {\n"
"     background-color: rgb(191, 184, 175);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgb(223, 219, 215);\n"
"}\n"
"")
        icon10 = QtGui.QIcon()
        icon10.addPixmap(QtGui.QPixmap(":/Images/grid_snap.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.gridSnapButtonSurface.setIcon(icon10)
        MainWindow.gridSnapButtonSurface.setIconSize(QtCore.QSize(24, 24))
        MainWindow.gridSnapButtonSurface.setCheckable(True)
        MainWindow.gridSnapButtonSurface.setChecked(True)
        MainWindow.gridSnapButtonSurface.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.gridSnapButtonSurface.setObjectName("gridSnapButtonSurface")
        MainWindow.verticalLayout_9.addWidget(MainWindow.gridSnapButtonSurface)
        MainWindow.gridButtonSurface = QtWidgets.QToolButton(MainWindow.layoutWidget)
        MainWindow.gridButtonSurface.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:checked {\n"
"     background-color: rgb(191, 184, 175);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgb(223, 219, 215);\n"
"}\n"
"")
        icon11 = QtGui.QIcon()
        icon11.addPixmap(QtGui.QPixmap(":/Images/grid.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.gridButtonSurface.setIcon(icon11)
        MainWindow.gridButtonSurface.setIconSize(QtCore.QSize(24, 24))
        MainWindow.gridButtonSurface.setCheckable(True)
        MainWindow.gridButtonSurface.setChecked(True)
        MainWindow.gridButtonSurface.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.gridButtonSurface.setObjectName("gridButtonSurface")
        MainWindow.verticalLayout_9.addWidget(MainWindow.gridButtonSurface)
        MainWindow.horizontalLayout_2.addLayout(MainWindow.verticalLayout_9)
        MainWindow.verticalLayout_8 = QtWidgets.QVBoxLayout()
        MainWindow.verticalLayout_8.setObjectName("verticalLayout_8")
        MainWindow.gridSpacingSpinBox = QtWidgets.QSpinBox(MainWindow.layoutWidget)
        MainWindow.gridSpacingSpinBox.setMinimum(5)
        MainWindow.gridSpacingSpinBox.setMaximum(1000)
        MainWindow.gridSpacingSpinBox.setSingleStep(10)
        MainWindow.gridSpacingSpinBox.setProperty("value", 20)
        MainWindow.gridSpacingSpinBox.setObjectName("gridSpacingSpinBox")
        MainWindow.verticalLayout_8.addWidget(MainWindow.gridSpacingSpinBox)
        MainWindow.gridSpacingButton = QtWidgets.QToolButton(MainWindow.layoutWidget)
        MainWindow.gridSpacingButton.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:pressed {\n"
"     background-color: rgb(191, 184, 175);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgb(223, 219, 215);\n"
"}\n"
"")
        icon12 = QtGui.QIcon()
        icon12.addPixmap(QtGui.QPixmap(":/Images/grid_settings.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.gridSpacingButton.setIcon(icon12)
        MainWindow.gridSpacingButton.setIconSize(QtCore.QSize(24, 24))
        MainWindow.gridSpacingButton.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.gridSpacingButton.setObjectName("gridSpacingButton")
        MainWindow.verticalLayout_8.addWidget(MainWindow.gridSpacingButton)
        MainWindow.horizontalLayout_2.addLayout(MainWindow.verticalLayout_8)
        MainWindow.line_6 = QtWidgets.QFrame(MainWindow.layoutWidget)
        MainWindow.line_6.setMinimumSize(QtCore.QSize(15, 0))
        MainWindow.line_6.setFrameShadow(QtWidgets.QFrame.Plain)
        MainWindow.line_6.setLineWidth(1)
        MainWindow.line_6.setFrameShape(QtWidgets.QFrame.VLine)
        MainWindow.line_6.setObjectName("line_6")
        MainWindow.horizontalLayout_2.addWidget(MainWindow.line_6)
        MainWindow.verticalLayout_11 = QtWidgets.QVBoxLayout()
        MainWindow.verticalLayout_11.setObjectName("verticalLayout_11")
        MainWindow.loadGeometryButton = QtWidgets.QToolButton(MainWindow.layoutWidget)
        MainWindow.loadGeometryButton.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:pressed {\n"
"     background-color: rgb(185, 211, 220);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgba(185, 211, 220,100);\n"
"}\n"
"")
        icon13 = QtGui.QIcon()
        icon13.addPixmap(QtGui.QPixmap(":/Images/open.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.loadGeometryButton.setIcon(icon13)
        MainWindow.loadGeometryButton.setIconSize(QtCore.QSize(24, 24))
        MainWindow.loadGeometryButton.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.loadGeometryButton.setObjectName("loadGeometryButton")
        MainWindow.verticalLayout_11.addWidget(MainWindow.loadGeometryButton)
        MainWindow.horizontalLayout_2.addLayout(MainWindow.verticalLayout_11)
        MainWindow.tabWidget.addTab(MainWindow.tab_2, "")
        MainWindow.tab = QtWidgets.QWidget()
        MainWindow.tab.setStyleSheet("QWidget{\n"
"    background-color: rgb(240, 240, 240);\n"
"}        ")
        MainWindow.tab.setObjectName("tab")
        MainWindow.horizontalLayoutWidget_2 = QtWidgets.QWidget(MainWindow.tab)
        MainWindow.horizontalLayoutWidget_2.setGeometry(QtCore.QRect(0, 0, 464, 81))
        MainWindow.horizontalLayoutWidget_2.setObjectName("horizontalLayoutWidget_2")
        MainWindow.horizontalLayout_3 = QtWidgets.QHBoxLayout(MainWindow.horizontalLayoutWidget_2)
        MainWindow.horizontalLayout_3.setContentsMargins(0, 0, 0, 0)
        MainWindow.horizontalLayout_3.setObjectName("horizontalLayout_3")
        spacerItem1 = QtWidgets.QSpacerItem(10, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum)
        MainWindow.horizontalLayout_3.addItem(spacerItem1)
        MainWindow.gridLayout_2 = QtWidgets.QGridLayout()
        MainWindow.gridLayout_2.setObjectName("gridLayout_2")
        MainWindow.arrowButtonBorder = QtWidgets.QToolButton(MainWindow.horizontalLayoutWidget_2)
        MainWindow.arrowButtonBorder.setCursor(QtGui.QCursor(QtCore.Qt.ArrowCursor))
        MainWindow.arrowButtonBorder.setFocusPolicy(QtCore.Qt.TabFocus)
        MainWindow.arrowButtonBorder.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:pressed {\n"
"     background-color: rgb(185, 211, 220);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgba(185, 211, 220,100);\n"
"}\n"
"")
        MainWindow.arrowButtonBorder.setIcon(icon1)
        MainWindow.arrowButtonBorder.setIconSize(QtCore.QSize(24, 24))
        MainWindow.arrowButtonBorder.setCheckable(True)
        MainWindow.arrowButtonBorder.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.arrowButtonBorder.setObjectName("arrowButtonBorder")
        MainWindow.gridLayout_2.addWidget(MainWindow.arrowButtonBorder, 0, 0, 1, 1)
        MainWindow.zoomOutButtonBorder = QtWidgets.QToolButton(MainWindow.horizontalLayoutWidget_2)
        MainWindow.zoomOutButtonBorder.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:pressed {\n"
"     background-color: rgb(185, 211, 220);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgba(185, 211, 220,100);\n"
"}\n"
"")
        MainWindow.zoomOutButtonBorder.setIcon(icon)
        MainWindow.zoomOutButtonBorder.setIconSize(QtCore.QSize(24, 24))
        MainWindow.zoomOutButtonBorder.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.zoomOutButtonBorder.setObjectName("zoomOutButtonBorder")
        MainWindow.gridLayout_2.addWidget(MainWindow.zoomOutButtonBorder, 0, 1, 1, 1)
        MainWindow.zoomInButtonBorder = QtWidgets.QToolButton(MainWindow.horizontalLayoutWidget_2)
        MainWindow.zoomInButtonBorder.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:pressed {\n"
"     background-color: rgb(185, 211, 220);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgba(185, 211, 220,100);\n"
"}\n"
"")
        MainWindow.zoomInButtonBorder.setIcon(icon3)
        MainWindow.zoomInButtonBorder.setIconSize(QtCore.QSize(24, 24))
        MainWindow.zoomInButtonBorder.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.zoomInButtonBorder.setObjectName("zoomInButtonBorder")
        MainWindow.gridLayout_2.addWidget(MainWindow.zoomInButtonBorder, 1, 1, 1, 1)
        MainWindow.panningButtonBorder = QtWidgets.QToolButton(MainWindow.horizontalLayoutWidget_2)
        MainWindow.panningButtonBorder.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:pressed {\n"
"     background-color: rgb(185, 211, 220);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgba(185, 211, 220,100);\n"
"}\n"
"")
        MainWindow.panningButtonBorder.setIcon(icon2)
        MainWindow.panningButtonBorder.setIconSize(QtCore.QSize(24, 24))
        MainWindow.panningButtonBorder.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.panningButtonBorder.setObjectName("panningButtonBorder")
        MainWindow.gridLayout_2.addWidget(MainWindow.panningButtonBorder, 1, 0, 1, 1)
        MainWindow.horizontalLayout_3.addLayout(MainWindow.gridLayout_2)
        MainWindow.line = QtWidgets.QFrame(MainWindow.horizontalLayoutWidget_2)
        MainWindow.line.setMinimumSize(QtCore.QSize(15, 0))
        MainWindow.line.setFrameShadow(QtWidgets.QFrame.Plain)
        MainWindow.line.setFrameShape(QtWidgets.QFrame.VLine)
        MainWindow.line.setObjectName("line")
        MainWindow.horizontalLayout_3.addWidget(MainWindow.line)
        MainWindow.verticalLayout_6 = QtWidgets.QVBoxLayout()
        MainWindow.verticalLayout_6.setObjectName("verticalLayout_6")
        MainWindow.setMarkerButton = QtWidgets.QToolButton(MainWindow.horizontalLayoutWidget_2)
        MainWindow.setMarkerButton.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:pressed {\n"
"     background-color: rgb(185, 211, 220);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgba(185, 211, 220,100);\n"
"}\n"
"")
        icon14 = QtGui.QIcon()
        icon14.addPixmap(QtGui.QPixmap(":/Images/add_marker.svg"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.setMarkerButton.setIcon(icon14)
        MainWindow.setMarkerButton.setIconSize(QtCore.QSize(24, 24))
        MainWindow.setMarkerButton.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.setMarkerButton.setObjectName("setMarkerButton")
        MainWindow.verticalLayout_6.addWidget(MainWindow.setMarkerButton)
        MainWindow.splitEdgeButton = QtWidgets.QToolButton(MainWindow.horizontalLayoutWidget_2)
        MainWindow.splitEdgeButton.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:pressed {\n"
"     background-color: rgb(185, 211, 220);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgba(185, 211, 220,100);\n"
"}\n"
"")
        icon15 = QtGui.QIcon()
        icon15.addPixmap(QtGui.QPixmap(":/Images/split_edge.svg"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.splitEdgeButton.setIcon(icon15)
        MainWindow.splitEdgeButton.setIconSize(QtCore.QSize(24, 24))
        MainWindow.splitEdgeButton.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.splitEdgeButton.setObjectName("splitEdgeButton")
        MainWindow.verticalLayout_6.addWidget(MainWindow.splitEdgeButton)
        MainWindow.horizontalLayout_3.addLayout(MainWindow.verticalLayout_6)
        MainWindow.line_3 = QtWidgets.QFrame(MainWindow.horizontalLayoutWidget_2)
        MainWindow.line_3.setFrameShadow(QtWidgets.QFrame.Plain)
        MainWindow.line_3.setFrameShape(QtWidgets.QFrame.VLine)
        MainWindow.line_3.setObjectName("line_3")
        MainWindow.horizontalLayout_3.addWidget(MainWindow.line_3)
        MainWindow.verticalLayout_13 = QtWidgets.QVBoxLayout()
        MainWindow.verticalLayout_13.setObjectName("verticalLayout_13")
        MainWindow.gridSnapButtonBorder = QtWidgets.QToolButton(MainWindow.horizontalLayoutWidget_2)
        MainWindow.gridSnapButtonBorder.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:checked {\n"
"     background-color: rgb(191, 184, 175);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgb(223, 219, 215);\n"
"}\n"
"")
        MainWindow.gridSnapButtonBorder.setIcon(icon10)
        MainWindow.gridSnapButtonBorder.setIconSize(QtCore.QSize(24, 24))
        MainWindow.gridSnapButtonBorder.setCheckable(True)
        MainWindow.gridSnapButtonBorder.setChecked(True)
        MainWindow.gridSnapButtonBorder.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.gridSnapButtonBorder.setObjectName("gridSnapButtonBorder")
        MainWindow.verticalLayout_13.addWidget(MainWindow.gridSnapButtonBorder)
        MainWindow.gridButtonBorder = QtWidgets.QToolButton(MainWindow.horizontalLayoutWidget_2)
        MainWindow.gridButtonBorder.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:checked {\n"
"     background-color: rgb(191, 184, 175);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgb(223, 219, 215);\n"
"}\n"
"")
        MainWindow.gridButtonBorder.setIcon(icon11)
        MainWindow.gridButtonBorder.setIconSize(QtCore.QSize(24, 24))
        MainWindow.gridButtonBorder.setCheckable(True)
        MainWindow.gridButtonBorder.setChecked(True)
        MainWindow.gridButtonBorder.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.gridButtonBorder.setObjectName("gridButtonBorder")
        MainWindow.verticalLayout_13.addWidget(MainWindow.gridButtonBorder)
        MainWindow.horizontalLayout_3.addLayout(MainWindow.verticalLayout_13)
        MainWindow.tabWidget.addTab(MainWindow.tab, "")
        MainWindow.tab_4 = QtWidgets.QWidget()
        MainWindow.tab_4.setStyleSheet("QWidget{\n"
"    background-color: rgb(240, 240, 240);\n"
"}        ")
        MainWindow.tab_4.setObjectName("tab_4")
        MainWindow.saveGeomButton = QtWidgets.QToolButton(MainWindow.tab_4)
        MainWindow.saveGeomButton.setGeometry(QtCore.QRect(10, 30, 141, 21))
        MainWindow.saveGeomButton.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:pressed {\n"
"     background-color: rgb(185, 211, 220);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgba(185, 211, 220,100);\n"
"}\n"
"")
        icon16 = QtGui.QIcon()
        icon16.addPixmap(QtGui.QPixmap(":/Images/save.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.saveGeomButton.setIcon(icon16)
        MainWindow.saveGeomButton.setIconSize(QtCore.QSize(24, 24))
        MainWindow.saveGeomButton.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.saveGeomButton.setObjectName("saveGeomButton")
        MainWindow.tabWidget.addTab(MainWindow.tab_4, "")
        MainWindow.tab_3 = QtWidgets.QWidget()
        MainWindow.tab_3.setStyleSheet("QWidget{\n"
"    background-color: rgb(240, 240, 240);\n"
"}        ")
        MainWindow.tab_3.setObjectName("tab_3")
        MainWindow.horizontalLayoutWidget = QtWidgets.QWidget(MainWindow.tab_3)
        MainWindow.horizontalLayoutWidget.setGeometry(QtCore.QRect(10, 0, 351, 71))
        MainWindow.horizontalLayoutWidget.setObjectName("horizontalLayoutWidget")
        MainWindow.horizontalLayout_5 = QtWidgets.QHBoxLayout(MainWindow.horizontalLayoutWidget)
        MainWindow.horizontalLayout_5.setContentsMargins(0, 0, 0, 0)
        MainWindow.horizontalLayout_5.setSpacing(0)
        MainWindow.horizontalLayout_5.setObjectName("horizontalLayout_5")
        MainWindow.label_9 = QtWidgets.QLabel(MainWindow.horizontalLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.label_9.sizePolicy().hasHeightForWidth())
        MainWindow.label_9.setSizePolicy(sizePolicy)
        MainWindow.label_9.setText("")
        MainWindow.label_9.setPixmap(QtGui.QPixmap(":/Images/save.png"))
        MainWindow.label_9.setObjectName("label_9")
        MainWindow.horizontalLayout_5.addWidget(MainWindow.label_9)
        MainWindow.line_4 = QtWidgets.QFrame(MainWindow.horizontalLayoutWidget)
        MainWindow.line_4.setFrameShape(QtWidgets.QFrame.VLine)
        MainWindow.line_4.setFrameShadow(QtWidgets.QFrame.Sunken)
        MainWindow.line_4.setObjectName("line_4")
        MainWindow.horizontalLayout_5.addWidget(MainWindow.line_4)
        MainWindow.verticalLayout_19 = QtWidgets.QVBoxLayout()
        MainWindow.verticalLayout_19.setSpacing(0)
        MainWindow.verticalLayout_19.setObjectName("verticalLayout_19")
        MainWindow.saveMeshButton = QtWidgets.QToolButton(MainWindow.horizontalLayoutWidget)
        MainWindow.saveMeshButton.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:pressed {\n"
"     background-color: rgb(185, 211, 220);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgba(185, 211, 220,100);\n"
"}\n"
"")
        icon17 = QtGui.QIcon()
        icon17.addPixmap(QtGui.QPixmap(":/Images/blank.svg"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.saveMeshButton.setIcon(icon17)
        MainWindow.saveMeshButton.setIconSize(QtCore.QSize(12, 12))
        MainWindow.saveMeshButton.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.saveMeshButton.setObjectName("saveMeshButton")
        MainWindow.verticalLayout_19.addWidget(MainWindow.saveMeshButton)
        MainWindow.saveArraysButton = QtWidgets.QToolButton(MainWindow.horizontalLayoutWidget)
        MainWindow.saveArraysButton.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:pressed {\n"
"     background-color: rgb(185, 211, 220);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgba(185, 211, 220,100);\n"
"}\n"
"")
        MainWindow.saveArraysButton.setIcon(icon17)
        MainWindow.saveArraysButton.setIconSize(QtCore.QSize(12, 12))
        MainWindow.saveArraysButton.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.saveArraysButton.setObjectName("saveArraysButton")
        MainWindow.verticalLayout_19.addWidget(MainWindow.saveArraysButton)
        MainWindow.saveArraysMatlabButton = QtWidgets.QToolButton(MainWindow.horizontalLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.saveArraysMatlabButton.sizePolicy().hasHeightForWidth())
        MainWindow.saveArraysMatlabButton.setSizePolicy(sizePolicy)
        MainWindow.saveArraysMatlabButton.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:pressed {\n"
"     background-color: rgb(185, 211, 220);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgba(185, 211, 220,100);\n"
"}\n"
"")
        MainWindow.saveArraysMatlabButton.setIcon(icon17)
        MainWindow.saveArraysMatlabButton.setIconSize(QtCore.QSize(12, 12))
        MainWindow.saveArraysMatlabButton.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        MainWindow.saveArraysMatlabButton.setObjectName("saveArraysMatlabButton")
        MainWindow.verticalLayout_19.addWidget(MainWindow.saveArraysMatlabButton)
        MainWindow.horizontalLayout_5.addLayout(MainWindow.verticalLayout_19)
        MainWindow.tabWidget.addTab(MainWindow.tab_3, "")
        MainWindow.verticalLayout_2.addWidget(MainWindow.tabWidget)
        MainWindow.horizontalLayout_4 = QtWidgets.QHBoxLayout()
        MainWindow.horizontalLayout_4.setObjectName("horizontalLayout_4")
        spacerItem2 = QtWidgets.QSpacerItem(20, 20, QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        MainWindow.horizontalLayout_4.addItem(spacerItem2)
        MainWindow.labelTooltip = QtWidgets.QLabel(MainWindow.centralwidget)
        MainWindow.labelTooltip.setStyleSheet("QLabel{\n"
"    background: rgba(255,255,255,0);\n"
"}")
        MainWindow.labelTooltip.setText("")
        MainWindow.labelTooltip.setObjectName("labelTooltip")
        MainWindow.horizontalLayout_4.addWidget(MainWindow.labelTooltip)
        spacerItem3 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        MainWindow.horizontalLayout_4.addItem(spacerItem3)
        MainWindow.label = QtWidgets.QLabel(MainWindow.centralwidget)
        MainWindow.label.setStyleSheet("QLabel{\n"
"    background: rgba(255,255,255,0);\n"
"}")
        MainWindow.label.setObjectName("label")
        MainWindow.horizontalLayout_4.addWidget(MainWindow.label)
        MainWindow.labelX = QtWidgets.QLabel(MainWindow.centralwidget)
        MainWindow.labelX.setStyleSheet("QLabel{\n"
"    background: rgba(255,255,255,0);\n"
"}")
        MainWindow.labelX.setObjectName("labelX")
        MainWindow.horizontalLayout_4.addWidget(MainWindow.labelX)
        MainWindow.label_3 = QtWidgets.QLabel(MainWindow.centralwidget)
        MainWindow.label_3.setStyleSheet("QLabel{\n"
"    background: rgba(255,255,255,0);\n"
"}")
        MainWindow.label_3.setObjectName("label_3")
        MainWindow.horizontalLayout_4.addWidget(MainWindow.label_3)
        MainWindow.labelY = QtWidgets.QLabel(MainWindow.centralwidget)
        MainWindow.labelY.setStyleSheet("QLabel{\n"
"    background: rgba(255,255,255,0);\n"
"}")
        MainWindow.labelY.setObjectName("labelY")
        MainWindow.horizontalLayout_4.addWidget(MainWindow.labelY)
        spacerItem4 = QtWidgets.QSpacerItem(220, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum)
        MainWindow.horizontalLayout_4.addItem(spacerItem4)
        MainWindow.verticalLayout_2.addLayout(MainWindow.horizontalLayout_4)
        MainWindow.horizontalLayout = QtWidgets.QHBoxLayout()
        MainWindow.horizontalLayout.setSpacing(0)
        MainWindow.horizontalLayout.setObjectName("horizontalLayout")
        MainWindow.stackedWidget = QtWidgets.QStackedWidget(MainWindow.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.stackedWidget.sizePolicy().hasHeightForWidth())
        MainWindow.stackedWidget.setSizePolicy(sizePolicy)
        MainWindow.stackedWidget.setMinimumSize(QtCore.QSize(661, 420))
        MainWindow.stackedWidget.setStyleSheet("QStackedWidget{\n"
"    background: rgba(255,255,255,0);\n"
"}")
        MainWindow.stackedWidget.setObjectName("stackedWidget")
        MainWindow.page = QtWidgets.QWidget()
        MainWindow.page.setStyleSheet("QWidget{\n"
"    background: rgba(255,255,255,0);\n"
"}")
        MainWindow.page.setObjectName("page")
        MainWindow.verticalLayout_4 = QtWidgets.QVBoxLayout(MainWindow.page)
        MainWindow.verticalLayout_4.setObjectName("verticalLayout_4")
        MainWindow.verticalLayout = QtWidgets.QVBoxLayout()
        MainWindow.verticalLayout.setObjectName("verticalLayout")
        MainWindow.verticalLayout_4.addLayout(MainWindow.verticalLayout)
        MainWindow.stackedWidget.addWidget(MainWindow.page)
        MainWindow.page_2 = QtWidgets.QWidget()
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.page_2.sizePolicy().hasHeightForWidth())
        MainWindow.page_2.setSizePolicy(sizePolicy)
        MainWindow.page_2.setStyleSheet("QWidget{\n"
"    background: rgba(255,255,255,0);\n"
"}")
        MainWindow.page_2.setObjectName("page_2")
        MainWindow.verticalLayout_3 = QtWidgets.QVBoxLayout(MainWindow.page_2)
        MainWindow.verticalLayout_3.setObjectName("verticalLayout_3")
        MainWindow.graphicsView = QtWidgets.QGraphicsView(MainWindow.page_2)
        MainWindow.graphicsView.setMinimumSize(QtCore.QSize(661, 401))
        MainWindow.graphicsView.viewport().setProperty("cursor", QtGui.QCursor(QtCore.Qt.ArrowCursor))
        MainWindow.graphicsView.setMouseTracking(True)
        MainWindow.graphicsView.setStyleSheet("QGraphicsView{\n"
"    background: rgb(255,255,255);\n"
"    border: 2px solid rgb(0,0,128);\n"
"}")
        MainWindow.graphicsView.setObjectName("graphicsView")
        MainWindow.verticalLayout_3.addWidget(MainWindow.graphicsView)
        MainWindow.stackedWidget.addWidget(MainWindow.page_2)
        MainWindow.horizontalLayout.addWidget(MainWindow.stackedWidget)
        MainWindow.stackedWidgetRight = QtWidgets.QStackedWidget(MainWindow.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.stackedWidgetRight.sizePolicy().hasHeightForWidth())
        MainWindow.stackedWidgetRight.setSizePolicy(sizePolicy)
        MainWindow.stackedWidgetRight.setMinimumSize(QtCore.QSize(200, 401))
        MainWindow.stackedWidgetRight.setStyleSheet("QStackedWidget{\n"
"    background: rgba(255,255,255,0);\n"
"}")
        MainWindow.stackedWidgetRight.setObjectName("stackedWidgetRight")
        MainWindow.page_5 = QtWidgets.QWidget()
        MainWindow.page_5.setStyleSheet("QWidget{\n"
"    background: rgba(255,255,255,0);\n"
"}")
        MainWindow.page_5.setObjectName("page_5")
        MainWindow.verticalLayout_14 = QtWidgets.QVBoxLayout(MainWindow.page_5)
        MainWindow.verticalLayout_14.setObjectName("verticalLayout_14")
        MainWindow.verticalLayout_15 = QtWidgets.QVBoxLayout()
        MainWindow.verticalLayout_15.setObjectName("verticalLayout_15")
        MainWindow.widget = QtWidgets.QWidget(MainWindow.page_5)
        MainWindow.widget.setStyleSheet("QWidget{\n"
"    background: rgb(255,255,255);\n"
"    border: 2px solid rgb(0,0,128);\n"
"}")
        MainWindow.widget.setObjectName("widget")
        MainWindow.gridLayout_3 = QtWidgets.QGridLayout(MainWindow.widget)
        MainWindow.gridLayout_3.setObjectName("gridLayout_3")
        MainWindow.label_4 = QtWidgets.QLabel(MainWindow.widget)
        MainWindow.label_4.setStyleSheet("QLabel{\n"
"    background: rgb(255,255,255);\n"
"    border: 0px solid rgb(0,0,128);\n"
"}")
        MainWindow.label_4.setObjectName("label_4")
        MainWindow.gridLayout_3.addWidget(MainWindow.label_4, 5, 0, 1, 1)
        MainWindow.elSizeSpinBox = QtWidgets.QDoubleSpinBox(MainWindow.widget)
        MainWindow.elSizeSpinBox.setStyleSheet("QDoubleSpinBox{\n"
"    background: rgb(255,255,255);\n"
"    border: 1px solid rgb(0,0,128);\n"
"}")
        MainWindow.elSizeSpinBox.setDecimals(1)
        MainWindow.elSizeSpinBox.setMinimum(0.0)
        MainWindow.elSizeSpinBox.setMaximum(1000.0)
        MainWindow.elSizeSpinBox.setSingleStep(5.0)
        MainWindow.elSizeSpinBox.setProperty("value", 25.0)
        MainWindow.elSizeSpinBox.setObjectName("elSizeSpinBox")
        MainWindow.gridLayout_3.addWidget(MainWindow.elSizeSpinBox, 5, 1, 1, 1)
        MainWindow.radioButtonQuadrangle = QtWidgets.QRadioButton(MainWindow.widget)
        MainWindow.radioButtonQuadrangle.setStyleSheet("QRadioButton{\n"
"    background: rgb(255,255,255);\n"
"    border: 0px solid rgb(0,0,128);\n"
"}")
        MainWindow.radioButtonQuadrangle.setChecked(False)
        MainWindow.radioButtonQuadrangle.setObjectName("radioButtonQuadrangle")
        MainWindow.gridLayout_3.addWidget(MainWindow.radioButtonQuadrangle, 2, 1, 1, 1)
        spacerItem5 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        MainWindow.gridLayout_3.addItem(spacerItem5, 10, 0, 1, 1)
        MainWindow.label_5 = QtWidgets.QLabel(MainWindow.widget)
        MainWindow.label_5.setStyleSheet("QLabel{\n"
"    background: rgb(255,255,255);\n"
"    border: 0px solid rgb(0,0,128);\n"
"}")
        MainWindow.label_5.setObjectName("label_5")
        MainWindow.gridLayout_3.addWidget(MainWindow.label_5, 7, 0, 1, 1)
        MainWindow.DOFSpinBox = QtWidgets.QSpinBox(MainWindow.widget)
        MainWindow.DOFSpinBox.setStyleSheet("QSpinBox{\n"
"    background: rgb(255,255,255);\n"
"    border: 1px solid rgb(0,0,128);\n"
"}")
        MainWindow.DOFSpinBox.setMinimum(1)
        MainWindow.DOFSpinBox.setObjectName("DOFSpinBox")
        MainWindow.gridLayout_3.addWidget(MainWindow.DOFSpinBox, 7, 1, 1, 1)
        MainWindow.radioButtonTriangle = QtWidgets.QRadioButton(MainWindow.widget)
        MainWindow.radioButtonTriangle.setStyleSheet("QRadioButton{\n"
"    background: rgb(255,255,255);\n"
"    border: 0px solid rgb(0,0,128);\n"
"}")
        MainWindow.radioButtonTriangle.setChecked(True)
        MainWindow.radioButtonTriangle.setObjectName("radioButtonTriangle")
        MainWindow.gridLayout_3.addWidget(MainWindow.radioButtonTriangle, 0, 1, 1, 1)
        spacerItem6 = QtWidgets.QSpacerItem(20, 10, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        MainWindow.gridLayout_3.addItem(spacerItem6, 4, 1, 1, 1)
        spacerItem7 = QtWidgets.QSpacerItem(20, 10, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        MainWindow.gridLayout_3.addItem(spacerItem7, 6, 1, 1, 1)
        MainWindow.label_2 = QtWidgets.QLabel(MainWindow.widget)
        MainWindow.label_2.setStyleSheet("QLabel{\n"
"    background: rgb(255,255,255);\n"
"    border: 0px solid rgb(0,0,128);\n"
"}")
        MainWindow.label_2.setObjectName("label_2")
        MainWindow.gridLayout_3.addWidget(MainWindow.label_2, 0, 0, 1, 1)
        spacerItem8 = QtWidgets.QSpacerItem(20, 10, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        MainWindow.gridLayout_3.addItem(spacerItem8, 8, 1, 1, 1)
        MainWindow.refreshMeshButton = QtWidgets.QPushButton(MainWindow.widget)
        MainWindow.refreshMeshButton.setObjectName("refreshMeshButton")
        MainWindow.gridLayout_3.addWidget(MainWindow.refreshMeshButton, 9, 1, 1, 1)
        MainWindow.verticalLayout_15.addWidget(MainWindow.widget)
        MainWindow.verticalLayout_14.addLayout(MainWindow.verticalLayout_15)
        MainWindow.stackedWidgetRight.addWidget(MainWindow.page_5)
        MainWindow.page_6 = QtWidgets.QWidget()
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.page_6.sizePolicy().hasHeightForWidth())
        MainWindow.page_6.setSizePolicy(sizePolicy)
        MainWindow.page_6.setStyleSheet("QWidget{\n"
"    background: rgba(255,255,255,0);\n"
"}")
        MainWindow.page_6.setObjectName("page_6")
        MainWindow.verticalLayout_16 = QtWidgets.QVBoxLayout(MainWindow.page_6)
        MainWindow.verticalLayout_16.setObjectName("verticalLayout_16")
        MainWindow.scrollArea = QtWidgets.QScrollArea(MainWindow.page_6)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.scrollArea.sizePolicy().hasHeightForWidth())
        MainWindow.scrollArea.setSizePolicy(sizePolicy)
        MainWindow.scrollArea.setMinimumSize(QtCore.QSize(200, 300))
        MainWindow.scrollArea.setStyleSheet("QScrollArea{\n"
"    background: rgb(255,255,255);\n"
"    border: 2px solid rgb(0,0,128);\n"
"}")
        MainWindow.scrollArea.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustIgnored)
        MainWindow.scrollArea.setWidgetResizable(True)
        MainWindow.scrollArea.setObjectName("scrollArea")
        MainWindow.scrollAreaWidgetContents_2 = QtWidgets.QWidget()
        MainWindow.scrollAreaWidgetContents_2.setGeometry(QtCore.QRect(0, 0, 196, 296))
        MainWindow.scrollAreaWidgetContents_2.setStyleSheet("QWidget{\n"
"    background: rgb(255,255,255);\n"
"}")
        MainWindow.scrollAreaWidgetContents_2.setObjectName("scrollAreaWidgetContents_2")
        MainWindow.scrollArea.setWidget(MainWindow.scrollAreaWidgetContents_2)
        MainWindow.verticalLayout_16.addWidget(MainWindow.scrollArea)
        spacerItem9 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        MainWindow.verticalLayout_16.addItem(spacerItem9)
        MainWindow.stackedWidgetRight.addWidget(MainWindow.page_6)
        MainWindow.page_7 = QtWidgets.QWidget()
        MainWindow.page_7.setObjectName("page_7")
        MainWindow.verticalLayout_5 = QtWidgets.QVBoxLayout(MainWindow.page_7)
        MainWindow.verticalLayout_5.setObjectName("verticalLayout_5")
        MainWindow.widget_2 = QtWidgets.QWidget(MainWindow.page_7)
        MainWindow.widget_2.setStyleSheet("QWidget{\n"
"    background: rgb(255,255,255);\n"
"    border: 2px solid rgb(0,0,128);\n"
"}")
        MainWindow.widget_2.setObjectName("widget_2")
        MainWindow.gridLayout_4 = QtWidgets.QGridLayout(MainWindow.widget_2)
        MainWindow.gridLayout_4.setObjectName("gridLayout_4")
        spacerItem10 = QtWidgets.QSpacerItem(20, 10, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        MainWindow.gridLayout_4.addItem(spacerItem10, 2, 1, 1, 1)
        MainWindow.radioButtonDisplayPointLabelsTrue = QtWidgets.QRadioButton(MainWindow.widget_2)
        MainWindow.radioButtonDisplayPointLabelsTrue.setStyleSheet("QRadioButton{\n"
"    background: rgb(255,255,255);\n"
"    border: 0px solid rgb(0,0,128);\n"
"}")
        MainWindow.radioButtonDisplayPointLabelsTrue.setChecked(True)
        MainWindow.radioButtonDisplayPointLabelsTrue.setObjectName("radioButtonDisplayPointLabelsTrue")
        MainWindow.buttonGroupDisplayPointLabels = QtWidgets.QButtonGroup(MainWindow)
        MainWindow.buttonGroupDisplayPointLabels.setObjectName("buttonGroupDisplayPointLabels")
        MainWindow.buttonGroupDisplayPointLabels.addButton(MainWindow.radioButtonDisplayPointLabelsTrue)
        MainWindow.gridLayout_4.addWidget(MainWindow.radioButtonDisplayPointLabelsTrue, 3, 1, 1, 1)
        MainWindow.radioButtonDisplayEdgeLabelsTrue = QtWidgets.QRadioButton(MainWindow.widget_2)
        MainWindow.radioButtonDisplayEdgeLabelsTrue.setStyleSheet("QRadioButton{\n"
"    background: rgb(255,255,255);\n"
"    border: 0px solid rgb(0,0,128);\n"
"}")
        MainWindow.radioButtonDisplayEdgeLabelsTrue.setChecked(True)
        MainWindow.radioButtonDisplayEdgeLabelsTrue.setObjectName("radioButtonDisplayEdgeLabelsTrue")
        MainWindow.buttonGroupDOFsPerNode = QtWidgets.QButtonGroup(MainWindow)
        MainWindow.buttonGroupDOFsPerNode.setObjectName("buttonGroupDOFsPerNode")
        MainWindow.buttonGroupDOFsPerNode.addButton(MainWindow.radioButtonDisplayEdgeLabelsTrue)
        MainWindow.gridLayout_4.addWidget(MainWindow.radioButtonDisplayEdgeLabelsTrue, 6, 1, 1, 1)
        MainWindow.radioButtonDisplayPointLabelsFalse = QtWidgets.QRadioButton(MainWindow.widget_2)
        MainWindow.radioButtonDisplayPointLabelsFalse.setStyleSheet("QRadioButton{\n"
"    background: rgb(255,255,255);\n"
"    border: 0px solid rgb(0,0,128);\n"
"}")
        MainWindow.radioButtonDisplayPointLabelsFalse.setObjectName("radioButtonDisplayPointLabelsFalse")
        MainWindow.buttonGroupDisplayPointLabels.addButton(MainWindow.radioButtonDisplayPointLabelsFalse)
        MainWindow.gridLayout_4.addWidget(MainWindow.radioButtonDisplayPointLabelsFalse, 4, 1, 1, 1)
        MainWindow.radioButtonDrawPointFalse = QtWidgets.QRadioButton(MainWindow.widget_2)
        MainWindow.radioButtonDrawPointFalse.setStyleSheet("QRadioButton{\n"
"    background: rgb(255,255,255);\n"
"    border: 0px solid rgb(0,0,128);\n"
"}")
        MainWindow.radioButtonDrawPointFalse.setObjectName("radioButtonDrawPointFalse")
        MainWindow.buttonGroupDrawPoints = QtWidgets.QButtonGroup(MainWindow)
        MainWindow.buttonGroupDrawPoints.setObjectName("buttonGroupDrawPoints")
        MainWindow.buttonGroupDrawPoints.addButton(MainWindow.radioButtonDrawPointFalse)
        MainWindow.gridLayout_4.addWidget(MainWindow.radioButtonDrawPointFalse, 1, 1, 1, 1)
        MainWindow.radioButtonDrawPointTrue = QtWidgets.QRadioButton(MainWindow.widget_2)
        MainWindow.radioButtonDrawPointTrue.setStyleSheet("QRadioButton{\n"
"    background: rgb(255,255,255);\n"
"    border: 0px solid rgb(0,0,128);\n"
"}")
        MainWindow.radioButtonDrawPointTrue.setChecked(True)
        MainWindow.radioButtonDrawPointTrue.setObjectName("radioButtonDrawPointTrue")
        MainWindow.buttonGroupDrawPoints.addButton(MainWindow.radioButtonDrawPointTrue)
        MainWindow.gridLayout_4.addWidget(MainWindow.radioButtonDrawPointTrue, 0, 1, 1, 1)
        MainWindow.label_8 = QtWidgets.QLabel(MainWindow.widget_2)
        MainWindow.label_8.setStyleSheet("QLabel{\n"
"    background: rgb(255,255,255);\n"
"    border: 0px solid rgb(0,0,128);\n"
"}")
        MainWindow.label_8.setObjectName("label_8")
        MainWindow.gridLayout_4.addWidget(MainWindow.label_8, 3, 0, 1, 1)
        MainWindow.label_6 = QtWidgets.QLabel(MainWindow.widget_2)
        MainWindow.label_6.setStyleSheet("QLabel{\n"
"    background: rgb(255,255,255);\n"
"    border: 0px solid rgb(0,0,128);\n"
"}")
        MainWindow.label_6.setObjectName("label_6")
        MainWindow.gridLayout_4.addWidget(MainWindow.label_6, 0, 0, 1, 1)
        MainWindow.label_7 = QtWidgets.QLabel(MainWindow.widget_2)
        MainWindow.label_7.setStyleSheet("QLabel{\n"
"    background: rgb(255,255,255);\n"
"    border: 0px solid rgb(0,0,128);\n"
"}")
        MainWindow.label_7.setObjectName("label_7")
        MainWindow.gridLayout_4.addWidget(MainWindow.label_7, 6, 0, 1, 1)
        spacerItem11 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        MainWindow.gridLayout_4.addItem(spacerItem11, 8, 0, 1, 1)
        MainWindow.radioButtonDisplayEdgeLabelsFalse = QtWidgets.QRadioButton(MainWindow.widget_2)
        MainWindow.radioButtonDisplayEdgeLabelsFalse.setStyleSheet("QRadioButton{\n"
"    background: rgb(255,255,255);\n"
"    border: 0px solid rgb(0,0,128);\n"
"}")
        MainWindow.radioButtonDisplayEdgeLabelsFalse.setObjectName("radioButtonDisplayEdgeLabelsFalse")
        MainWindow.buttonGroupDOFsPerNode.addButton(MainWindow.radioButtonDisplayEdgeLabelsFalse)
        MainWindow.gridLayout_4.addWidget(MainWindow.radioButtonDisplayEdgeLabelsFalse, 7, 1, 1, 1)
        spacerItem12 = QtWidgets.QSpacerItem(20, 10, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        MainWindow.gridLayout_4.addItem(spacerItem12, 5, 1, 1, 1)
        MainWindow.verticalLayout_5.addWidget(MainWindow.widget_2)
        MainWindow.stackedWidgetRight.addWidget(MainWindow.page_7)
        MainWindow.horizontalLayout.addWidget(MainWindow.stackedWidgetRight)
        MainWindow.verticalLayout_2.addLayout(MainWindow.horizontalLayout)
        MainWindow.horizontalLayout_7 = QtWidgets.QHBoxLayout()
        MainWindow.horizontalLayout_7.setObjectName("horizontalLayout_7")
        spacerItem13 = QtWidgets.QSpacerItem(10, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum)
        MainWindow.horizontalLayout_7.addItem(spacerItem13)
        MainWindow.verticalLayout_17 = QtWidgets.QVBoxLayout()
        MainWindow.verticalLayout_17.setObjectName("verticalLayout_17")
        MainWindow.textBrowser = QtWidgets.QTextBrowser(MainWindow.centralwidget)
        MainWindow.textBrowser.setMaximumSize(QtCore.QSize(16777215, 150))
        MainWindow.textBrowser.setStyleSheet("QTextBrowser{\n"
"border: 2px solid rgb(0,0,128);\n"
"background: rgb(255,255,255);\n"
"}")
        MainWindow.textBrowser.setObjectName("textBrowser")
        MainWindow.verticalLayout_17.addWidget(MainWindow.textBrowser)
        MainWindow.horizontalLayout_7.addLayout(MainWindow.verticalLayout_17)
        MainWindow.verticalLayout_18 = QtWidgets.QVBoxLayout()
        MainWindow.verticalLayout_18.setObjectName("verticalLayout_18")
        MainWindow.hideTextBrowserButton = QtWidgets.QToolButton(MainWindow.centralwidget)
        MainWindow.hideTextBrowserButton.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:pressed {\n"
"     background-color: rgb(185, 211, 220);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgba(185, 211, 220,100);\n"
"}")
        MainWindow.hideTextBrowserButton.setText("")
        icon18 = QtGui.QIcon()
        icon18.addPixmap(QtGui.QPixmap(":/Images/minus.svg"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.hideTextBrowserButton.setIcon(icon18)
        MainWindow.hideTextBrowserButton.setIconSize(QtCore.QSize(24, 24))
        MainWindow.hideTextBrowserButton.setObjectName("hideTextBrowserButton")
        MainWindow.verticalLayout_18.addWidget(MainWindow.hideTextBrowserButton)
        MainWindow.showTextBrowserButton = QtWidgets.QToolButton(MainWindow.centralwidget)
        MainWindow.showTextBrowserButton.setEnabled(True)
        MainWindow.showTextBrowserButton.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:pressed {\n"
"     background-color: rgb(185, 211, 220);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgba(185, 211, 220,100);\n"
"}")
        MainWindow.showTextBrowserButton.setText("")
        icon19 = QtGui.QIcon()
        icon19.addPixmap(QtGui.QPixmap(":/Images/plus.svg"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.showTextBrowserButton.setIcon(icon19)
        MainWindow.showTextBrowserButton.setIconSize(QtCore.QSize(24, 24))
        MainWindow.showTextBrowserButton.setObjectName("showTextBrowserButton")
        MainWindow.verticalLayout_18.addWidget(MainWindow.showTextBrowserButton)
        MainWindow.clearTextBrowserButton = QtWidgets.QToolButton(MainWindow.centralwidget)
        MainWindow.clearTextBrowserButton.setStyleSheet("QToolButton{\n"
"    border: 0px;\n"
"}\n"
"\n"
"QToolButton:pressed {\n"
"     background-color: rgb(185, 211, 220);\n"
"}\n"
"\n"
"QToolButton:hover {\n"
"    background-color: rgba(185, 211, 220,100);\n"
"}\n"
"")
        MainWindow.clearTextBrowserButton.setText("")
        MainWindow.clearTextBrowserButton.setIcon(icon9)
        MainWindow.clearTextBrowserButton.setIconSize(QtCore.QSize(24, 24))
        MainWindow.clearTextBrowserButton.setObjectName("clearTextBrowserButton")
        MainWindow.verticalLayout_18.addWidget(MainWindow.clearTextBrowserButton)
        MainWindow.horizontalLayout_7.addLayout(MainWindow.verticalLayout_18)
        MainWindow.verticalLayout_2.addLayout(MainWindow.horizontalLayout_7)
        MainWindow.setCentralWidget(MainWindow.centralwidget)
        MainWindow.menubar = QtWidgets.QMenuBar(MainWindow)
        MainWindow.menubar.setGeometry(QtCore.QRect(0, 0, 1493, 21))
        MainWindow.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(MainWindow.menubar)
        MainWindow.statusbar = QtWidgets.QStatusBar(MainWindow)
        MainWindow.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(MainWindow.statusbar)

        self.retranslateUi(MainWindow)
        MainWindow.tabWidget.setCurrentIndex(0)
        MainWindow.stackedWidget.setCurrentIndex(1)
        MainWindow.stackedWidgetRight.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        MainWindow.zoomOutButtonSurface.setText(_translate("MainWindow", "Zoom Out"))
        MainWindow.arrowButtonSurface.setText(_translate("MainWindow", "Arrow"))
        MainWindow.panningButtonSurface.setText(_translate("MainWindow", "Panning"))
        MainWindow.zoomInButtonSurface.setText(_translate("MainWindow", "Zoom In"))
        MainWindow.polyButton.setText(_translate("MainWindow", "Draw Polygon"))
        MainWindow.rectButton.setText(_translate("MainWindow", "Draw Rectangle"))
        MainWindow.addPolyHoleButton.setText(_translate("MainWindow", "Draw Polygon Hole"))
        MainWindow.addRectHoleButton.setText(_translate("MainWindow", "Draw Rectangle Hole"))
        MainWindow.mergeButton.setText(_translate("MainWindow", "Merge Overlapping Geometries"))
        MainWindow.deleteButton.setText(_translate("MainWindow", "Delete Item"))
        MainWindow.gridSnapButtonSurface.setText(_translate("MainWindow", "Toggle Gridsnap"))
        MainWindow.gridButtonSurface.setText(_translate("MainWindow", "Toggle Grid"))
        MainWindow.gridSpacingButton.setText(_translate("MainWindow", "Set Grid Spacing"))
        MainWindow.loadGeometryButton.setText(_translate("MainWindow", "Load Geometry"))
        MainWindow.tabWidget.setTabText(MainWindow.tabWidget.indexOf(MainWindow.tab_2), _translate("MainWindow", "Page"))
        MainWindow.arrowButtonBorder.setText(_translate("MainWindow", "Arrow"))
        MainWindow.zoomOutButtonBorder.setText(_translate("MainWindow", "Zoom Out"))
        MainWindow.zoomInButtonBorder.setText(_translate("MainWindow", "Zoom In"))
        MainWindow.panningButtonBorder.setText(_translate("MainWindow", "Panning"))
        MainWindow.setMarkerButton.setText(_translate("MainWindow", "Set Marker"))
        MainWindow.splitEdgeButton.setText(_translate("MainWindow", "Split Edge"))
        MainWindow.gridSnapButtonBorder.setText(_translate("MainWindow", "Toggle Gridsnap"))
        MainWindow.gridButtonBorder.setText(_translate("MainWindow", "Toggle Grid"))
        MainWindow.tabWidget.setTabText(MainWindow.tabWidget.indexOf(MainWindow.tab), _translate("MainWindow", "Tab 1"))
        MainWindow.saveGeomButton.setText(_translate("MainWindow", "Save Geometry"))
        MainWindow.tabWidget.setTabText(MainWindow.tabWidget.indexOf(MainWindow.tab_4), _translate("MainWindow", "Page"))
        MainWindow.saveMeshButton.setText(_translate("MainWindow", "Save Mesh"))
        MainWindow.saveArraysButton.setText(_translate("MainWindow", "Save Mesh Arrays"))
        MainWindow.saveArraysMatlabButton.setText(_translate("MainWindow", "Save Mesh Arrays to Matlab"))
        MainWindow.tabWidget.setTabText(MainWindow.tabWidget.indexOf(MainWindow.tab_3), _translate("MainWindow", "Page"))
        MainWindow.label.setText(_translate("MainWindow", "x: "))
        MainWindow.labelX.setText(_translate("MainWindow", "TextLabel"))
        MainWindow.label_3.setText(_translate("MainWindow", "y:"))
        MainWindow.labelY.setText(_translate("MainWindow", "TextLabel"))
        MainWindow.label_4.setText(_translate("MainWindow", "Element Size:"))
        MainWindow.radioButtonQuadrangle.setText(_translate("MainWindow", "Quadrangle"))
        MainWindow.label_5.setText(_translate("MainWindow", "DOFs per Node:"))
        MainWindow.radioButtonTriangle.setText(_translate("MainWindow", "Triangle"))
        MainWindow.label_2.setText(_translate("MainWindow", "Element Type:"))
        MainWindow.refreshMeshButton.setText(_translate("MainWindow", "Refresh Mesh"))
        MainWindow.radioButtonDisplayPointLabelsTrue.setText(_translate("MainWindow", "True"))
        MainWindow.radioButtonDisplayEdgeLabelsTrue.setText(_translate("MainWindow", "True"))
        MainWindow.radioButtonDisplayPointLabelsFalse.setText(_translate("MainWindow", "False"))
        MainWindow.radioButtonDrawPointFalse.setText(_translate("MainWindow", "False"))
        MainWindow.radioButtonDrawPointTrue.setText(_translate("MainWindow", "True"))
        MainWindow.label_8.setText(_translate("MainWindow", "Display Point Labels:"))
        MainWindow.label_6.setText(_translate("MainWindow", "Draw Points:"))
        MainWindow.label_7.setText(_translate("MainWindow", "Display Edge Labels"))
        MainWindow.radioButtonDisplayEdgeLabelsFalse.setText(_translate("MainWindow", "False"))
        MainWindow.textBrowser.setHtml(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))

#import resources_rc
