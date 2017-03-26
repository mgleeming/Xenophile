# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui/mzML_viewer.ui'
#
# Created: Mon Feb  6 19:55:19 2017
#      by: PyQt4 UI code generator 4.10.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName(_fromUtf8("Dialog"))
        Dialog.resize(1305, 718)
        self.gridLayout_2 = QtGui.QGridLayout(Dialog)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.gridLayout = QtGui.QGridLayout()
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.load_button = QtGui.QPushButton(Dialog)
        self.load_button.setObjectName(_fromUtf8("load_button"))
        self.verticalLayout.addWidget(self.load_button)
        self.list = QtGui.QTreeWidget(Dialog)
        self.list.setObjectName(_fromUtf8("list"))
        self.verticalLayout.addWidget(self.list)
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        self.go_to_rt_field = QtGui.QLineEdit(Dialog)
        self.go_to_rt_field.setObjectName(_fromUtf8("go_to_rt_field"))
        self.horizontalLayout_3.addWidget(self.go_to_rt_field)
        self.go_to_rt_button = QtGui.QPushButton(Dialog)
        self.go_to_rt_button.setObjectName(_fromUtf8("go_to_rt_button"))
        self.horizontalLayout_3.addWidget(self.go_to_rt_button)
        self.verticalLayout.addLayout(self.horizontalLayout_3)
        self.plot_chromatogram = QtGui.QPushButton(Dialog)
        self.plot_chromatogram.setObjectName(_fromUtf8("plot_chromatogram"))
        self.verticalLayout.addWidget(self.plot_chromatogram)
        self.closeHighlighted = QtGui.QPushButton(Dialog)
        self.closeHighlighted.setObjectName(_fromUtf8("closeHighlighted"))
        self.verticalLayout.addWidget(self.closeHighlighted)
        self.gridLayout.addLayout(self.verticalLayout, 0, 0, 1, 1)
        self.traces = GraphicsLayoutWidget(Dialog)
        self.traces.setObjectName(_fromUtf8("traces"))
        self.gridLayout.addWidget(self.traces, 0, 1, 1, 1)
        self.gridLayout.setColumnStretch(1, 4)
        self.gridLayout_2.addLayout(self.gridLayout, 0, 0, 1, 1)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "Dialog", None))
        self.load_button.setText(_translate("Dialog", "Load Data FIles", None))
        self.list.headerItem().setText(0, _translate("Dialog", "File", None))
        self.go_to_rt_button.setText(_translate("Dialog", "Go to RT", None))
        self.plot_chromatogram.setText(_translate("Dialog", "Plot Chromatogram", None))
        self.closeHighlighted.setText(_translate("Dialog", "Close Highlighted Datafile", None))

from pyqtgraph import GraphicsLayoutWidget

if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Dialog = QtGui.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())

