# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'EIC_dialog.ui'
#
# Created: Tue Aug 30 15:22:29 2016
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
        Dialog.resize(733, 362)
        self.gridLayout_3 = QtGui.QGridLayout(Dialog)
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.label_3 = QtGui.QLabel(Dialog)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.verticalLayout_2.addWidget(self.label_3)
        self.EIC_dialog_select_files = QtGui.QListWidget(Dialog)
        self.EIC_dialog_select_files.setObjectName(_fromUtf8("EIC_dialog_select_files"))
        self.verticalLayout_2.addWidget(self.EIC_dialog_select_files)
        self.gridLayout_3.addLayout(self.verticalLayout_2, 0, 0, 1, 1)
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        spacerItem = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.gridLayout_2 = QtGui.QGridLayout()
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.EIC_dialog_targets = QtGui.QLineEdit(Dialog)
        self.EIC_dialog_targets.setObjectName(_fromUtf8("EIC_dialog_targets"))
        self.gridLayout_2.addWidget(self.EIC_dialog_targets, 1, 0, 1, 1)
        self.label_2 = QtGui.QLabel(Dialog)
        self.label_2.setMaximumSize(QtCore.QSize(16777215, 100))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout_2.addWidget(self.label_2, 0, 0, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout_2)
        spacerItem1 = QtGui.QSpacerItem(20, 18, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem1)
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.label = QtGui.QLabel(Dialog)
        self.label.setObjectName(_fromUtf8("label"))
        self.horizontalLayout_2.addWidget(self.label)
        self.EIC_dialog_width = QtGui.QLineEdit(Dialog)
        self.EIC_dialog_width.setObjectName(_fromUtf8("EIC_dialog_width"))
        self.horizontalLayout_2.addWidget(self.EIC_dialog_width)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        spacerItem2 = QtGui.QSpacerItem(20, 18, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem2)
        self.gridLayout = QtGui.QGridLayout()
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.label_4 = QtGui.QLabel(Dialog)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.gridLayout.addWidget(self.label_4, 0, 0, 1, 1)
        spacerItem3 = QtGui.QSpacerItem(58, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem3, 0, 1, 1, 1)
        spacerItem4 = QtGui.QSpacerItem(58, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem4, 0, 2, 1, 1)
        spacerItem5 = QtGui.QSpacerItem(68, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem5, 1, 0, 1, 1)
        self.msLevel_MS1 = QtGui.QCheckBox(Dialog)
        self.msLevel_MS1.setChecked(True)
        self.msLevel_MS1.setObjectName(_fromUtf8("msLevel_MS1"))
        self.gridLayout.addWidget(self.msLevel_MS1, 1, 1, 1, 1)
        spacerItem6 = QtGui.QSpacerItem(68, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem6, 1, 2, 1, 1)
        spacerItem7 = QtGui.QSpacerItem(68, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem7, 2, 0, 1, 1)
        self.msLevel_MS2 = QtGui.QCheckBox(Dialog)
        self.msLevel_MS2.setChecked(True)
        self.msLevel_MS2.setObjectName(_fromUtf8("msLevel_MS2"))
        self.gridLayout.addWidget(self.msLevel_MS2, 2, 1, 1, 1)
        spacerItem8 = QtGui.QSpacerItem(68, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem8, 2, 2, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout)
        spacerItem9 = QtGui.QSpacerItem(20, 18, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem9)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.EIC_dialog_OK = QtGui.QPushButton(Dialog)
        self.EIC_dialog_OK.setObjectName(_fromUtf8("EIC_dialog_OK"))
        self.horizontalLayout.addWidget(self.EIC_dialog_OK)
        self.EIC_dialog_cancel = QtGui.QPushButton(Dialog)
        self.EIC_dialog_cancel.setObjectName(_fromUtf8("EIC_dialog_cancel"))
        self.horizontalLayout.addWidget(self.EIC_dialog_cancel)
        self.verticalLayout.addLayout(self.horizontalLayout)
        spacerItem10 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem10)
        self.gridLayout_3.addLayout(self.verticalLayout, 0, 1, 1, 1)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "Dialog", None))
        self.label_3.setText(_translate("Dialog", "Select FIles", None))
        self.label_2.setText(_translate("Dialog", "Targets", None))
        self.label.setText(_translate("Dialog", "EIC width (+/-- m/z)", None))
        self.label_4.setText(_translate("Dialog", "MS Level", None))
        self.msLevel_MS1.setText(_translate("Dialog", "MS1", None))
        self.msLevel_MS2.setText(_translate("Dialog", "MS2", None))
        self.EIC_dialog_OK.setText(_translate("Dialog", "OK", None))
        self.EIC_dialog_cancel.setText(_translate("Dialog", "Cancel", None))


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Dialog = QtGui.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())

