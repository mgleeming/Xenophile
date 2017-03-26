# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'launch.ui'
#
# Created: Tue Mar 22 17:36:23 2016
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
        Dialog.resize(747, 347)
        self.NTPS_button = QtGui.QPushButton(Dialog)
        self.NTPS_button.setGeometry(QtCore.QRect(30, 40, 261, 44))
        self.NTPS_button.setObjectName(_fromUtf8("NTPS_button"))
        self.label = QtGui.QLabel(Dialog)
        self.label.setGeometry(QtCore.QRect(320, 40, 361, 41))
        self.label.setObjectName(_fromUtf8("label"))
        self.targeted_search_button = QtGui.QPushButton(Dialog)
        self.targeted_search_button.setGeometry(QtCore.QRect(30, 110, 261, 44))
        self.targeted_search_button.setObjectName(_fromUtf8("targeted_search_button"))
        self.label_2 = QtGui.QLabel(Dialog)
        self.label_2.setGeometry(QtCore.QRect(320, 110, 361, 41))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.HT_search_button = QtGui.QPushButton(Dialog)
        self.HT_search_button.setGeometry(QtCore.QRect(30, 180, 261, 44))
        self.HT_search_button.setObjectName(_fromUtf8("HT_search_button"))
        self.label_3 = QtGui.QLabel(Dialog)
        self.label_3.setGeometry(QtCore.QRect(320, 180, 361, 41))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.mzML_file_viewer_button = QtGui.QPushButton(Dialog)
        self.mzML_file_viewer_button.setGeometry(QtCore.QRect(30, 250, 261, 44))
        self.mzML_file_viewer_button.setObjectName(_fromUtf8("mzML_file_viewer_button"))
        self.label_4 = QtGui.QLabel(Dialog)
        self.label_4.setGeometry(QtCore.QRect(320, 250, 361, 41))
        self.label_4.setObjectName(_fromUtf8("label_4"))

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "Dialog", None))
        self.NTPS_button.setText(_translate("Dialog", "Non-Target Reactive\n"
" Metabolite Identification", None))
        self.label.setText(_translate("Dialog", "Conduct a non-targeted search of LC-MS/MS data to\n"
" determine the identity of reactive metabolites", None))
        self.targeted_search_button.setText(_translate("Dialog", "Targeted CRM protein \n"
"adduct identification", None))
        self.label_2.setText(_translate("Dialog", "Run a targeted search for protein\n"
" adducts of a specified CRM", None))
        self.HT_search_button.setText(_translate("Dialog", "HiTIME search", None))
        self.label_3.setText(_translate("Dialog", "Run a HiTIME search", None))
        self.mzML_file_viewer_button.setText(_translate("Dialog", "mzML file vieser", None))
        self.label_4.setText(_translate("Dialog", "Interactive viewer for mzML files", None))


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Dialog = QtGui.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())

