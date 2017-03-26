# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'NTPS_rangefinder.ui'
#
# Created: Sun Aug 21 14:33:51 2016
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
        Dialog.resize(943, 778)
        self.gridLayout_3 = QtGui.QGridLayout(Dialog)
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.gridLayout = QtGui.QGridLayout()
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.label = QtGui.QLabel(Dialog)
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.label_4 = QtGui.QLabel(Dialog)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.gridLayout.addWidget(self.label_4, 0, 1, 1, 1)
        self.label_3 = QtGui.QLabel(Dialog)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.gridLayout.addWidget(self.label_3, 0, 2, 1, 1)
        self.HT_input_button = QtGui.QPushButton(Dialog)
        self.HT_input_button.setObjectName(_fromUtf8("HT_input_button"))
        self.gridLayout.addWidget(self.HT_input_button, 1, 0, 1, 1)
        self.HT_input_field = QtGui.QLineEdit(Dialog)
        self.HT_input_field.setObjectName(_fromUtf8("HT_input_field"))
        self.gridLayout.addWidget(self.HT_input_field, 1, 1, 1, 1)
        self.HT_show_checkbox = QtGui.QCheckBox(Dialog)
        self.HT_show_checkbox.setText(_fromUtf8(""))
        self.HT_show_checkbox.setObjectName(_fromUtf8("HT_show_checkbox"))
        self.gridLayout.addWidget(self.HT_show_checkbox, 1, 2, 1, 1)
        self.Mascot_input_button = QtGui.QPushButton(Dialog)
        self.Mascot_input_button.setObjectName(_fromUtf8("Mascot_input_button"))
        self.gridLayout.addWidget(self.Mascot_input_button, 2, 0, 1, 1)
        self.Mascot_input_field = QtGui.QLineEdit(Dialog)
        self.Mascot_input_field.setObjectName(_fromUtf8("Mascot_input_field"))
        self.gridLayout.addWidget(self.Mascot_input_field, 2, 1, 1, 1)
        self.Mascot_show_checkbox = QtGui.QCheckBox(Dialog)
        self.Mascot_show_checkbox.setText(_fromUtf8(""))
        self.Mascot_show_checkbox.setObjectName(_fromUtf8("Mascot_show_checkbox"))
        self.gridLayout.addWidget(self.Mascot_show_checkbox, 2, 2, 1, 1)
        self.horizontalLayout_2.addLayout(self.gridLayout)
        self.gridLayout_2 = QtGui.QGridLayout()
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.label_8 = QtGui.QLabel(Dialog)
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.gridLayout_2.addWidget(self.label_8, 3, 0, 1, 1)
        self.label_5 = QtGui.QLabel(Dialog)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.gridLayout_2.addWidget(self.label_5, 1, 0, 1, 1)
        self.RT_range = QtGui.QLineEdit(Dialog)
        self.RT_range.setObjectName(_fromUtf8("RT_range"))
        self.gridLayout_2.addWidget(self.RT_range, 1, 1, 1, 1)
        self.label_6 = QtGui.QLabel(Dialog)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.gridLayout_2.addWidget(self.label_6, 2, 0, 1, 1)
        self.mz_range = QtGui.QLineEdit(Dialog)
        self.mz_range.setObjectName(_fromUtf8("mz_range"))
        self.gridLayout_2.addWidget(self.mz_range, 2, 1, 1, 1)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.HT_timeunit_mins = QtGui.QRadioButton(Dialog)
        self.HT_timeunit_mins.setObjectName(_fromUtf8("HT_timeunit_mins"))
        self.horizontalLayout.addWidget(self.HT_timeunit_mins)
        self.HT_timeunit_secs = QtGui.QRadioButton(Dialog)
        self.HT_timeunit_secs.setObjectName(_fromUtf8("HT_timeunit_secs"))
        self.horizontalLayout.addWidget(self.HT_timeunit_secs)
        self.gridLayout_2.addLayout(self.horizontalLayout, 0, 1, 1, 1)
        self.label_7 = QtGui.QLabel(Dialog)
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.gridLayout_2.addWidget(self.label_7, 0, 0, 1, 1)
        self.HT_threshold = QtGui.QLineEdit(Dialog)
        self.HT_threshold.setObjectName(_fromUtf8("HT_threshold"))
        self.gridLayout_2.addWidget(self.HT_threshold, 3, 1, 1, 1)
        self.horizontalLayout_2.addLayout(self.gridLayout_2)
        self.RF_done = QtGui.QPushButton(Dialog)
        self.RF_done.setObjectName(_fromUtf8("RF_done"))
        self.horizontalLayout_2.addWidget(self.RF_done)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.RF_plot = PlotWidget(Dialog)
        self.RF_plot.setObjectName(_fromUtf8("RF_plot"))
        self.verticalLayout.addWidget(self.RF_plot)
        self.gridLayout_3.addLayout(self.verticalLayout, 0, 0, 1, 1)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "Dialog", None))
        self.label.setText(_translate("Dialog", "Select File", None))
        self.label_4.setText(_translate("Dialog", "File Name", None))
        self.label_3.setText(_translate("Dialog", "Show", None))
        self.HT_input_button.setText(_translate("Dialog", "HiTIME input", None))
        self.Mascot_input_button.setText(_translate("Dialog", "Mascot Input", None))
        self.label_8.setText(_translate("Dialog", "HT threshold", None))
        self.label_5.setText(_translate("Dialog", "RT range", None))
        self.label_6.setText(_translate("Dialog", "m/z range", None))
        self.HT_timeunit_mins.setText(_translate("Dialog", "Minutes", None))
        self.HT_timeunit_secs.setText(_translate("Dialog", "Seconds", None))
        self.label_7.setText(_translate("Dialog", "HT Time Units", None))
        self.RF_done.setText(_translate("Dialog", "Done", None))

from pyqtgraph import PlotWidget

if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Dialog = QtGui.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())

